#include <fstream>
#include <iostream>
#include <math.h>

#include "lib/argparse/argparse.h"
#include "src/includes.h"

#include "mcmc_chain.pb.h"
#include "hierarchies/poisson_gamma_hierarchy.h"
#include "hierarchies/empty_hierarchy.h"
#include "mixing/spp_mixing.h"
#include "cmc/spatialcmc_utils.h"

#define EMPTYSTR std::string("\"\"")

bool check_args(const argparse::ArgumentParser &args) {
	spatialcmc::check_file_is_readable(args.get<std::string>("--data-file"));
  if (args.get<std::string>("--adj-matrix-file") != EMPTYSTR) {
	  spatialcmc::check_file_is_readable(args.get<std::string>("--adj-matrix-file"));
  }
	spatialcmc::check_file_is_readable(args.get<std::string>("--algo-params-file"));
	spatialcmc::check_file_is_readable(args.get<std::string>("--hier-prior-file"));
	spatialcmc::check_file_is_readable(args.get<std::string>("--mix-prior-file"));
  if (args.get<std::string>("--chain-file") != EMPTYSTR) {
    bayesmix::check_file_is_writeable(args.get<std::string>("--chain-file"));
  }
  return true;
}

int main(int argc, char *argv[]) {

	// Set up argparse arguments in input
  argparse::ArgumentParser args("spatialCMC::run_mcmc");
	args.add_argument("--data-file")
		.required()
		.help("Path to a .csv file containing the observations (one per row)");
  args.add_argument("--adj-matrix-file")
    .required()
    .default_value(EMPTYSTR)
    .help("Path to a .csv file with the adjacency matrix used in the mixing");
  args.add_argument("--algo-params-file")
    .required()
    .help("Path to .asciipb file with the parameters of the algorithm");
  args.add_argument("--hier-type")
    .required()
    .help("Enum string of the hierarchy");
  args.add_argument("--hier-prior-file")
		.required()
		.help("Path to .asciipb file with the parameters of the hierarchy");
  args.add_argument("--mix-type")
    .required()
    .default_value(std::string("sPP"))
    .help("Enum string of the mixing");
  args.add_argument("--mix-prior-file")
		.required()
		.help("Path to .asciipb file with the parameters of the mixing");
  args.add_argument("--chain-file")
    .required()
    .default_value(EMPTYSTR)
    .help("Path to a .recordio file where to store the MCMC chain, serialized via Google Protocol Buffer");

	// Parse arguments
  try {
    args.parse_args(argc, argv);
  }
  catch (const std::runtime_error &err) {
    std::cerr << err.what() << std::endl;
    std::cout << std::endl;
		std::cout << args << std::endl;
    std::exit(1);
  }

	// Check parsed arguments
  check_args(args);

	// Startup message
  std::cout << std::endl;
	std::cout << "Running run_mcmc.cc" << std::endl;

  // Read data files
  auto data = bayesmix::read_eigen_matrix(args.get<std::string>("--data-file"));
  Eigen::MatrixXd adj_matrix;
  if (args.get<std::string>("--adj-matrix-file") != EMPTYSTR){
    adj_matrix = bayesmix::read_eigen_matrix(args.get<std::string>("--adj-matrix-file"));
  }

  // Create hierarchy object
  std::shared_ptr<AbstractHierarchy> hierarchy;
  if(args.get<std::string>("--hier-type") == "PoissonGamma") {
    hierarchy = std::make_shared<PoissonGammaHierarchy>();
  } else if (args.get<std::string>("--hier-type") == "Empty") {
    hierarchy = std::make_shared<EmptyHierarchy>();
  } else {
    auto & factory_hier = HierarchyFactory::Instance();
    hierarchy = factory_hier.create_object(args.get<std::string>("--hier-type"));
  }
  bayesmix::read_proto_from_file(args.get<std::string>("--hier-prior-file"), hierarchy->get_mutable_prior());

  // Create mixing object
  std::shared_ptr<AbstractMixing> mixing;
  if(args.get<std::string>("--mix-type") == "sPP") {
    mixing = std::make_shared<sPPMixing>();
  } else {
    auto & factory_mixing = MixingFactory::Instance();
    mixing = factory_mixing.create_object(args.get<std::string>("--mix-type"));
  }
  bayesmix::read_proto_from_file(args.get<std::string>("--mix-prior-file"), mixing->get_mutable_prior());

  // Read algorithm settings proto
  bayesmix::AlgorithmParams algo_proto;
  bayesmix::read_proto_from_file(args.get<std::string>("--algo-params-file"), &algo_proto);

  // Create algorithm object and set up
  auto &factory_algo = AlgorithmFactory::Instance();
  auto algo = factory_algo.create_object(algo_proto.algo_id());
  algo->read_params_from_proto(algo_proto);

  // Proper algorithm initialization
  algo->set_mixing(mixing);
  algo->set_data(data);
  algo->set_hierarchy(hierarchy);

  // Set covariates for mixing
  if (mixing->is_dependent()){
    algo->set_mix_covariates(adj_matrix);
  }

  // Define proper collector
  BaseCollector *coll = new MemoryCollector();

  // Run algorithm
  algo->run(coll);

  // Save serialized MCMC chain to file
  if (args["--chain-file"] != EMPTYSTR) {    
    
    // Create container for the MCMC chain
    spatialcmc::MCMCChain chain;

    // Fill the chain container
    bayesmix::AlgorithmState state;
    for (int i = 0; i < coll->get_size(); i++) {
      coll->get_next_state(&state);
      chain.add_state()->CopyFrom(state);
    }

    // Write to file in serialized format
    std::ofstream outfile(args.get<std::string>("--chain-file"));
    outfile << chain.SerializeAsString();
		std::cout << "Successfully wrote serialized MCMC chain to "
              << args.get<std::string>("--chain-file") << std::endl;
  }	

	// Final message
  std::cout << "End of run_mcmc.cc" << std::endl;
  std::cout << std::endl;

	// Memory management
  delete coll;
  google::protobuf::ShutdownProtobufLibrary();

	// Terminate program with no error
  return 0;
}
