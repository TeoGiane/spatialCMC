#include <math.h>

#include <fstream>
#include <iostream>

#include "lib/argparse/argparse.h"

#include "algorithm/poisson_regression_sampler.h"
#include "poisson_regression_algorithm_state.pb.h"
// #include "utils/io_utils.h"

#include "src/includes.h"

#define EMPTYSTR std::string("\"\"")

bool check_args(const argparse::ArgumentParser &args) {
  spatialcmc::check_file_is_readable(args.get<std::string>("--data-file"));
  spatialcmc::check_file_is_readable(args.get<std::string>("--hier-cov-file"));
  spatialcmc::check_file_is_readable(args.get<std::string>("--adj-matrix-file"));
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
  argparse::ArgumentParser args("spatialCMC::run_poisson_regression_mcmc");
  args.add_argument("--data-file")
		.required()
    .help("Path to a .csv file containing the observations (one per row)");
  args.add_argument("--hier-cov-file")
    .required()
    .help("Path to a .csv file with the covariates used in the hierarchy");
  args.add_argument("--adj-matrix-file")
    .required()
    .help("Path to a .csv file with the adjacency matrix used in the mixing");
  args.add_argument("--algo-params-file")
      .required()
      .help("Path to .asciipb file with the parameters of the algorithm");
  args.add_argument("--hier-prior-file")
      .required()
      .help("Path to .asciipb file with the parameters of the hierarchy");
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
  } catch (const std::runtime_error &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << args;
    std::exit(1);
  }

  // Check parsed arguments
  check_args(args);


  std::cout << std::endl;
  std::cout << "Running run_poisson_regression_mcmc.cc" << std::endl;

  // Read data files
  auto data = bayesmix::read_eigen_matrix(args.get<std::string>("--data-file"));
  auto cov_matrix = bayesmix::read_eigen_matrix(args.get<std::string>("--hier-cov-file"));
  auto adj_matrix = bayesmix::read_eigen_matrix(args.get<std::string>("--adj-matrix-file"));

  // Read algorithm parameters
  bayesmix::AlgorithmParams algo_proto;
  bayesmix::read_proto_from_file(args.get<std::string>("--algo-params-file"), &algo_proto);

  // Create algorithm object
  auto algo = std::make_shared<PoissonRegressionSampler>();

  // Set-up Poisson Regression Algorithm
  algo->set_data(data);
  algo->set_hier_covariates(cov_matrix);
  algo->set_adjacency_matrix(adj_matrix);
  algo->set_hierarchy_prior(args.get<std::string>("--hier-prior-file"));
  algo->set_mixing_prior(args.get<std::string>("--mix-prior-file"));
  algo->read_params_from_proto(algo_proto);

  // Define proper collector
  BaseCollector *coll = new MemoryCollector();

  // Run algorithm
  algo->run(coll);

  // Save serialized MCMC chain to file
  if (args["--chain-file"] != EMPTYSTR) {    
    
    // Create container for the MCMC chain
    spatialcmc::PoisRegMCMC chain;

    // Fill the chain container
    spatialcmc::PoisRegAlgorithmState state;
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
  std::cout << "End of run_poisreg_mcmc.cc" << std::endl;
  std::cout << std::endl;

	// Memory management
  delete coll;
  google::protobuf::ShutdownProtobufLibrary();

	// Terminate program with no error
  return 0;
}