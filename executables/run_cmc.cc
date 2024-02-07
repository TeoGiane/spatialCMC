// Include STL headers
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <random>

// Include openmp headers
#include <omp.h>

// Include bayesmix headers
#include "src/includes.h"

// Include argparse: input files
#include "lib/argparse/argparse.h"

// Local Inclusions
#include "hierarchies/poisson_gamma_hierarchy.h"
#include "mixing/spp_mixing.h"
#include "shard.h"
#include "shard_merger.h"
#include "shard_partition.h"
#include "spatialcmc_utils.h"

// Check ArgParse input arguments
void check_args(const argparse::ArgumentParser &args) {
  try {
    // Check if input files exist
    spatialcmc::check_file_is_readable(args.get<std::string>("--data-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--adj-matrix-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--shard-assignment-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--algo-params-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--hier-prior-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--mix-prior-file"));
    if (args["--chain-file"] != EMPTYSTR) {
      bayesmix::check_file_is_writeable(args.get<std::string>("--chain-file"));
    }
  }
  catch(const std::runtime_error & err) {
    std::cerr << err.what() << std::endl;
    std::exit(1);
  }
  return;
}

/* Main Function - SpatialCMC Sampler */

int main(int argc, char const *argv[]) {

	// Set up argparse arguments in input
  argparse::ArgumentParser args("spatialCMC::run_cmc");
  args.add_argument("--data-file")
    .required()
    .help("Path to a .csv file containing the observations (one per row)");
  args.add_argument("--adj-matrix-file")
  .required()
  .help("Path to a .csv file containing the adjacency matrix");
  args.add_argument("--shard-assignment-file")
    .required()
    .help("Path to a .csv file containing the shard assignment for each observation (one per row)");
  args.add_argument("--algo-params-file")
    .required()
    .help("Path to .asciipb file with the parameters of the algorithm");
  args.add_argument("--hier-type")
    .required()
    .help("Enum string of the hierarchy");
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
  }
  catch (const std::runtime_error &err) {
    std::cerr << err.what() << std::endl;
    std::cout << std::endl;
		std::cout << args << std::endl;
    std::exit(1);
  }

	// Check parsed arguments
  check_args(args);

  // Set number of threads - 75% of available threads
  int num_threads;
  #pragma omp parallel
  {
    num_threads =  omp_get_num_threads();
  }
  num_threads *= 0.75;

  // Startup message
  std::cout << std::endl;
  std::cout << "Running run_cmc.cc" << std::endl;

  // Read data files
  auto data = bayesmix::read_eigen_matrix(args.get<std::string>("--data-file"));
  auto adj_matrix = bayesmix::read_eigen_matrix(args.get<std::string>("--adj-matrix-file"));
  auto shard_allocation = spatialcmc::read_shard_allocation_file(args.get<std::string>("--shard-assignment-file"));

  // Create hierarchy object
  std::shared_ptr<AbstractHierarchy> hierarchy;
  if(args.get<std::string>("--hier-type") == "PoissonGamma") {
    hierarchy = std::make_shared<PoissonGammaHierarchy>();
  } else {
    auto & factory_hier = HierarchyFactory::Instance();
    hierarchy = factory_hier.create_object(args.get<std::string>("--hier-type"));
  }
  
  // Create mixing object
  auto mixing = std::make_shared<sPPMixing>();

  // Split data in shards
  unsigned int num_shards = std::set<int>(shard_allocation.cbegin(), shard_allocation.cend()).size();
  std::vector<std::deque<int>> global_numbering(num_shards);
  std::vector<Eigen::MatrixXd> data_in_shards(num_shards);
  std::vector<Shard> shards(num_shards);
  for (size_t i = 0; i < data.size(); i++) {
    global_numbering[shard_allocation[i]].push_back(i);
  }

  // Generating shards in parallel
  #pragma omp parallel for num_threads(num_threads)
	for (size_t i = 0; i < num_shards; i++) {
    shards[i].set_data(data(global_numbering[i], 0));
    shards[i].set_adjacency_matrix(adj_matrix(global_numbering[i], global_numbering[i]));
    shards[i].set_sampler_parameters(args.get<std::string>("--algo-params-file"));
    shards[i].set_hierarchy(hierarchy->clone(), args.get<std::string>("--hier-prior-file"));
    shards[i].set_mixing(mixing->clone(), args.get<std::string>("--mix-prior-file"));
    if(i != 0) 
      shards[i].set_verbose(false);
  }

  // Check - OK
  // std::cout << "N° of shards: " << shards.size() << std::endl;
  // for (auto && elem : shards) { elem.print(); std::cout << std::endl; }

  // Run MCMC samplers in each shard - in parallel --> ConsensusMCSampler::parallel_sample();
  std::vector<spatialcmc::MCMCChain> chains(num_shards);
  #pragma omp parallel for num_threads(num_threads)
  for (size_t i = 0; i < num_shards; i++) {
    chains[i] = shards[i].sample();
  }

  // Prepare buffers for merging MCMC chains in shards
  ShardMerger shard_merger(shards, chains, global_numbering, adj_matrix);
  std::vector<bayesmix::AlgorithmState> state_vect(shard_merger.get_num_iter());
  progresscpp::ProgressBar* bar = new progresscpp::ProgressBar(shard_merger.get_num_iter(), 60);

  // Merging shards in parallel
  std::cout << "Merging MCMC chains..." << std::endl;
  #pragma omp parallel for num_threads(num_threads)
  for (size_t i = 0; i < state_vect.size(); i++) {
    state_vect[i].CopyFrom(shard_merger.merge(i));
		++(*bar);
    #pragma omp critical
		bar->display();
  }
	bar->done();
	std::cout << "Done" << std::endl;

  // Serialization of the final chain
  spatialcmc::MCMCChain merged_chain;
  * merged_chain.mutable_state() = { state_vect.begin(), state_vect.end() };
  std::cout << "Successfully serialized final MCMC chain" << std::endl;

  // Write to file in serialized format
  if (args["--chain-file"] != EMPTYSTR) {
    std::ofstream outfile(args.get<std::string>("--chain-file"));
    outfile << merged_chain.SerializeAsString();
    std::cout << "Successfully wrote serialized MCMC chain to "
              << args.get<std::string>("--chain-file") << std::endl;
  }

  // Final message
  std::cout << "End of run_cmc.cc" << std::endl;
  std::cout << std::endl;

  // Memory management
  google::protobuf::ShutdownProtobufLibrary();
  delete bar;

  // Terminate program with no error
  return 0;
}
