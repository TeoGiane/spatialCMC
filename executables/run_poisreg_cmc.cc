// STL
#include <math.h>
#include <fstream>
#include <iostream>

// OpenMP
#include <omp.h>
#include <stan/math.hpp>
#include <Eigen/Dense>

// argparse
#include "lib/argparse/argparse.h"

// bayesmix
#include "src/includes.h"

// spatialcmc
#include "algorithm/poisson_regression_sampler.h"
#include "cmc/poisson_regression_local_cluster_merger.h"
#include "cmc/regression_coefficients_merger.h"
#include "cmc/spatialcmc_utils.h"
#include "hierarchies/poisson_regression_hierarchy.h"
#include "mixing/spp_mixing.h"
#include "poisson_regression_algorithm_state.pb.h"

#define EMPTYSTR std::string("\"\"")

// Check ArgParse input arguments
void check_args(const argparse::ArgumentParser &args) {
  try {
    // Check if input files exist
    spatialcmc::check_file_is_readable(args.get<std::string>("--data-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--hier-cov-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--adj-matrix-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--shard-assignment-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--algo-params-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--hier-prior-file"));
    spatialcmc::check_file_is_readable(args.get<std::string>("--mix-prior-file"));
    if (args.get<std::string>("--chain-file") != EMPTYSTR) {
      bayesmix::check_file_is_writeable(args.get<std::string>("--chain-file"));
    }
  }
  catch(const std::runtime_error & err) {
    std::cerr << err.what() << std::endl;
    std::exit(1);
  }
  return;
};

void run_mcmc(std::shared_ptr<PoissonRegressionSampler> algo,
              std::vector<bayesmix::AlgorithmState> & partition_chain,
              std::vector<bayesmix::Vector> & reg_coeffs_chain) {
  // Run the MCMC algorithm
  BaseCollector * coll = new MemoryCollector();
  algo->run(coll);
  // Create chain buffers
  partition_chain.resize(coll->get_size());
  reg_coeffs_chain.resize(coll->get_size());
  // Fill the chains containers
  spatialcmc::PoissonRegAlgorithmState curr_state;
  for (int i = 0; i < coll->get_size(); i++) {
    coll->get_next_state(&curr_state);
    partition_chain[i].CopyFrom(curr_state.partition());
    reg_coeffs_chain[i].CopyFrom(curr_state.regression_coefficients());
  }
  // Deallocate collector
  delete coll;
};

/* Main Function - SpatialCMC Sampler */

int main(int argc, char const *argv[]) {

	// Set up argparse arguments in input
  argparse::ArgumentParser args("spatialCMC::run_poisson_regression_cmc");
  args.add_argument("--data-file")
    .required()
    .help("Path to a .csv file containing the observations (one per row)");
  args.add_argument("--hier-cov-file")
    .required()
    .help("Path to a .csv file with the covariates used in the hierarchy");
  args.add_argument("--adj-matrix-file")
  .required()
  .help("Path to a .csv file containing the adjacency matrix");
  args.add_argument("--shard-assignment-file")
    .required()
    .help("Path to a .csv file containing the shard assignment for each observation (one per row)");
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
  std::cout << "Running run_poisson_regression_cmc.cc" << std::endl;

  // Read data files
  auto data = bayesmix::read_eigen_matrix(args.get<std::string>("--data-file"));
  auto cov_matrix = bayesmix::read_eigen_matrix(args.get<std::string>("--hier-cov-file"));
  auto adj_matrix = bayesmix::read_eigen_matrix(args.get<std::string>("--adj-matrix-file"));
  auto shard_allocation = spatialcmc::read_shard_allocation_file(args.get<std::string>("--shard-assignment-file"));

  // Split data in shards
  unsigned int num_shards = std::set<int>(shard_allocation.cbegin(), shard_allocation.cend()).size();
  std::vector<std::deque<int>> global_numbering(num_shards);
  for (size_t i = 0; i < data.size(); i++) {
    global_numbering[shard_allocation[i]].push_back(i);
  }

  // Read algorithm parameters
  bayesmix::AlgorithmParams algo_proto;
  bayesmix::read_proto_from_file(args.get<std::string>("--algo-params-file"), &algo_proto);

  // Generating shards
  std::vector<std::shared_ptr<PoissonRegressionSampler>> shards(num_shards);
  std::cout << "Generating shards ... ";
  #pragma omp parallel for num_threads(num_threads)
	for (size_t i = 0; i < num_shards; i++) {
    shards[i] = std::make_shared<PoissonRegressionSampler>();
    shards[i]->set_data(data(global_numbering[i], 0));
    shards[i]->set_hier_covariates(cov_matrix(global_numbering[i], Eigen::all));
    shards[i]->set_adjacency_matrix(adj_matrix(global_numbering[i], global_numbering[i]));
    shards[i]->set_hierarchy_prior(args.get<std::string>("--hier-prior-file"));
    shards[i]->set_mixing_prior(args.get<std::string>("--mix-prior-file"));
    shards[i]->read_params_from_proto(algo_proto);
    if(i != 0) 
      shards[i]->set_verbose(false);
  }
  std::cout << "Done" << std::endl;

  // Run MCMC samplers in each shard - in parallel
  std::vector<std::vector<bayesmix::AlgorithmState>> sharded_partitions(num_shards);
  std::vector<std::vector<bayesmix::Vector>> sharded_regression_coefficients(num_shards);
  #pragma omp parallel for num_threads(num_threads)
  for (size_t i = 0; i < num_shards; i++) {
    run_mcmc(shards[i], sharded_partitions[i], sharded_regression_coefficients[i]);
  }

  // Set up hierarchy for partition merger
  std::shared_ptr<AbstractHierarchy> hierarchy = std::make_shared<PoissonRegHierarchy>();
  bayesmix::read_proto_from_file(args.get<std::string>("--hier-prior-file"), hierarchy->get_mutable_prior());

  // Set-up partition merger
  PoissonRegLocalClusterMerger partition_merger(sharded_partitions);
  partition_merger.set_hierarchy(hierarchy->clone());
  partition_merger.set_data(&data);
  partition_merger.set_cov_matrix(&cov_matrix);
  partition_merger.set_adj_matrix(&adj_matrix);
  partition_merger.set_global_numbering(global_numbering);
  partition_merger.set_seed(algo_proto.rng_seed());
  partition_merger.set_hierarchy_prior(args.get<std::string>("--hier-prior-file"));

  // Set-up regression coefficients merger
  RegressionCoefficientsMerger reg_coeffs_merger(sharded_regression_coefficients);

  // Merging shards in parallel
  std::cout << "Merging MCMC chains..." << std::endl;
  std::vector<spatialcmc::PoissonRegAlgorithmState> merged_states(partition_merger.get_num_iter());
  progresscpp::ProgressBar* bar = new progresscpp::ProgressBar(partition_merger.get_num_iter(), 60);
  #pragma omp parallel for num_threads(num_threads)
  for (size_t i = 0; i < merged_states.size(); i++) {
    auto merged_reg_coeffs_proto = reg_coeffs_merger.merge(i);
    merged_states[i].mutable_regression_coefficients()->CopyFrom(merged_reg_coeffs_proto);
    partition_merger.set_reg_coeffs(i, bayesmix::to_eigen(merged_reg_coeffs_proto));
    merged_states[i].mutable_partition()->CopyFrom(partition_merger.merge(i));
		++(*bar);
    #pragma omp critical
		bar->display();
  }
	bar->done();
	std::cout << "Done" << std::endl;

  // Serialization of the final chain
  spatialcmc::PoissonRegMCMC merged_chain;
  * merged_chain.mutable_state() = { merged_states.begin(), merged_states.end() };
  std::cout << "Successfully serialized final MCMC chain" << std::endl;

  // Write to file in serialized format
  if (args["--chain-file"] != EMPTYSTR) {
    std::ofstream outfile(args.get<std::string>("--chain-file"));
    outfile << merged_chain.SerializeAsString();
    std::cout << "Successfully wrote serialized MCMC chain to "
              << args.get<std::string>("--chain-file") << std::endl;
  }

  // Final message
  std::cout << "End of run_poisson_regression_cmc.cc" << std::endl;
  std::cout << std::endl;

  // Memory management
  google::protobuf::ShutdownProtobufLibrary();
  delete bar;

  // Terminate program with no error
  return 0;
};