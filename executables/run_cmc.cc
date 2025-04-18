// Include STL headers
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <random>
#include <typeinfo>

// Include openmp headers
#include <omp.h>

// Include bayesmix headers
#include "src/includes.h"

// Include argparse: input files
#include "lib/argparse/argparse.h"

// Local Inclusions
#include "hierarchies/poisson_gamma_hierarchy.h"
#include "hierarchies/empty_hierarchy.h"
#include "mixing/spp_mixing.h"

#include "mcmc_chain.pb.h"

#include "cmc/local_cluster_merger.h"
// #include "cmc/shard.h"
// #include "cmc/shard_merger.h"
// #include "cmc/shard_partition.h"
#include "cmc/spatialcmc_utils.h"

// Check ArgParse input arguments
void check_args(const argparse::ArgumentParser &args) {
  try {
    // Check if input files exist
    spatialcmc::check_file_is_readable(args.get<std::string>("--data-file"));
    if (args.get<std::string>("--hier-cov-file") != EMPTYSTR) {
      spatialcmc::check_file_is_readable(args.get<std::string>("--hier-cov-file"));
    }
    if (args.get<std::string>("--adj-matrix-file") != EMPTYSTR) {
      spatialcmc::check_file_is_readable(args.get<std::string>("--adj-matrix-file"));
    }
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
}

void run_mcmc(std::shared_ptr<BaseAlgorithm> algo,
              std::vector<bayesmix::AlgorithmState> & partition_chain) {
  // Run the MCMC algorithm
  BaseCollector * coll = new MemoryCollector();
  algo->run(coll);
  // Create chain buffers
  partition_chain.resize(coll->get_size());
  // Fill the chains containers
  bayesmix::AlgorithmState curr_state;
  for (int i = 0; i < coll->get_size(); i++) {
    coll->get_next_state(&curr_state);
    partition_chain[i].CopyFrom(curr_state);
  }
  // Deallocate collector
  delete coll;
};

/* Main Function - SpatialCMC Sampler */

int main(int argc, char const *argv[]) {

	// Set up argparse arguments in input
  argparse::ArgumentParser args("spatialCMC::run_cmc");
  args.add_argument("--data-file")
    .required()
    .help("Path to a .csv file containing the observations (one per row)");
  args.add_argument("--hier-cov-file")
    .required()
    .default_value(EMPTYSTR)
    .help("Path to a .csv file with the covariates used in the hierarchy");
  args.add_argument("--adj-matrix-file")
    .required()
    .default_value(EMPTYSTR)
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
  Eigen::MatrixXd cov_matrix;
  if (args.get<std::string>("--hier-cov-file") != EMPTYSTR){
    cov_matrix = bayesmix::read_eigen_matrix(args.get<std::string>("--hier-cov-file"));
  }
  Eigen::MatrixXd adj_matrix;
  if (args.get<std::string>("--adj-matrix-file") != EMPTYSTR) {
    adj_matrix = bayesmix::read_eigen_matrix(args.get<std::string>("--adj-matrix-file"));
  }
  auto shard_allocation = spatialcmc::read_shard_allocation_file(args.get<std::string>("--shard-assignment-file"));

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
  // std::cout << "Hierarchy Type: " << typeof(*hierarchy) << std::endl;
  
  // Create mixing object
  std::shared_ptr<AbstractMixing> mixing;
  if(args.get<std::string>("--mix-type") == "sPP") {
    mixing = std::make_shared<sPPMixing>();
  } else {
    auto & factory_mixing = MixingFactory::Instance();
    mixing = factory_mixing.create_object(args.get<std::string>("--mix-type"));
  }
  bayesmix::read_proto_from_file(args.get<std::string>("--mix-prior-file"), mixing->get_mutable_prior());

  // Split data in shards
  unsigned int num_shards = std::set<int>(shard_allocation.cbegin(), shard_allocation.cend()).size();
  std::vector<std::deque<int>> global_numbering(num_shards);
  // std::vector<Eigen::MatrixXd> data_in_shards(num_shards);
  for (size_t i = 0; i < data.size(); i++) {
    global_numbering[shard_allocation[i]].push_back(i);
  }

  // Read algorithm parameters
  bayesmix::AlgorithmParams algo_proto;
  bayesmix::read_proto_from_file(args.get<std::string>("--algo-params-file"), &algo_proto);

  // Generating shards
  std::vector<std::shared_ptr<BaseAlgorithm>> shards(num_shards);
  auto & factory_algo = AlgorithmFactory::Instance();
  std::cout << "Generating shards ... ";
  for (size_t i = 0; i < num_shards; i++) {
    shards[i] = factory_algo.create_object(algo_proto.algo_id());
    shards[i]->set_data(data(global_numbering[i], 0));
    shards[i]->set_hier_covariates(cov_matrix(global_numbering[i], Eigen::all));
    shards[i]->set_mix_covariates(adj_matrix(global_numbering[i], global_numbering[i]));
    shards[i]->set_hierarchy(hierarchy->clone());
    // shards[i]->set_hierarchy_prior(args.get<std::string>("--hier-prior-file"));
    shards[i]->set_mixing(mixing->clone());
    // shards[i]->set_mixing_prior(args.get<std::string>("--mix-prior-file"));
    shards[i]->read_params_from_proto(algo_proto);
    if(i != 0) 
      shards[i]->set_verbose(false);
  }
  std::cout << "Done" << std::endl;

  // Generating shards in parallel
  // std::vector<Shard> shards(num_shards);
  // #pragma omp parallel for num_threads(num_threads)
	// for (size_t i = 0; i < num_shards; i++) {
  //   shards[i].set_data(data(global_numbering[i], 0));
  //   if(hierarchy->is_dependent()){
  //     algo->set_hier_covariates(cov_matrix);
  //   }
  //   if(mixing->is_dependent()){
  //     shards[i].set_adjacency_matrix(adj_matrix(global_numbering[i], global_numbering[i]));
  //   }
  //   shards[i].set_sampler_parameters(args.get<std::string>("--algo-params-file"));
  //   shards[i].set_hierarchy(hierarchy->clone(), args.get<std::string>("--hier-prior-file"));
  //   shards[i].set_mixing(mixing->clone(), args.get<std::string>("--mix-prior-file"));
  //   if(i != 0) 
  //     shards[i].set_verbose(false);
  // }

  // Check - OK
  // std::cout << "N° of shards: " << shards.size() << std::endl;
  // for (auto && elem : shards) { elem.print(); std::cout << std::endl; }

  // Run MCMC samplers in each shard - in parallel
  std::vector<std::vector<bayesmix::AlgorithmState>> sharded_partitions(num_shards);
  #pragma omp parallel for num_threads(num_threads)
  for (size_t i = 0; i < num_shards; i++) {
    run_mcmc(shards[i], sharded_partitions[i]);
  }
  
  // std::vector<spatialcmc::MCMCChain> chains(num_shards);
  // #pragma omp parallel for num_threads(num_threads)
  // for (size_t i = 0; i < num_shards; i++) {
  //   chains[i] = shards[i].sample();
  // }

  // Set-up shard merger
  LocalClusterMerger<AbstractHierarchy> partition_merger(sharded_partitions);
  partition_merger.set_hierarchy(hierarchy->clone());
  partition_merger.set_data(&data);
  partition_merger.set_cov_matrix(&cov_matrix);
  partition_merger.set_adj_matrix(&adj_matrix);
  partition_merger.set_global_numbering(global_numbering);
  partition_merger.set_seed(algo_proto.rng_seed());
  partition_merger.set_hierarchy_prior(args.get<std::string>("--hier-prior-file"));
  // ShardMerger shard_merger(shards, chains, global_numbering, adj_matrix);
  // std::vector<bayesmix::AlgorithmState> state_vect(shard_merger.get_num_iter());
  // progresscpp::ProgressBar* bar = new progresscpp::ProgressBar(shard_merger.get_num_iter(), 60);

  // Merging shards in parallel
  std::cout << "Merging MCMC chains..." << std::endl;
  std::vector<bayesmix::AlgorithmState> merged_states(partition_merger.get_num_iter());
  progresscpp::ProgressBar* bar = new progresscpp::ProgressBar(partition_merger.get_num_iter(), 60);
  #pragma omp parallel for num_threads(num_threads)
  for (size_t i = 0; i < merged_states.size(); i++) {
    // std::cout << "ITERATION: " << i << std::endl;
    merged_states[i].CopyFrom(partition_merger.merge(i));
		++(*bar);
    #pragma omp critical
		bar->display();
  }
	bar->done();
	std::cout << "Done" << std::endl;

  // Serialization of the final chain
  spatialcmc::MCMCChain merged_chain;
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
  std::cout << "End of run_cmc.cc" << std::endl;
  std::cout << std::endl;

  // Memory management
  google::protobuf::ShutdownProtobufLibrary();
  delete bar;

  // Terminate program with no error
  return 0;
};
