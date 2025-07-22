#ifndef ALGORITHMS_SAMPLER_H_
#define ALGORITHMS_SAMPLER_H_

#include <lib/progressbar/progressbar.h>

#include <memory>
#include <stan/math/rev.hpp>
#include <vector>

// #include "algorithm_id.pb.h"
#include "algorithm_params.pb.h"
#include "poisson_regression_algorithm_state.pb.h"
#include "poisson_gamma_prior.pb.h"
#include "spp_prior.pb.h"
// #include "algorithm_state.pb.h"


#include "src/algorithms/neal2_algorithm.h"
#include "src/collectors/base_collector.h"
//#include "src/mixings/dirichlet_mixing.h"

#include "hierarchies/poisson_regression_hierarchy.h" 
#include "mixing/spp_mixing.h"
// #include "utils/mixing_utils.h"

// #include "src/hierarchies/base_hierarchy.h"
// #include "src/hierarchies/nnig_hierarchy.h"
// #include "src/mixings/abstract_mixing.h"
// #include "src/algorithms/base_algorithm.h"

class PoissonRegressionSampler {
 private:
  
  // ALGORITHM PARAMS
  unsigned int maxiter = 1000;
  unsigned int burnin = 100;
  unsigned int init_num_clusters = 3;

  // Metropolis Hastings parameters
  int iter = 0;//, n_accepted;
  double cov_size = 0;
  // Eigen::MatrixXd Sigma;
  // Eigen::MatrixXd cov_n_beta;
  // Eigen::VectorXd beta_mean;

  // double n_accepted = 0;
  // double step_size = 0.05;
  Eigen::RowVectorXd n_accepted;
  // Eigen::VectorXd accepted_in_batch;
  Eigen::RowVectorXd log_sigmas;
  int batch_size = 50;
  
  // will be useless
  // Eigen::VectorXd sigmas;
  // Eigen::VectorXd beta_sum;
  // Eigen::VectorXd beta_sum_sq;
  // Eigen::VectorXd sigma_n_beta;

  // PROTO FILES WITH PRIOR PARAMETERS
  std::string hier_prior_file;
  std::string mix_prior_file;

  // DATA
  Eigen::MatrixXd data;
  Eigen::MatrixXd hier_covariates;
  Eigen::MatrixXd adj_matrix;

  // CURRENT STATE
  spatialcmc::PoissonRegAlgorithmState curr_state;

  // ALGORITHM STATE
  Eigen::VectorXd regression_coefficients;
  std::shared_ptr<Neal2Algorithm> partition_updater = std::make_shared<Neal2Algorithm>();

  // MISCELLANEOUS
  //! Turns on or off the descriptive output of the class methods
  bool verbose = true;
  bool mh_verbose = true;
 
 public:
  
  PoissonRegressionSampler() = default;
  ~PoissonRegressionSampler() = default;

  //! Runs the algorithm and saves the MCMC chain to the given `Collector`
  void run(BaseCollector *collector);

  unsigned int get_maxiter() const { return maxiter; };

  unsigned int get_burnin() const { return burnin; };

  //! Returns the Protobuf ID associated to this class
  // virtual bayesmix::AlgorithmId get_id() const = 0;

  void set_maxiter(const unsigned int maxiter_) { maxiter = maxiter_; };

  void set_burnin(const unsigned int burnin_) { burnin = burnin_; };

  void set_init_num_clusters(const unsigned int init_) {
    init_num_clusters = init_;
  };

  void set_hierarchy_prior(const std::string & _hier_prior_file) {
    hier_prior_file = _hier_prior_file;
  };

  void set_mixing_prior(const std::string & _mix_prior_file) {
    mix_prior_file = _mix_prior_file;
  };

  void set_data(const Eigen::MatrixXd &data_) {
    data = data_;
    // partition_updater->set_data(data_);
  }

  void set_hier_covariates(const Eigen::MatrixXd & _hier_covariates) {
    hier_covariates = _hier_covariates;
    // partition_updater->set_hier_covariates(_hier_covariates);
  };

  void set_adjacency_matrix(const Eigen::MatrixXd & _adj_matrix) {
    adj_matrix = _adj_matrix;
    // partition_updater->set_mix_covariates(_adj_matrix);
  };

  // void set_mixing(const std::shared_ptr<AbstractMixing> mixing_) {
  //   mixing = mixing_;
  // }

  // void set_hierarchy(const std::shared_ptr<AbstractHierarchy> hier_) {
  //   unique_values.clear();
  //   unique_values.push_back(hier_);
  // }

  void set_verbose(const bool verbose_) {
    verbose = verbose_;
    mh_verbose = verbose_;
    partition_updater->set_verbose(verbose_);
  };

  //! Reads and sets algorithm parameters from an appropriate Protobuf message
  void read_params_from_proto(const bayesmix::AlgorithmParams &params) {
    // Generic parameters
    maxiter = params.iterations();
    burnin = params.burnin();
    init_num_clusters = params.init_num_clusters();
    auto &rng = bayesmix::Rng::Instance().get();
    rng.seed(params.rng_seed());
  };

  void set_state_proto(const std::shared_ptr<google::protobuf::Message> state);

  //! Evaluates the estimated log-lpdf on the state contained in `curr_state`
  Eigen::VectorXd lpdf_from_state(const Eigen::MatrixXd &grid, unsigned int data_idx) {
  
    // Prepare buffer
    Eigen::VectorXd lpdf_final(grid.rows());

    // Extract datum cluster allocation and compute corresponding rate
    int data_alloc = curr_state.partition().cluster_allocs(data_idx);
    double rate = spatialcmc::unpack_to<spatialcmc::PoissonState>(
                    curr_state.partition().cluster_states(data_alloc).custom_state()).rate();
    double offset = hier_covariates.row(data_idx)(0);
    Eigen::RowVectorXd cov_vector = hier_covariates.row(data_idx)(Eigen::seq(1,Eigen::last));
    Eigen::VectorXd beta = bayesmix::to_eigen(curr_state.regression_coefficients());
    rate *= offset * std::exp(cov_vector * beta);
    
    // Get lpdf evaluation
    for (size_t i = 0; i < grid.rows(); i++) {
      lpdf_final(i) = stan::math::poisson_lpmf(grid(i), rate);
    }
    return lpdf_final;
  }

  // virtual std::shared_ptr<BaseAlgorithm> clone() const = 0;

  // std::vector<std::shared_ptr<AbstractHierarchy>> get_unique_values() const {
  //   return partition_updater->get_unique_values();
  // };

 protected:
 
  double reg_coeffs_target_lpdf(const Eigen::VectorXd & curr_state) {
    double lpdf = 0;
    // Get unique values and cluster allocations
    auto unique_vals = partition_updater->get_unique_values();
    auto clust_allocs = partition_updater->get_allocations();
    // Joint likelihood
    for (size_t i = 0; i < data.size(); i++) {
      double offset = hier_covariates.row(i)(0);
      Eigen::RowVectorXd cov_vector = hier_covariates.row(i)(Eigen::seq(1,Eigen::last));
      double rate = spatialcmc::unpack_to<spatialcmc::PoissonState>(unique_vals[clust_allocs[i]]->get_state_proto()->custom_state()).rate();
      lpdf += stan::math::poisson_log_lpmf(data.row(i)(0),log(offset)+log(rate)+cov_vector*curr_state);
    }
    // Prior
    for (size_t l = 0; l < cov_size; l++) {
      lpdf += stan::math::normal_lpdf(curr_state(l), 0, sqrt(10));
    }
    // Return
    return lpdf;
  }

  // double reg_coeffs_lpdf(const Eigen::VectorXd & curr_state) {
  //   // std::cout << "reg_coeffs_lpdf()" << std::endl;
  //   double lpdf = 0;
  //   auto allocs = partition_updater->get_allocations();
  //   auto unique_vals = get_unique_values();

  //   for (size_t i = 0; i < data.size(); i++) {
  //     // std::cout << "Province: " << i << std::endl;
  //     Eigen::VectorXi hr_eigen = high_resolution_data(province_idx[i],0);
  //     std::vector<int> hr_prov(hr_eigen.data(), hr_eigen.data()+hr_eigen.size());
  //     Eigen::VectorXd thetas(province_idx[i].size());
  //     for (size_t k = 0; k < province_idx[i].size(); k++) {
  //       // std::cout << "Municipality: " << k << std::endl;
  //       unsigned int idx = province_idx[i][k];
  //       thetas(k) = log(unique_vals[allocs[idx]]->get_state_proto()->general_state().data(0));
  //       thetas(k) += hier_covariates->row(idx) * curr_state;
  //     }
  //     // thetas /= thetas.sum();
  //     // std::cout << "prop_lambda: " << exp(stan::math::log_sum_eget_cxp(thetas)) << std::endl;
  //     // std::cout << "exp_version: " << stan::math::poisson_lpmf(data[i], exp(stan::math::log_sum_exp(thetas))) << std::endl;
  //     // lpdf += stan::math::poisson_log_lpmf(data[i], stan::math::log_sum_exp(thetas));
  //     lpdf += stan::math::multinomial_lpmf(hr_prov, stan::math::softmax(thetas));
  //   }
  //   for (size_t l = 0; l < curr_state.size(); l++) {
  //     lpdf += stan::math::normal_lpdf(curr_state(l), 0, sqrt(10));
  //   }
  //   return lpdf;

    // for (size_t j = 0; j < high_resolution_data.rows(); j++) {
    //   double lambda = log(unique_vals[allocs[j]]->get_state_proto()->general_state().data(0));
    //   lambda += hier_covariates->row(j) * curr_state;
    //   lpdf += stan::math::poisson_lpmf(high_resolution_data(j,0), exp(lambda));
    // }
    // for (size_t l = 0; l < curr_state.size(); l++) {
    //   lpdf += stan::math::normal_lpdf(curr_state(l), 0, 10);
    // }
  // };

  double reg_coeffs_prop_lpdf(const Eigen::VectorXd & prop_state,
                              const Eigen::VectorXd & curr_state) {
    // Computing lpdf
    double lpdf = 0;
    for (size_t l = 0; l < cov_size; l++) {
      lpdf += stan::math::normal_lpdf(prop_state(l), curr_state(l), exp(log_sigmas(l)));
    }
    return lpdf;
  };

  void compute_sigma() {
    // std::cout << "ITER: " << iter << std::endl;
    // Starting point
    if (iter == 0)
      log_sigmas = std::log(0.05) * Eigen::VectorXd::Ones(cov_size);
    // Increment iteration counter
    iter++;
    // Adapt variances every batch_size iterations
    if (iter % batch_size == 0) {
      double adapt = std::min(0.05, 1/sqrt(iter));
      for (size_t l = 0; l < cov_size; l++) {
        if (n_accepted(l) / iter >= 0.44)
          log_sigmas(l) += adapt;
        else
          log_sigmas(l) -= adapt;
      }
    }
    // Check
    // std::cout << "sigmas: " << log_sigmas.array().exp() << std::endl;
  };

  // ALGORITHM FUNCTIONS
  //! Initializes all members of the class before running the algorithm
  void initialize() {

    // Get random seed
    auto &rng = bayesmix::Rng::Instance().get();

    // Check input data were provided
    if (data.rows() == 0) {
      throw std::invalid_argument("Data was not provided to algorithm");
    }
    if (hier_covariates.rows() == 0) {
      throw std::invalid_argument("Covariates was not provided to algorithm");
    } else {
      if(hier_covariates.rows() != data.rows()) {
        throw std::invalid_argument("Sizes of data and hierarchy covariates do not match");
      }
    }
    if(adj_matrix.rows() == 0) {
      throw std::invalid_argument("Adjacency Matrix was not provided to algorithm");
    } else {
      if(adj_matrix.rows() != adj_matrix.cols()){
        throw std::invalid_argument("Adjacency Matrix is not a square matrix");
      }
      if(adj_matrix.rows() != data.rows()){
        throw std::invalid_argument("Sizes of data and adjacency matrix do not match");
      }
    }

    // Initialize regression coefficients
    cov_size = hier_covariates.cols() - 1;
    regression_coefficients = Eigen::VectorXd::Zero(cov_size);
    
    // Resize MH parameter
    log_sigmas.resize(cov_size);
    n_accepted = Eigen::VectorXd::Zero(cov_size);

    // Initialize partition updater
    // 1. Hierarchy set-up
    auto hier = std::make_shared<PoissonRegHierarchy>();
    bayesmix::read_proto_from_file(hier_prior_file, hier->get_mutable_prior());
    hier->set_covariates(hier_covariates);
    hier->set_reg_coeffs(regression_coefficients);
    partition_updater->set_hierarchy(hier);
    partition_updater->set_hier_covariates(hier_covariates);
    // 2. Mixing Set-up
    auto mixing = std::make_shared<sPPMixing>();
    bayesmix::read_proto_from_file(mix_prior_file, mixing->get_mutable_prior());
    partition_updater->set_mixing(mixing);
    partition_updater->set_mix_covariates(adj_matrix);
    // 3. Algorithm Set-up
    partition_updater->set_data(data);
    partition_updater->set_init_num_clusters(init_num_clusters);
    // 4. Initialize partition_updater
    partition_updater->initialize();
  };

  //! Prints a message at the beginning of `run()`
  void print_startup_message() const {
    std::cout << "Running Poisson Regression sampler... " << std::endl;
  };

  //! Performs Gibbs sampling sub-step for partition update (usign updater)
  void sample_partition() { partition_updater->step(); };

  void sample_regression_coefficients() {
    // Define random seed and increase iteration counter
    auto &rng = bayesmix::Rng::Instance().get();
    // Store copy of current and next state
    Eigen::VectorXd curr_state = regression_coefficients;
    Eigen::VectorXd next_state = curr_state;
    // Compute adaptive standard deviations
    compute_sigma();
    // Conditional adaptive MH for each component of the vector
    for (size_t l = 0; l < cov_size; l++) {
      // Proposal for l-th component
      Eigen::VectorXd prop_state = curr_state;
      prop_state(l) += stan::math::normal_rng(0, exp(log_sigmas(l)), rng);
      // Compute log acceptance rate
      double log_arate = reg_coeffs_target_lpdf(prop_state) -
                         reg_coeffs_target_lpdf(curr_state) +
                         reg_coeffs_prop_lpdf(curr_state, prop_state) -
                         reg_coeffs_prop_lpdf(prop_state, curr_state);
      // Test for acceptance + update
      if (std::log(stan::math::uniform_rng(0, 1, rng)) < log_arate) {
        n_accepted(l)++; next_state(l) = prop_state(l);
      }
    }
    // Update accordingly
    regression_coefficients = next_state;
  };

  // void update_reg_coeffs(const Eigen::VectorXd & _reg_coeffs) {
  //   // (*regression_coefficients) = _reg_coeffs;
  //   // std::static_pointer_cast<PoissonGammaHierarchy>(partition_updater->get_unique_values()[0])->set_reg_coeffs(_reg_coeffs);
  // }

  //! Prints a message at the end of `run()`
  void print_ending_message() const {
    if (mh_verbose) {
      std::cout << "Accept. Rates: " << n_accepted / iter << std::endl;
    }
    std::cout << "Done" << std::endl;
  };

  //! Saves the current iteration's state in Protobuf form to a `Collector`
  void save_state(BaseCollector *collector, unsigned int iter) {
    collector->collect(get_state_as_proto(iter));
  }

  //! Performs a single step of algorithm
  void step() {
    sample_regression_coefficients();
    for (auto & clust : partition_updater->get_unique_values()) {
      auto clustcast = std::static_pointer_cast<PoissonRegHierarchy>(clust);
      clustcast->set_reg_coeffs(regression_coefficients);
    }
    sample_partition();
  }

  // AUXILIARY TOOLS
  //! Returns Protobuf object containing current state values and iter number
  spatialcmc::PoissonRegAlgorithmState get_state_as_proto(const unsigned int iter) {
    spatialcmc::PoissonRegAlgorithmState iter_out;
    bayesmix::to_proto(regression_coefficients, iter_out.mutable_regression_coefficients());
    iter_out.mutable_partition()->CopyFrom(partition_updater->get_state_as_proto(iter));
    return iter_out;
  };

  //! Advances `Collector` reading pointer by one, and returns 1 if successful
  bool update_state_from_collector(BaseCollector *const coll) {
    bool success = coll->get_next_state(&curr_state);
    return success;
  };

  /*
    Returns whether the posterior parameters for the hierarchies should be
     updated each time an observation is added or removed from the cluster.
     This can potentially reduce computational effort in algorithms which do
     not need this kind of update. If false, this update is usually performed
     instead during the `sample_unique_values()` substep, and viceversa.
  */
  //virtual bool update_hierarchy_params() { return false; }
};

#endif  // ALGORITHMS_SAMPLER_H_

// Cosa inutile ma già scritta
  // // DENSITY ESTIMATE FUNCTIONS
  // //! Evaluates the estimated log-pdf on all iterations over a grid of points
  // //! @param collector   `Collector` object from which the MCMC is read
  // //! @param grid   Grid of row points on which the density is to be evaluated
  // //! @param hier_covariate   Optional covariates related to the `Hierarchy`
  // //! @param mix_covariate   Optional covariates related to the `Mixing`
  // //! @return   The estimation over all iterations (rows) and points (cols)
  // Eigen::MatrixXd eval_lpdf(BaseCollector *const collector,
  //                           const Eigen::MatrixXd &grid,
  //                           unsigned int data_idx) {
  //   std::deque<Eigen::VectorXd> lpdf;
  //   bool keep = true;
  //   progresscpp::ProgressBar *bar = nullptr;
  //   if (verbose) {
  //     bar = new progresscpp::ProgressBar(collector->get_size(), 60);
  //   }
  //   while (keep) {
  //     keep = update_state_from_collector(collector);
  //     if (!keep) {
  //       break;
  //     }
  //     lpdf.push_back(lpdf_from_state(grid, data_idx));
  //     if (verbose) {
  //       ++(*bar);
  //       bar->display();
  //     }
  //   }
  //   collector->reset();
  //   if (verbose) {
  //     bar->done();
  //     delete bar;
  //     print_ending_message();
  //   }
  //   return bayesmix::stack_vectors(lpdf);
  // }