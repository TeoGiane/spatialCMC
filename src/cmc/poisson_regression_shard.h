#ifndef SPATIALCMC_CMC_POISSON_REGRESSION_SHARD_H_
#define SPATIALCMC_CMC_POISSON_REGRESSION_SHARD_H_

#include <memory>

#include <src/includes.h>
#include <stan/math.hpp>

#include "spatialcmc_utils.h"

#include "hierarchies/poisson_regression_hierarchy.h"
#include "mixing/spp_mixing.h"
#include "algorithm/poisson_regression_sampler.h"
#include "poisson_regression_algorithm_state.pb.h"

class PoissonRegShard {
 private:

	// What is required for the proper initialization of the algorithm
	bayesmix::AlgorithmParams sampler_parameters;

	// What is needed to run mcmc in shard (per ora fisso, opportuna estensione bayesmix)
	std::shared_ptr<PoissonRegHierarchy> hierarchy = std::make_shared<PoissonRegHierarchy>();
	std::shared_ptr<sPPMixing> mixing;
	std::shared_ptr<PoissonRegressionSampler> sampler;

	// Data we need to store in each shard
	Eigen::MatrixXd data;
	Eigen::MatrixXd adj_matrix;
	Eigen::MatrixXd hier_covariates;

	// Chain output
	spatialcmc::PoisRegMCMC chain;

	// Verbose parameter
	bool verbose = true;

 public:

	// Constructor and destructor
	PoissonRegShard() = default;
	~PoissonRegShard() = default;

	// Setters
	void set_data(const Eigen::MatrixXd & _data);

	void set_adjacency_matrix(const Eigen::MatrixXd & _adj_matrix) {
		adj_matrix = _adj_matrix;
		return;
	};

	void set_hierarchy_covariates(const Eigen::MatrixXd & _hier_covariates);
	void set_hierarchy_prior(const std::string & _hierarchy_prior);
	void set_mixing_prior(const std::string & _mixing_prior);
	void set_sampler_parameters(const std::string & _sampler_parameters_proto);
	void set_verbose(bool _verbose);

	// DEBUG
	std::string print_hier_prior() const;
	std::string print_mix_prior() const;

	// Getters
	Eigen::MatrixXd get_data() const;
	unsigned int get_seed() const { return sampler_parameters.rng_seed(); };
	std::shared_ptr<AbstractHierarchy> get_hierarchy() const;
	Eigen::MatrixXd get_adjacency_matrix() const { return adj_matrix; };

	// Utilities
	void initialize_sampler();
	spatialcmc::PoisRegMCMC sample();

	void print() {
		std::cout << "Shard debug info" << std::endl;
		std::cout << "data: " << data.transpose() << std::endl;
		std::cout << "adj_matrix: (" << adj_matrix.rows() << "," << adj_matrix.cols() << "):" << std::endl;
		std::cout << adj_matrix << std::endl;
		std::cout << "hier_id: " << hierarchy->get_id() << std::endl;
		std::cout << "hier_prior:\n" << print_hier_prior() << std::endl;
		std::cout << "mix_prior:\n" << print_mix_prior() << std::endl;
		std::cout << "mix_id: " << mixing->get_id() << std::endl;
		std::cout << "Algorithm Params:\n" << sampler_parameters.DebugString() << std::endl;
	};
};

#endif // SPATIALCMC_CMC_POISSON_REGRESSION_SHARD_H_
