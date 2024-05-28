#include <gtest/gtest.h>

#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>

#include "hierarchies/poisson_regression_hierarchy.h"

TEST(pois_reg_hierarchy, update_beta) {
    // Create hierachy object
    auto hier = std::make_shared<PoissonRegHierarchy>();
    spatialcmc::PoissonGammaPrior prior;
    prior.mutable_fixed_values()->set_shape(3.0);
    prior.mutable_fixed_values()->set_rate(2.0);
    hier->get_mutable_prior()->CopyFrom(prior);
    hier->initialize();
    
    // Set common regression coefficients
    Eigen::VectorXd common_beta = Eigen::VectorXd::Zero(2);
    hier->set_reg_coeffs(&common_beta);

    // Clone hierarchy
    auto hier2 = hier->clone();

    // Check regression coefficients
    std::cout << "First: " << std::endl;
    std::cout << "common_beta: " << common_beta.transpose() << std::endl;
    std::cout << "hier->reg_coeffs(): " << hier->get_reg_coeffs().transpose() << std::endl;
    std::cout << "hier2->reg_coeffs(): " << std::static_pointer_cast<PoissonRegHierarchy>(hier2)->get_reg_coeffs().transpose() << std::endl;
    // hier->get_mutable_prior()->PrintDebugString();

    common_beta = Eigen::VectorXd::Ones(2);

    std::cout << "Second: " << std::endl;
    std::cout << "common_beta: " << common_beta.transpose() << std::endl;
    std::cout << "hier->reg_coeffs(): " << hier->get_reg_coeffs().transpose() << std::endl;
    std::cout << "hier2->reg_coeffs(): " << std::static_pointer_cast<PoissonRegHierarchy>(hier2)->get_reg_coeffs().transpose() << std::endl;


    // mutable_custom_state()->PackFrom(unp_hypers);
    // hier->set_reg_coeffs(*regression_coefficients);
    // partition_updater->set_hierarchy(hier);
    // partition_updater->set_hier_covariates(hier_covariates);
    // sistemare roba
    bool success = true;
    ASSERT_TRUE(success);
}