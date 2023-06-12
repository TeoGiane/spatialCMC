#ifndef SPP_MIXING_H_
#define SPP_MIXING_H_

#include <google/protobuf/message.h>
#include <google/protobuf/stubs/casts.h>

#include <memory>
#include <stan/math/rev.hpp>
#include <stan/math/prim/prob.hpp>
#include <vector>

#include "spp_prior.pb.h"
#include "src/hierarchies/abstract_hierarchy.h"
#include "src/mixings/base_mixing.h"
#include "src/utils/rng.h"

namespace sPP {
  struct State {
    double totalmass, logtotmass, lambda;//, alpha, beta, gamma_e, gamma_ne;
    // Eigen::MatrixXd prox_matrix;
  };
};  // namespace sPP

class sPPMixing
    : public BaseMixing<sPPMixing, sPP::State, spatialcmc::sPPPrior> {
 public:

  sPPMixing() = default;
  ~sPPMixing() = default;

  //! Initializes state parameters to appropriate values
  void initialize_state() override;

  //! Performs conditional update of state, given allocations and unique values
  //! @param unique_values  A vector of (pointers to) Hierarchy objects
  //! @param allocations    A vector of allocations label
  void update_state(
      const std::vector<std::shared_ptr<AbstractHierarchy>> &unique_values,
      const std::vector<unsigned int> &allocations) override;

  //! Read and set state values from a given Protobuf message
  void set_state_from_proto(const google::protobuf::Message &state_) override;

  //! Writes current state to a Protobuf message and return a shared_ptr
  //! New hierarchies have to first modify the field 'oneof val' in the
  //! MixingState message by adding the appropriate type
  std::shared_ptr<bayesmix::MixingState> get_state_proto() const override;

  //! Returns the Protobuf ID associated to this class
  bayesmix::MixingId get_id() const override { return bayesmix::MixingId::UNKNOWN_MIXING; }

  //! Returns whether the mixing is conditional or marginal
  bool is_conditional() const override { return false; }

  //! Returns whether the mixing depends on covariates
  bool is_dependent() const override { return true; }

 protected:
  //! Returns probability mass for an old cluster (for marginal mixings only)
  //! @param n          Total dataset size
  //! @param n_clust    Number of clusters
  //! @param log        Whether to return logarithm-scale values or not
  //! @param propto     Whether to include normalizing constants or not
  //! @param hier       `Hierarchy` object representing the cluster
  //! @param covariate  Covariate vector
  //! @return           Probability value
  double mass_existing_cluster(
      const unsigned int n, const unsigned int n_clust, const bool log,
      const bool propto, const std::shared_ptr<AbstractHierarchy> hier,
      const Eigen::RowVectorXd &covariate) const override;

  //! Returns probability mass for a new cluster (for marginal mixings only)
  //! @param n          Total dataset size
  //! @param log        Whether to return logarithm-scale values or not
  //! @param propto     Whether to include normalizing constants or not
  //! @param n_clust    Current number of clusters
  //! @param covariate  Covariate vector
  //! @return           Probability value
  double mass_new_cluster(const unsigned int n, const unsigned int n_clust,
                          const bool log, const bool propto,
                          const Eigen::RowVectorXd &covariate) const override;

};

#endif  // SPP_MIXING_H_
