#ifndef BAYESMIX_HIERARCHIES_UPDATERS_EMPTY_UPDATER_H_
#define BAYESMIX_HIERARCHIES_UPDATERS_EMPTY_UPDATER_H_

#include "semi_conjugate_updater.h"
#include "src/hierarchies/likelihoods/empty_likelihood.h"
#include "src/hierarchies/priors/empty_prior_model.h"

/**
 * Empty Updater
 */

class EmptyUpdater
    : public SemiConjugateUpdater<EmptyLikelihood, EmptyPriorModel> {
 public:
  EmptyUpdater() = default;
  ~EmptyUpdater() = default;

  bool is_conjugate() const override { return true; };

  ProtoHypersPtr compute_posterior_hypers(AbstractLikelihood &like,
                                          AbstractPriorModel &prior) override;

  std::shared_ptr<AbstractUpdater> clone() const override {
    auto out = std::make_shared<EmptyUpdater>(static_cast<EmptyUpdater const &>(*this));
    out->clear_hypers();
    return out;
  }
};

#endif  // BAYESMIX_HIERARCHIES_UPDATERS_NNIG_UPDATER_H_
