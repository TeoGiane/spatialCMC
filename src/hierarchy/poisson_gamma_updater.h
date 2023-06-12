#ifndef POISSON_GAMMA_UPDATER_H_
#define POISSON_GAMMA_UPDATER_H_

#include "poisson_likelihood.h"
#include "gamma_prior_model.h"
#include "src/hierarchies/updaters/semi_conjugate_updater.h"

class PoissonGammaUpdater : public SemiConjugateUpdater<PoissonLikelihood, GammaPriorModel> {
 public:
  PoissonGammaUpdater() = default;
  ~PoissonGammaUpdater() = default;

  bool is_conjugate() const override { return true; };

  ProtoHypersPtr compute_posterior_hypers(AbstractLikelihood& like, AbstractPriorModel& prior) override;

  std::shared_ptr<AbstractUpdater> clone() const override;
};

#endif // POISSON_GAMMA_UPDATER_H_