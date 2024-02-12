#include "empty_updater.h"

// #include "src/hierarchies/likelihoods/states/includes.h"
// #include "src/hierarchies/priors/hyperparams.h"

AbstractUpdater::ProtoHypersPtr EmptyUpdater::compute_posterior_hypers(AbstractLikelihood& like, AbstractPriorModel& prior) {
  // Proto conversion
  ProtoHypers out;
  spatialcmc::EmptyDistribution unp_hypers;
  out.mutable_custom_state()->PackFrom(unp_hypers);
  return std::make_shared<ProtoHypers>(out);
}
