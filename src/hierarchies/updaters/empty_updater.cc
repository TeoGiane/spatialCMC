#include "empty_updater.h"

#include "src/hierarchies/likelihoods/states/includes.h"
#include "src/hierarchies/priors/hyperparams.h"

AbstractUpdater::ProtoHypersPtr EmptyUpdater::compute_posterior_hypers(AbstractLikelihood& like, AbstractPriorModel& prior) {
  // Proto conversion
  ProtoHypers out;
  out.mutable_empty_state()->set_fake_field(0);
  return std::make_shared<ProtoHypers>(out);
}
