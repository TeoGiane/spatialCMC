#include "poisson_gamma_updater.h"

std::shared_ptr<AbstractUpdater> PoissonGammaUpdater::clone() const {
    auto out = std::make_shared<PoissonGammaUpdater>(static_cast<PoissonGammaUpdater const &>(*this));
    out->clear_hypers();
    return out;
  }

AbstractUpdater::ProtoHypersPtr PoissonGammaUpdater::compute_posterior_hypers(
    AbstractLikelihood& like, AbstractPriorModel& prior) {
  // Likelihood and Prior downcast
  auto& likecast = downcast_likelihood(like);
  auto& priorcast = downcast_prior(prior);
  // Getting required quantities from likelihood and prior
  int card = likecast.get_card();
  double data_sum = likecast.get_data_sum();
  auto hypers = priorcast.get_hypers();
  // No update possible
  if (card == 0) {
    return priorcast.get_hypers_proto();
  }
  // Compute posterior hyperparameters
  double shape_post = hypers.shape + data_sum;
  double rate_post = hypers.rate + card;
  // Proto conversion
  ProtoHypers out;
  bayesmix::GammaDistribution unp_hypers;
  unp_hypers.set_shape(shape_post);
  unp_hypers.set_rate(rate_post);
  out.mutable_custom_state()->PackFrom(unp_hypers);
  return std::make_shared<ProtoHypers>(out);
}