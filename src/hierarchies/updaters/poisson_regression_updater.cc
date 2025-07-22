#include "poisson_regression_updater.h"

std::shared_ptr<AbstractUpdater> PoissonRegUpdater::clone() const {
    auto out = std::make_shared<PoissonRegUpdater>(static_cast<PoissonRegUpdater const &>(*this));
    out->clear_hypers();
    return out;
  }

AbstractUpdater::ProtoHypersPtr PoissonRegUpdater::compute_posterior_hypers(
    AbstractLikelihood& like, AbstractPriorModel& prior) {
  // std::cout << "PoissonRegUpdater::compute_posterior_hypers()" << std::endl;
  // Likelihood and Prior downcast
  auto& likecast = downcast_likelihood(like);
  auto& priorcast = downcast_prior(prior);
  // Getting required quantities from likelihood and prior
  int card = likecast.get_card();
  // std::cout << "card: " << card << std::endl;
  double data_sum = likecast.get_data_sum();
  // std::cout << "data_sum: " << data_sum << std::endl;
  double cov_exp_sum = likecast.get_cov_exp_sum();
  // std::cout << "cov_exp_sum: " << cov_exp_sum << std::endl;
  auto hypers = priorcast.get_hypers();
  // No update possible
  if (card == 0) {
    return priorcast.get_hypers_proto();
  }
  // Compute posterior hyperparameters
  double shape_post = hypers.shape + data_sum;
  double rate_post = hypers.rate + cov_exp_sum;
  // Proto conversion
  ProtoHypers out;
  bayesmix::GammaDistribution unp_hypers;
  unp_hypers.set_shape(shape_post);
  unp_hypers.set_rate(rate_post);
  out.mutable_custom_state()->PackFrom(unp_hypers);
  return std::make_shared<ProtoHypers>(out);
}