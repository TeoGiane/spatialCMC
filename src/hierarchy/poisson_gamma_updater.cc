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
  // double exp_cov_sum = likecast.get_exp_cov_sum();
  // std::cout << "exp_cov_sum: " << exp_cov_sum << std::endl;
  auto hypers = priorcast.get_hypers();

  // No update possible
  if (card == 0) {
    return priorcast.get_hypers_proto();
  }
  
  // Compute posterior hyperparameters
  double shape_post = hypers.shape + data_sum;
  double rate_post = hypers.rate + card; //exp_cov_sum;

  // Proto conversion
  ProtoHypers out;
  out.mutable_general_state()->add_data(shape_post);
  out.mutable_general_state()->add_data(rate_post);
  return std::make_shared<ProtoHypers>(out);
}