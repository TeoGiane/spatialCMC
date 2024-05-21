#include "sampler.h"

void BNPDownscalingSampler::run(BaseCollector *collector) {
  initialize();
  if (verbose) {
    print_startup_message();
  }
  unsigned int iter = 0;
  collector->start_collecting();
  progresscpp::ProgressBar *bar = nullptr;
  if (verbose) {
    bar = new progresscpp::ProgressBar(maxiter, 60);
  }
  // Main loop
  while (iter < maxiter) {
    // std::cout << "Iteration: " << iter << std::endl;
    step();    
    if (iter >= burnin) {
      save_state(collector, iter);
    }
    iter++;
    if (verbose) {
      ++(*bar);
      bar->display();
    }
  }
  collector->finish_collecting();
  if (verbose) {
    bar->done();
    delete bar;
    print_ending_message();
    mh_verbose = false;
  }
}

void BNPDownscalingSampler::set_data(const Eigen::VectorXi & data_) { data = data_; }

void BNPDownscalingSampler::set_state_proto(const std::shared_ptr<google::protobuf::Message> state) {
  curr_state.CopyFrom(google::protobuf::internal::down_cast<bayesmix::BNPDownscalingState &>(*state.get()));
}