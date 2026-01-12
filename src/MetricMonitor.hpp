#ifndef METRIC_MONITOR_HPP
#define METRIC_MONITOR_HPP

#include <simgrid/s4u.hpp>
#include <string>
#include <vector>

class MetricMonitor {
  double interval_;
  std::string filename_;
  std::vector<double> cumulative_load_;

public:
  explicit MetricMonitor(double interval, std::string filename);
  void operator()();

  void setup_callbacks();
};

#endif // METRIC_MONITOR_HPP
