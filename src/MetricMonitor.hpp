#ifndef METRIC_MONITOR_HPP
#define METRIC_MONITOR_HPP

#include <simgrid/s4u.hpp>
#include <string>

class MetricMonitor {
  double interval_;
  std::string filename_;

public:
  explicit MetricMonitor(double interval, std::string filename);
  void operator()();
};

#endif // METRIC_MONITOR_HPP
