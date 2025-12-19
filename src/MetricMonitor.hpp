#ifndef METRIC_MONITOR_HPP
#define METRIC_MONITOR_HPP

#include <map>
#include <mutex>
#include <simgrid/s4u.hpp>
#include <string>

class MetricMonitor {
  double interval_;
  std::string filename_;

  // Cumulative bytes transferred per link name
  std::map<std::string, double> cumulative_load_;
  std::mutex mutex_;

public:
  explicit MetricMonitor(double interval, std::string filename);
  void operator()();

  void setup_callbacks();
};

#endif // METRIC_MONITOR_HPP
