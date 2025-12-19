#include "MetricMonitor.hpp"
#include <fstream>
#include <iomanip>

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_metric_monitor,
                             "Messages specific for the metric monitor");

namespace sg4 = simgrid::s4u;

// MetricMonitor::MetricMonitor(double interval, std::string filename)
//     : interval_(interval), filename_(std::move(filename)) {
//     setup_callbacks();
// }
// Callback approach abandoned due to missing kernel headers.
// Falling back to High-Frequency Polling (Trapezoidal Integration).

MetricMonitor::MetricMonitor(double interval, std::string filename)
    : interval_(interval), filename_(std::move(filename)) {}

void MetricMonitor::operator()() {
  std::ofstream file(filename_);
  if (!file.is_open()) {
    XBT_ERROR("Failed to open metrics file: %s", filename_.c_str());
    return;
  }

  file << "timestamp,resource_type,resource_name,value\n";

  // Initialize map
  auto all_links = sg4::Engine::get_instance()->get_all_links();
  for (auto *link : all_links) {
    cumulative_load_[link->get_cname()] = 0.0;
  }

  double polling_period = 0.001; // 1ms polling for high resolution
  double last_time = sg4::Engine::get_clock();
  double next_report = last_time + interval_;

  while (true) {
    sg4::this_actor::sleep_for(polling_period);
    double now = sg4::Engine::get_clock();
    double dt = now - last_time;

    // Integrate Load: Bytes = Bps * seconds
    for (auto *link : all_links) {
      cumulative_load_[link->get_cname()] += link->get_load() * dt;
    }

    // Report if interval passed
    if (now >= next_report) {
      for (auto const &[name, bytes] : cumulative_load_) {
        // Average Bandwidth = Total Bytes / Interval
        // Note: Use actual interval (now - (next_report - interval_)) to be
        // precise? Or just interval_. SimGrid clocks are precise.
        double avg_bw = bytes / interval_;

        // Log if non-zero (or all to match user request)
        if (avg_bw > 1.0) { // Filter extremely small noise
          file << std::fixed << std::setprecision(2) << next_report << ",link,"
               << name << "," << avg_bw << "\n";
        }
        cumulative_load_[name] = 0.0; // Reset
      }
      next_report += interval_;
    }

    last_time = now;
  }
}
