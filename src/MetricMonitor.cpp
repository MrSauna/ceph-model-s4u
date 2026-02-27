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

  auto all_links = sg4::Engine::get_instance()->get_all_links();

  double dt = 0.01; // Sleep for 10ms (High frequency)
  double time_passed = 0.0;
  std::map<const simgrid::s4u::Link *, double> accumulated_loads;

  while (true) {
    sg4::this_actor::sleep_for(dt);
    time_passed += dt;

    for (const auto *link : all_links) {
      accumulated_loads[link] += link->get_load() * dt;
    }

    if (time_passed >= interval_) {
      double now = sg4::Engine::get_clock();
      for (const auto *link : all_links) {
        // Log average utilization over the time window
        double avg_load = accumulated_loads[link] / interval_;
        std::string name = link->get_cname();
        file << now << ",link," << name << "," << avg_load << "\n";
        accumulated_loads[link] = 0.0; // Reset
      }
      time_passed = 0.0;
    }
  }
}
