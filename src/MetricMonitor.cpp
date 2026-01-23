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

  while (true) {
    sg4::this_actor::sleep_for(interval_);
    double now = sg4::Engine::get_clock();

    for (const auto *link : all_links) {
      double current_load = link->get_load();

      // Log if non-zero
      std::string name = link->get_cname();
      file << now << ",link," << name << "," << current_load << "\n";
    }
  }
}
