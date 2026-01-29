#include "SimContext.hpp"
#include <fstream>
#include <functional>
#include <sstream>

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_sim_context,
                             "Messages specific for the simulation context");

namespace sg4 = simgrid::s4u;

// --- Helper Functions ---

std::vector<std::string> split_str(const std::string &s, char delim) {
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
    elems.push_back(item);
  if (elems.empty())
    elems.push_back("");
  return elems;
}

double resolve_val(std::string spec, int self_id, int parent_id,
                   double fallback) {
  if (spec.empty())
    return fallback;
  bool use_parent = false;
  if (spec[0] == '@') {
    use_parent = true;
    spec.erase(0, 1);
  }
  std::vector<double> values;
  std::stringstream ss(spec);
  std::string val_str;
  while (std::getline(ss, val_str, ',')) {
    try {
      if (!val_str.empty())
        values.push_back(std::stod(val_str));
    } catch (...) {
      return fallback;
    }
  }
  if (values.empty())
    return fallback;
  int idx = use_parent ? parent_id : self_id;
  return values[idx % values.size()];
}

// --- Validator ---

ShapeValidator::ShapeValidator() {
  name_ = "SHAPE (A:B:C)";
  func_ = [](std::string &str) {
    auto parts = split_str(str, ':');
    if (parts.size() != 3) {
      return std::string("Format must be 'A:B:C' (3 levels). Received: " + str);
    }
    return std::string();
  };
}

// --- SimContext ---

std::string SimContext::get_hierarchy_spec(const std::vector<std::string> &vec,
                                           int dc_idx, int level) {
  if (dc_idx >= static_cast<int>(vec.size()))
    return "";
  auto parts = split_str(vec[dc_idx], ':');
  // Safety for index access
  if (level >= static_cast<int>(parts.size()))
    return "";
  return parts[level];
}

// --- Builders (Internal) ---

// Level 3: Host (Index 2)
static void build_host(SimContext &ctx, sg4::NetZone *rack_zone, int dc_idx,
                       int rack_idx, int host_idx, double parent_speed) {

  std::string hostname = "host-" + std::to_string(ctx.global_host_id++);

  // Create Host
  auto *host = rack_zone->add_host(hostname, "100Gf");
  ctx.host_zones[hostname] = rack_zone;
  ctx.ordered_host_names.push_back(hostname);

  // Resolve Link Speed (Index 2: Host Link)
  std::string speed_spec = ctx.get_hierarchy_spec(ctx.speeds, dc_idx, 2);
  double speed_val = resolve_val(speed_spec, host_idx, rack_idx, parent_speed);
  std::string bw_str = std::to_string(speed_val) + "Gbps";

  // Create Link
  auto *link = rack_zone->add_split_duplex_link(hostname + "_link", bw_str)
                   ->set_latency("5us");

  // Route: Host -> Rack Router (ToR)
  rack_zone->add_route(host, nullptr, {{link, sg4::LinkInRoute::Direction::UP}},
                       true);

  // OSD Configuration (Index 2: OSD Count & Weight)
  std::string count_spec = ctx.get_hierarchy_spec(ctx.shapes, dc_idx, 2);
  int osd_count = std::stoi(count_spec.empty() ? "0" : count_spec);

  for (int o = 0; o < osd_count; ++o) {
    host->add_disk("disk" + std::to_string(o), ctx.disk_write_bandwidth,
                   ctx.disk_write_bandwidth);
  }
}

// Level 2: Rack (Index 1)
static void build_rack(SimContext &ctx, sg4::NetZone *dc_zone, int dc_idx,
                       int rack_idx, double parent_speed) {

  std::string rack_name = dc_zone->get_name() + "_r" + std::to_string(rack_idx);
  auto *rack_zone = dc_zone->add_netzone_star(rack_name);

  // Router: Top of Rack (ToR)
  auto *rack_router = rack_zone->add_router(rack_name + "_router");
  rack_zone->set_gateway(rack_router);

  // Resolve Uplink Speed (Index 1: Rack Uplink)
  std::string speed_spec = ctx.get_hierarchy_spec(ctx.speeds, dc_idx, 1);
  double speed_val = resolve_val(speed_spec, rack_idx, dc_idx, parent_speed);
  std::string bw_str = std::to_string(speed_val) + "Gbps";

  auto *uplink = rack_zone->add_split_duplex_link(rack_name + "_uplink", bw_str)
                     ->set_latency("10us");
  dc_zone->add_route(rack_zone, nullptr,
                     {{uplink, sg4::LinkInRoute::Direction::UP}}, true);

  // Host Count (Index 1)
  std::string count_spec = ctx.get_hierarchy_spec(ctx.shapes, dc_idx, 1);
  int host_count = std::stoi(count_spec.empty() ? "0" : count_spec);

  // Weight (Index 1)
  for (int h = 0; h < host_count; ++h) {
    build_host(ctx, rack_zone, dc_idx, rack_idx, h, speed_val);
  }
}

// --- Top-Level Builder ---

// Level 1: Data Center (Index 0)
void build_dc(SimContext &ctx, sg4::NetZone *world, int dc_config_idx) {
  std::string dc_name = "dc-" + std::to_string(dc_config_idx);
  auto *dc_zone = world->add_netzone_star(dc_name);

  // Router: DC Switch
  auto *dc_router = dc_zone->add_router(dc_name + "_router");
  dc_zone->set_gateway(dc_router);

  // Resolve DC Uplink Speed (Index 0: DC Uplink)
  std::string speed_spec = ctx.get_hierarchy_spec(ctx.speeds, dc_config_idx, 0);
  double speed_val = resolve_val(speed_spec, dc_config_idx, 0, 1000.0);
  std::string bw_str = std::to_string(speed_val) + "Gbps";

  auto *uplink = dc_zone->add_split_duplex_link(dc_name + "_uplink", bw_str)
                     ->set_latency("100us");

  // Route: DC Switch -> Outside (DC Boundary)
  // Route: DC Switch -> World (in World Zone)
  world->add_route(dc_zone, nullptr,
                   {{uplink, sg4::LinkInRoute::Direction::UP}}, true);

  // Rack Count (Index 0)
  std::string count_spec = ctx.get_hierarchy_spec(ctx.shapes, dc_config_idx, 0);
  int rack_count = std::stoi(count_spec.empty() ? "0" : count_spec);

  // Weight (Index 0)
  for (int r = 0; r < rack_count; ++r) {
    build_rack(ctx, dc_zone, dc_config_idx, r, speed_val);
  }
}

std::string SimContext::to_string() {
  std::stringstream ss;
  ss << "  DCs: " << shapes.size() << std::endl;
  ss << "  Racks: " << global_rack_id << std::endl;
  ss << "  Hosts: " << global_host_id << std::endl;
  ss << "  OSDs: " << global_osd_id << std::endl;
  ss << "  Start Up Delay: " << start_up_delay << std::endl;
  ss << "  Shut Down Delay: " << shut_down_delay << std::endl;
  std::string profile_str = "unknown";
  switch (profile) {
  case SchedulerProfile::BALANCED:
    profile_str = "balanced";
    break;
  case SchedulerProfile::HIGH_CLIENT_OPS:
    profile_str = "high_client_ops";
    break;
  case SchedulerProfile::HIGH_RECOVERY_OPS:
    profile_str = "high_recovery_ops";
    break;
  }
  ss << "  Profile: " << profile_str << std::endl;
  ss << "  Disk Read/Write Bandwidth: " << disk_write_bandwidth << std::endl;
  ss << "  Pool ID: " << pool_id << std::endl;
  ss << "  PG Objects: " << pg_objects << std::endl;
  ss << "  Object Size: " << object_size << std::endl;
  ss << "  Clients: ";
  for (const auto &c : clients)
    ss << c << " ";
  ss << std::endl;
  ss << "  Client Read Queue: " << client_read_queue << std::endl;
  ss << "  Client Write Queue: " << client_write_queue << std::endl;
  ss << "  Engine Config: ";
  for (const auto &s : cfg)
    ss << s << " ";
  ss << std::endl;
  ss << "  Shapes: ";
  for (const auto &s : shapes)
    ss << s << " ";
  ss << std::endl;
  ss << "  Speeds: ";
  for (const auto &s : speeds)
    ss << s << " ";
  ss << std::endl;

  return ss.str();
}

// --- Topology Export ---

void SimContext::serialize_topology(
    const std::string &filename,
    const std::map<std::string, std::vector<std::string>> &host_actors) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    XBT_ERROR("Failed to open topology file: %s", filename.c_str());
    return;
  }

  simgrid::s4u::Engine *e = simgrid::s4u::Engine::get_instance();
  auto *root = e->get_netzone_root();

  // Pre-process hosts to organize by zone
  std::map<simgrid::s4u::NetZone *, std::vector<std::string>> zone_to_hosts;
  for (const auto &[hostname, zone_ptr] : host_zones) {
    zone_to_hosts[zone_ptr].push_back(hostname);
  }

  // Helper to strip parent prefix from name
  auto get_short_name = [](const std::string &full_name,
                           const std::string &parent_name) {
    if (full_name.rfind(parent_name + "_", 0) == 0) {
      return full_name.substr(parent_name.length() + 1);
    }
    return full_name;
  };

  // Recursive Lambda
  std::function<void(simgrid::s4u::NetZone *, int)> dump_zone;
  dump_zone = [&](simgrid::s4u::NetZone *zone, int depth) {
    if (!zone)
      return;
    XBT_INFO("Dumping zone: %s", zone->get_name().c_str());
    std::string indent(depth * 2, ' ');

    std::string parent_name = "";
    if (zone->get_name() != "_world_") {
      if (auto *parent = zone->get_parent()) {
        parent_name = parent->get_name();
      }
    }
    std::string short_name = get_short_name(zone->get_name(), parent_name);

    // Special case for root/star
    if (zone->get_name() == "star" && parent_name == "_world_")
      short_name = "star";

    file << indent << "{\n";
    file << indent << "  \"id\": \"" << zone->get_name() << "\",\n";
    file << indent << "  \"name\": \"" << short_name << "\",\n";
    file << indent << "  \"type\": \"zone\",\n";

    // Try to find an uplink for this zone
    std::string uplink_name = "";
    double bandwidth = 0.0;

    // Convention: zone_name + "_uplink"
    std::string uplink_base = zone->get_name() + "_uplink";
    // Check for simple link or split-duplex _UP
    auto *uplink_link = simgrid::s4u::Link::by_name_or_null(uplink_base);
    if (!uplink_link) {
      uplink_link = simgrid::s4u::Link::by_name_or_null(uplink_base + "_UP");
    }

    if (uplink_link) {
      uplink_name = uplink_base; // Use base name for cleaner topology? Or the
                                 // actual link name found?
      // Use base name if it was split, or the name found? User wants to match
      // metrics. Metrics use "base_UP" and "base_DOWN". If we found "base_UP",
      // we should probably report "base" as the uplink ID and maybe specify it
      // is split duplex? Or just report the "uplink" as a logical object. Let's
      // stick to the base name if possible, or the one found if simple.

      // If we found _UP, implies split duplex.
      // Bandwidth: get_bandwidth() returns bytes/sec. Convert to bits/sec.
      bandwidth = uplink_link->get_bandwidth() * 8.0;
      // Wait, Link::get_latency() returns double.
    }

    file << indent << "  \"uplink\": \"" << uplink_name << "\",\n";
    file << indent << "  \"bandwidth\": " << bandwidth << ",\n";
    file << indent << "  \"latency\": "
         << (uplink_link ? uplink_link->get_latency() : 0.0)
         << ",\n"; // in seconds

    file << indent << "  \"children\": [\n";

    bool first_child = true;

    // 1. Sub-zones
    auto children_vec = zone->get_children();
    for (auto *child : children_vec) {
      if (!first_child)
        file << ",\n";
      dump_zone(child, depth + 1);
      first_child = false;
    }

    // 2. Hosts in this zone
    if (zone_to_hosts.count(zone)) {
      const auto &hosts = zone_to_hosts[zone];
      for (const auto &hostname : hosts) {
        if (!first_child)
          file << ",\n";
        first_child = false;

        std::string host_indent = indent + "  ";
        file << host_indent << "{\n";
        file << host_indent << "  \"id\": \"" << hostname << "\",\n";
        file << host_indent << "  \"name\": \"" << hostname << "\",\n";
        file << host_indent << "  \"type\": \"host\",\n";

        // Host link
        std::string host_link_name = "";
        double h_bw = 0.0;
        double h_lat = 0.0;

        std::string h_link_base = hostname + "_link";
        auto *host_link = simgrid::s4u::Link::by_name_or_null(h_link_base);
        if (!host_link) {
          host_link = simgrid::s4u::Link::by_name_or_null(h_link_base + "_UP");
        }

        if (host_link) {
          host_link_name = h_link_base;
          h_bw = host_link->get_bandwidth() * 8.0;
          h_lat = host_link->get_latency();
        }

        file << host_indent << "  \"uplink\": \"" << host_link_name << "\",\n";
        file << host_indent << "  \"bandwidth\": " << h_bw << ",\n";
        file << host_indent << "  \"latency\": " << h_lat << ",\n";

        // Actors on this host
        file << host_indent << "  \"children\": [\n";
        if (host_actors.count(hostname)) {
          const auto &actors = host_actors.at(hostname);
          for (size_t k = 0; k < actors.size(); ++k) {
            if (k > 0)
              file << ",\n";
            file << host_indent << "    {\n";
            // Actor entry in map is currently "name, disk" or just "name"
            // Let's parse it if needed, or just use it as name.
            // valid JSON string
            file << host_indent << "      \"id\": \"" << hostname << "_" << k
                 << "\",\n"; // Synthetic ID
            file << host_indent << "      \"name\": \"" << actors[k] << "\",\n";
            file << host_indent << "      \"type\": \"actor\"\n";
            file << host_indent << "    }";
          }
        }
        file << "\n" << host_indent << "  ]\n";
        file << host_indent << "}";
      }
    }

    file << "\n" << indent << "  ]\n";
    file << indent << "}";
  };

  dump_zone(root, 0);
  file.close();
  XBT_INFO("Topology exported to %s", filename.c_str());
}