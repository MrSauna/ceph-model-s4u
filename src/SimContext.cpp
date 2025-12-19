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
                       int rack_idx, int host_idx, double parent_speed,
                       double parent_weight) {

  std::string hostname = "host-" + std::to_string(ctx.global_host_id++);

  // Create Host
  auto *host = rack_zone->add_host(hostname, "100Gf");
  ctx.host_zones[hostname] = rack_zone;

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

  std::string osd_w_spec = ctx.get_hierarchy_spec(ctx.weights, dc_idx, 2);

  // Resolve Host Weight to use as parent for OSDs
  // Note: We use the same index (2) for weight inheritance base here
  double host_weight =
      resolve_val(osd_w_spec, host_idx, rack_idx, parent_weight);

  for (int o = 0; o < osd_count; ++o) {
    // Resolve final OSD weight
    double osd_weight = resolve_val(osd_w_spec, o, host_idx, host_weight);
    (void)osd_weight; // Unused for now, but resolved
    host->add_disk("disk" + std::to_string(o), ctx.disk_read_bandwidth,
                   ctx.disk_write_bandwidth);
  }
}

// Level 2: Rack (Index 1)
static void build_rack(SimContext &ctx, sg4::NetZone *dc_zone, int dc_idx,
                       int rack_idx, double parent_speed,
                       double parent_weight) {

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
  std::string weight_spec = ctx.get_hierarchy_spec(ctx.weights, dc_idx, 1);
  double rack_weight =
      resolve_val(weight_spec, rack_idx, dc_idx, parent_weight);

  for (int h = 0; h < host_count; ++h) {
    build_host(ctx, rack_zone, dc_idx, rack_idx, h, speed_val, rack_weight);
  }
}

// --- Top-Level Builder ---

// Level 1: Data Center (Index 0)
void build_dc(SimContext &ctx, sg4::NetZone *world, int dc_config_idx) {
  std::string dc_name = "dc-" + std::to_string(dc_config_idx);
  auto *dc_zone = world->add_netzone_star(dc_name);

  // Router: DC Switch
  // auto *dc_router = dc_zone->add_router(dc_name + "_router");
  // dc_zone->set_gateway(dc_router);

  // Resolve DC Uplink Speed (Index 0: DC Uplink)
  std::string speed_spec = ctx.get_hierarchy_spec(ctx.speeds, dc_config_idx, 0);
  double speed_val = resolve_val(speed_spec, dc_config_idx, 0, 1000.0);
  std::string bw_str = std::to_string(speed_val) + "Gbps";

  auto *uplink =
      dc_zone->add_link(dc_name + "_uplink", bw_str)->set_latency("100us");

  // Route: DC Switch -> Outside (DC Boundary)
  // DISABLED: dc_zone->add_route(dc_router, nullptr, {{uplink,
  // sg4::LinkInRoute::Direction::UP}}, true); Route: DC Switch -> World (in
  // World Zone) DISABLED: world->add_route(dc_router, nullptr, {}, true);

  // Rack Count (Index 0)
  std::string count_spec = ctx.get_hierarchy_spec(ctx.shapes, dc_config_idx, 0);
  int rack_count = std::stoi(count_spec.empty() ? "0" : count_spec);

  // Weight (Index 0)
  std::string weight_spec =
      ctx.get_hierarchy_spec(ctx.weights, dc_config_idx, 0);
  double dc_weight = resolve_val(weight_spec, dc_config_idx, 0, 4.0);

  for (int r = 0; r < rack_count; ++r) {
    build_rack(ctx, dc_zone, dc_config_idx, r, speed_val, dc_weight);
  }
}

std::string SimContext::to_string() {
  std::stringstream ss;
  ss << "  DCs: " << shapes.size() << std::endl;
  ss << "  Racks: " << global_rack_id << std::endl;
  ss << "  Hosts: " << global_host_id << std::endl;
  ss << "  OSDs: " << global_osd_id << std::endl;
  ss << "  Disk Read Bandwidth: " << disk_read_bandwidth << std::endl;
  ss << "  Disk Write Bandwidth: " << disk_write_bandwidth << std::endl;
  ss << "  Pool ID: " << pool_id << std::endl;
  ss << "  PG Objects: " << pg_objects << std::endl;
  ss << "  Object Size: " << object_size << std::endl;
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
  ss << "  Weights: ";
  for (const auto &s : weights)
    ss << s << " ";
  ss << std::endl;
  return ss.str();
}

// --- Topology Export ---

void SimContext::serialize_topology(const std::string &filename) {
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

  // Recursive Lambda
  std::function<void(simgrid::s4u::NetZone *, int)> dump_zone;
  dump_zone = [&](simgrid::s4u::NetZone *zone, int depth) {
    std::string indent(depth * 2, ' ');
    file << indent << "{\n";
    file << indent << "  \"name\": \"" << zone->get_name() << "\",\n";

    file << indent << "  \"hosts\": [";
    if (zone_to_hosts.count(zone)) {
      const auto &hosts = zone_to_hosts[zone];
      for (size_t i = 0; i < hosts.size(); ++i) {
        if (i > 0)
          file << ", ";
        file << "\"" << hosts[i] << "\"";
      }
    }
    file << "],\n";

    file << indent << "  \"children\": [\n";
    auto children_vec = zone->get_children();
    bool first = true;
    for (auto *child : children_vec) {
      if (!first)
        file << ",\n";
      dump_zone(child, depth + 1);
      first = false;
    }
    file << "\n" << indent << "  ]\n";
    file << indent << "}";
  };

  dump_zone(root, 0);
  file.close();
  XBT_INFO("Topology exported to %s", filename.c_str());
}