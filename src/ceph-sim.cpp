#include "CLI11.hpp"
#include "ClientActor.hpp"
#include "MonActor.hpp"
#include "OsdActor.hpp"
#include <simgrid/s4u.hpp>

#define READ_BANDWIDTH (120 * 1024 * 1024) // 120 MiB/s
#define WRITE_BANDWIDTH (80 * 1024 * 1024) // 80 MiB/s
XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_sim, "Messages specific for ceph-sim");

#include "MetricMonitor.hpp"
#include "SimContext.hpp"

// --- Helper Functions and Topology Builders have been moved to SimContext.cpp
// ---

// Recursive function to print the NetZone tree
std::string get_tree_str(
    simgrid::s4u::NetZone *zone, std::string prefix,
    const std::map<std::string, std::vector<std::string>> &host_actors,
    const std::map<simgrid::s4u::NetZone *, std::vector<simgrid::s4u::Host *>>
        &zone_hosts) {
  std::stringstream result_ss;
  result_ss << prefix << zone->get_name() << "\n";

  // Check for uplink
  auto *uplink =
      simgrid::s4u::Link::by_name_or_null(zone->get_name() + "_uplink");
  if (uplink) {
    result_ss << prefix << "  Uplink: "
              << uplink->get_bandwidth() * 8.0 / 1000 / 1000 / 1000
              << " Gbps\n";
  }

  // Print hosts in this zone
  if (zone_hosts.count(zone)) {
    for (auto *host : zone_hosts.at(zone)) {
      std::string host_name = host->get_name();
      result_ss << prefix << "  Host: " << host_name << "\n";

      // Check for host uplink
      auto *host_link =
          simgrid::s4u::Link::by_name_or_null(host_name + "_link");
      if (host_link) {
        result_ss << prefix << "    Uplink: "
                  << host_link->get_bandwidth() * 8.0 / 1000 / 1000 / 1000
                  << " Gbps\n";
      }

      if (host_actors.count(host_name)) {
        for (const auto &actor_info : host_actors.at(host_name)) {
          result_ss << prefix << "    " << actor_info << "\n";
        }
      }
    }
  }

  // Recurse into children zones
  for (auto *child : zone->get_children()) {
    result_ss << get_tree_str(child, prefix + "  ", host_actors, zone_hosts);
  }
  return result_ss.str();
}

void seal_netzone(simgrid::s4u::NetZone *zone) {
  for (auto *child : zone->get_children()) {
    seal_netzone(child);
  }
  zone->seal();
}

int main(int argc, char *argv[]) {

  // CLI11 argument parsing
  CLI::App app{"Ceph simulator"};

  // footer
  std::string footer_text = R"(
  <SHAPE>:
  The shape/speed/weight arguments use a colon-separated hierarchy:
     Rack : Host : OSD
  
  Expansion Rules:
  - Uniform:   Use a single value (e.g., '25') to apply to all items.
  - Pattern:   Use a comma list (e.g., '10,20') to cycle through values.
  - Specific:  Use '@' to cycle based on parent (e.g., '@10,100' alternates racks).
  - Inherit:   Leave empty (e.g., '::') to inherit from parent level.

  Hierarchy Inheritance:
  - shape:     Absolute hierarchy
  - speed:     Each level has unique links
  - weight:    Only OSDs have weights

  Examples:
  --dc-shape 4:10:24        (4 Racks, 10 Hosts, 24 Disks)
  --dc-speed 200:200:25       (200G Uplink, 200G Rack, 25G Host)
  --dc-weight 14::            (Mix of 10G and 100G Disks) // todo I don't need this
  --dc-weight :10,14:         (Mix of 10G and 100G Disks)
  )";
  app.footer(footer_text);
  app.set_help_flag("--help2", "Print this help message and exit");

  SimContext ctx;

  // --dc-shape
  auto *dc_shape_opt =
      app.add_option("--dc-shape", ctx.shapes, "Add a data center");
  dc_shape_opt->required()->expected(1)->type_name("SHAPE");

  // --dc-speed
  auto *dc_speed_opt = app.add_option("--dc-speed", ctx.speeds,
                                      "Configure data center link speeds");
  dc_speed_opt->required()->expected(1)->type_name("SHAPE");

  // --dc-weight
  std::vector<std::string> dc_weights;
  auto *dc_weight_opt = app.add_option("--dc-weight", ctx.weights,
                                       "Configure data center link speeds");
  dc_weight_opt->required()->expected(1)->type_name("SHAPE");

  // --pgdump
  std::vector<std::string> pgdump_files;
  auto *pgdump_opt = app.add_option("--pgdump", pgdump_files,
                                    "Path to pgdump output files (min 2)");
  pgdump_opt->required()->check(CLI::ExistingFile)->expected(2);

  // --pg-objects
  auto *pg_objects_opt = app.add_option("--pg-objects", ctx.pg_objects,
                                        "Number of objects per PG");
  pg_objects_opt->required()->expected(1)->type_name("INT");

  // --object-size
  auto *object_size_opt =
      app.add_option("--object-size", ctx.object_size, "Size of objects");
  object_size_opt->required()->expected(1)->type_name("INT");

  // --pool-id
  auto *pool_id_opt = app.add_option("--pool-id", ctx.pool_id, "Pool ID");
  pool_id_opt->default_val(1)->expected(1)->type_name("INT");

  // --disk-read-bandwidth
  auto *disk_read_bandwidth_opt = app.add_option(
      "--disk-read-bandwidth", ctx.disk_read_bandwidth, "Disk read bandwidth");
  disk_read_bandwidth_opt->default_val(READ_BANDWIDTH)
      ->expected(1)
      ->type_name("INT");

  // --disk-write-bandwidth
  auto *disk_write_bandwidth_opt =
      app.add_option("--disk-write-bandwidth", ctx.disk_write_bandwidth,
                     "Disk write bandwidth");
  disk_write_bandwidth_opt->default_val(WRITE_BANDWIDTH)
      ->expected(1)
      ->type_name("INT");

  // --cfg (passthrough to simgrid engine)
  auto *cfg_opt = app.add_option("--cfg", ctx.cfg, "Simgrid configuration");
  cfg_opt->type_name("CFG");

  // --help (passthrough to simgrid engine)
  app.add_flag("--help", "Print this help message and exit")
      ->default_val(false);

  // --help-tracing
  app.add_flag("--help-tracing", "Print this help message and exit")
      ->default_val(false);

  // Do the parsing
  CLI11_PARSE(app, argc, argv);

  // Create the simgrid engine
  simgrid::s4u::Engine e(&argc, argv);
  // todo: Engine has set_config(string, string)

  // initialize pg map
  PGMap *pgmap = new PGMap(ctx.pool_id, pgdump_files[0], ctx.object_size,
                           ctx.pg_objects); // acting
  PGMap up_pgmap(ctx.pool_id, pgdump_files[1], ctx.object_size,
                 ctx.pg_objects); // up
  for (size_t i = 0; i < pgmap->size(); i++) {
    std::vector<int> up_members = up_pgmap.get_pg(i)->get_up_ids();
    pgmap->get_pg(i)->init_up(up_members);
  }
  pgmap->_update_primary_osd_to_pg_index();
  XBT_INFO("%s", pgmap->to_string().c_str());

  // Build the platform
  auto *world_star = e.get_netzone_root()->add_netzone_star("star");

  // Iterate over defined shapes and build data centers
  for (size_t i = 0; i < ctx.shapes.size(); ++i) {
    build_dc(ctx, world_star, i);
  }

  // Deploy the application
  std::map<std::string, std::vector<std::string>> host_actors;
  std::map<simgrid::s4u::NetZone *, std::vector<simgrid::s4u::Host *>>
      zone_hosts;

  // deploy the osds
  for (auto host : e.get_all_hosts()) {
    // Populate zone map
    if (ctx.host_zones.count(host->get_name())) {
      zone_hosts[ctx.host_zones[host->get_name()]].push_back(host);
    } else {
      xbt_die("Host %s not found in host_zones", host->get_name().c_str());
    }

    std::vector<sg4::Disk *> disks = host->get_disks();
    for (auto disk : disks) {
      int osd_id = ctx.global_osd_id++;
      std::string osd_name = "osd." + std::to_string(osd_id);

      // Store info for tree view
      host_actors[host->get_name()].push_back(osd_name + ", " +
                                              disk->get_name());
      e.add_actor(osd_name, host, [pgmap, osd_id, disk]() {
        Osd osd(pgmap, osd_id, disk->get_name());
        osd();
      });
    }
  }

  // deploy a mon (to first rack on a new host)
  auto nz = ctx.host_zones["host-0"];
  sg4::Host *mon = nz->add_host("mon", "100Gf");

  // Create Uplink for Mon (25Gbps default)
  auto mon_link =
      nz->add_split_duplex_link("mon_link", "25Gbps")->set_latency("5us");
  nz->add_route(mon, nullptr, {{mon_link, sg4::LinkInRoute::Direction::UP}},
                true);

  // Add keys to maps
  zone_hosts[nz].push_back(mon);
  host_actors["mon"].push_back("mon");

  e.add_actor("mon", mon, [pgmap]() {
    Mon mon(pgmap);
    mon();
  });

  // seal racks and data centers
  world_star->seal();
  // seal_netzone(world_star);

  // Print the deployment tree
  std::string tree_str =
      get_tree_str(world_star, "  ", host_actors, zone_hosts);
  XBT_INFO("Deployment Tree:\n%s", tree_str.c_str());

  XBT_INFO("Simulation context:\n%s", ctx.to_string().c_str());

  // Export Topology
  ctx.serialize_topology("topology.json");

  // Start Metric Monitor
  e.add_actor("metric_monitor", e.get_all_hosts()[0], []() {
    sg4::Actor::self()->daemonize();
    MetricMonitor monitor(1.0, "metrics.csv");
    monitor();
  });

  auto start_time = std::chrono::steady_clock::now();
  /* Run the simulation */
  e.run();
  auto end_time = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  XBT_INFO("Simulation wall time: %lld s", duration.count());

  XBT_INFO("Simulation is over");

  return 0;
}