#include "CLI11.hpp"
#include "ClientActor.hpp"
#include "MonActor.hpp"
#include "OsdActor.hpp"
#include <filesystem>
#include <simgrid/s4u.hpp>
namespace fs = std::filesystem;

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

  // Pre-scan argv for --help before anything else, because SimGrid's
  // Engine constructor will print its own help and exit() on --help,
  // leaving no chance for CLI11 help to be shown afterward.
  bool help_requested = false;
  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") {
      help_requested = true;
      break;
    }
  }

  // CLI11 argument parsing
  // Disable built-in --help so SimGrid can handle it
  CLI::App app{"Ceph simulator"};
  app.set_help_flag();

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

  Examples:
  --dc-shape 4:10:24        (4 Racks, 10 Hosts, 24 Disks)
  --dc-shape 2:3,5:24       (first rack, 3 Hosts, second rack, 5 Hosts, 24 Disks)
  --dc-speed 200:25:10      (200G Uplink, 25G Rack, 10G Host)
  )";
  app.footer(footer_text);

  SimContext ctx;

  // --dc-shape
  auto *dc_shape_opt =
      app.add_option("--dc-shape", ctx.shapes, "Add a data center");
  dc_shape_opt->required()->type_name("SHAPE")->delimiter(',');

  // --dc-speed
  auto *dc_speed_opt = app.add_option("--dc-speed", ctx.speeds,
                                      "Configure data center link speeds");
  dc_speed_opt->required()->type_name("SHAPE")->delimiter(',');

  // --pgdump
  std::vector<std::string> pgdump_files;
  auto *pgdump_opt = app.add_option("--pgdump", pgdump_files,
                                    "Path to pgdump output files (min 2)");
  pgdump_opt->required()->check(CLI::ExistingFile)->expected(2);

  // --start-up-delay
  auto *start_up_delay_opt =
      app.add_option("--start-up-delay", ctx.start_up_delay, "Start up delay");
  start_up_delay_opt->default_val(0)->expected(1)->type_name("INT");

  // --shut-down-delay
  auto *shut_down_delay_opt = app.add_option(
      "--shut-down-delay", ctx.shut_down_delay, "Shut down delay");
  shut_down_delay_opt->default_val(0)->expected(1)->type_name("INT");

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

  // --disk-write-bandwidth
  auto *disk_write_bandwidth_opt =
      app.add_option("--disk-write-bandwidth", ctx.disk_write_bandwidth,
                     "Disk write bandwidth");
  disk_write_bandwidth_opt->default_val(WRITE_BANDWIDTH)
      ->expected(1)
      ->type_name("INT");

  // --iops
  auto *iops_opt = app.add_option("--iops", ctx.iops, "IOPS");
  iops_opt->default_val(100)->expected(1)->type_name("INT");

  // --profile
  std::map<std::string, SchedulerProfile> profile_map{
      {"balanced", SchedulerProfile::BALANCED},
      {"high_client_ops", SchedulerProfile::HIGH_CLIENT_OPS},
      {"high_recovery_ops", SchedulerProfile::HIGH_RECOVERY_OPS}};

  auto *profile_opt = app.add_option("--profile", ctx.profile, "Profile");
  profile_opt->transform(CLI::CheckedTransformer(profile_map, CLI::ignore_case))
      ->default_val(SchedulerProfile::BALANCED)
      ->expected(1)
      ->type_name("PROFILE");

  // --output-dir
  auto *output_dir_opt = app.add_option("--output-dir", ctx.output_dir,
                                        "Output directory for artifacts");
  output_dir_opt->required()->type_name("DIR");

  // --dc-clients
  auto *clients_opt = app.add_option("--dc-clients", ctx.clients,
                                     "Number of clients (repeatable)");

  // --client-read-queue-depth
  auto *client_read_queue_opt =
      app.add_option("--client-read-queue-depth", ctx.client_read_queue,
                     "Client read queue depth");
  client_read_queue_opt->default_val(1);

  // --client-write-queue-depth
  auto *client_write_queue_opt =
      app.add_option("--client-write-queue-depth", ctx.client_write_queue,
                     "Client write queue depth");
  client_write_queue_opt->default_val(1);

  // --client-op-size
  auto *client_op_size_opt = app.add_option(
      "--client-op-size", ctx.client_op_size, "Client operation size");
  client_op_size_opt->default_val(4096);

  // If --help was requested, print CLI11 help first (before Engine init,
  // because Engine will print SimGrid help and exit on --help).
  if (help_requested) {
    std::cout << app.help() << std::endl;
    std::cout << "--- SimGrid Engine Options ---" << std::endl;
  }

  // Initialize SimGrid engine. Handles --cfg, --help, and --help-tracing.
  // On --help, SimGrid prints its own help and exits (after CLI11 help above).
  // Otherwise, it strips --cfg from argv so CLI11 only sees app-specific flags.
  simgrid::s4u::Engine e(&argc, argv);

  // Do the parsing
  CLI11_PARSE(app, argc, argv);
  XBT_INFO("Context:\n%s", ctx.to_string().c_str());

  // Create output directory
  try {
    if (!ctx.output_dir.empty() && !fs::exists(ctx.output_dir)) {
      fs::create_directories(ctx.output_dir);
      XBT_INFO("Created output directory: %s", ctx.output_dir.c_str());
    }
  } catch (const std::exception &e) {
    XBT_ERROR("Failed to create output directory: %s", e.what());
    return 1;
  }

  // initialize pg map
  PGMap *pgmap = new PGMap(ctx.pool_id, pgdump_files[0], ctx.object_size,
                           ctx.pg_objects); // acting
  PGMap up_pgmap(ctx.pool_id, pgdump_files[1], ctx.object_size,
                 ctx.pg_objects); // up
  for (size_t i = 0; i < pgmap->size(); i++) {
    std::vector<int> up_members = up_pgmap.get_pg(i)->get_up_ids();
    pgmap->get_pg(i)->update_up(up_members);
  }
  pgmap->update_primary_osd_to_pg_index();
  XBT_INFO("%s", pgmap->to_string().c_str());
  XBT_INFO("Initial map's reverse index\n%s",
           pgmap->primary_osds_to_pgs_string().c_str());

  // Build the platform
  auto *world_star = e.get_netzone_root()->add_netzone_star("star");
  auto *world_gw = world_star->add_router("world_gw");
  world_star->set_gateway(world_gw);

  // Iterate over defined shapes and build data centers
  for (size_t i = 0; i < ctx.shapes.size(); ++i) {
    build_dc(ctx, world_star, i);
  }

  // Deploy the application
  std::map<std::string, std::vector<std::string>> host_actors;
  std::map<simgrid::s4u::NetZone *, std::vector<simgrid::s4u::Host *>>
      zone_hosts;

  // Calculate total OSDs for validation
  int total_osds = 0;
  for (auto host : e.get_all_hosts()) {
    total_osds += host->get_disks().size();
  }
  xbt_assert(pgmap->get_max_osd_id() < total_osds,
             "PGMap expects OSD %d, but only %d OSDs are configured in the "
             "simulation (meaning max_id %d)! Deadlock inevitable.",
             pgmap->get_max_osd_id(), total_osds, total_osds - 1);

  // deploy the osds
  for (const auto &host_name : ctx.ordered_host_names) {
    auto *host = e.host_by_name(host_name);
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
      e.add_actor(osd_name, host, [pgmap, osd_id, disk, ctx]() {
        Osd osd(pgmap, osd_id, disk->get_name(), ctx.iops, ctx.profile);
        osd();
      });
    }
  }

  // add client
  std::vector<std::string> client_names;
  // Configure ClientMetrics
  std::string client_metrics_path =
      (fs::path(ctx.output_dir) / "client_metrics.csv").string();
  Client::set_metrics_output(client_metrics_path);

  // Enable aggregated metrics (memory-efficient alternative to per-op CSV)
  Client::set_aggregate_output(ctx.output_dir);

  if (!ctx.clients.empty()) {
    int client_global_counter = 1;
    for (size_t dc_idx = 0; dc_idx < ctx.clients.size(); ++dc_idx) {
      int count = ctx.clients[dc_idx];
      std::string dc_name = "dc-" + std::to_string(dc_idx);
      // Attempt to find the DC zone
      simgrid::s4u::NetZone *dc_zone = nullptr;
      for (auto *child : world_star->get_children()) {
        if (child->get_name() == dc_name) {
          dc_zone = child;
          break;
        }
      }
      // Fallback to world_star if DC doesn't exist (e.g. more client args than
      // DCs)
      if (!dc_zone) {
        XBT_WARN(
            "DC %s not found for client group %lu. Attaching to star zone.",
            dc_name.c_str(), dc_idx);
      }

      // Create or find client rack if we have a valid DC zone
      simgrid::s4u::NetZone *client_rack_zone = nullptr;
      if (dc_zone) {
        std::string client_rack_name = dc_name + "_client_rack";
        // Check if it exists
        for (auto *child : dc_zone->get_children()) {
          if (child->get_name() == client_rack_name) {
            client_rack_zone = child;
            break;
          }
        }

        // Create if not exists
        if (!client_rack_zone) {
          client_rack_zone = dc_zone->add_netzone_star(client_rack_name);
          auto *router =
              client_rack_zone->add_router(client_rack_name + "_router");
          client_rack_zone->set_gateway(router);

          auto *uplink = dc_zone
                             ->add_split_duplex_link(
                                 client_rack_name + "_uplink", "100Gbps")
                             ->set_latency("0us");
          dc_zone->add_route(client_rack_zone, nullptr,
                             {{uplink, sg4::LinkInRoute::Direction::UP}}, true);
        }
      }

      auto *target_zone = client_rack_zone ? client_rack_zone
                                           : (dc_zone ? dc_zone : world_star);

      for (int i = 0; i < count; ++i) {
        std::string client_name =
            "client." + std::to_string(client_global_counter);
        client_names.push_back(client_name);

        sg4::Host *client = target_zone->add_host(client_name, "100Gf");
        auto client_link =
            target_zone->add_split_duplex_link(client_name + "_link", "25Gbps")
                ->set_latency("5us");
        target_zone->add_route(client, nullptr,
                               {{client_link, sg4::LinkInRoute::Direction::UP}},
                               true);

        int client_id = -client_global_counter;
        e.add_actor(client_name, client, [pgmap, ctx, client_id]() {
          Client client(pgmap, client_id, ctx.client_read_queue,
                        ctx.client_write_queue, ctx.client_op_size);
          client();
        });
        zone_hosts[target_zone].push_back(client);
        host_actors[client_name].push_back(client_name);
        // Register for topology export
        ctx.host_zones[client_name] = target_zone;

        client_global_counter++;
      }
    }
  }

  // deploy a mon (to first rack on a new host)
  auto nz = ctx.host_zones["host-0"];
  sg4::Host *mon = nz->add_host("mon", "100Gf");
  auto mon_link =
      nz->add_split_duplex_link("mon_link", "25Gbps")->set_latency("0us");
  nz->add_route(mon, nullptr, {{mon_link, sg4::LinkInRoute::Direction::UP}},
                true);
  zone_hosts[nz].push_back(mon);
  host_actors["mon"].push_back("mon");
  // Register for topology export
  ctx.host_zones["mon"] = nz;

  // Configure MonMetrics
  std::string mon_metrics_path =
      (fs::path(ctx.output_dir) / "mon_metrics.csv").string();
  Mon::set_metrics_output(mon_metrics_path);

  e.add_actor("mon", mon, [pgmap, client_names, ctx]() {
    Mon mon(pgmap, client_names, ctx.start_up_delay, ctx.shut_down_delay);
    mon();
  });

  // seal racks and data centers
  world_star->seal();

  // Print the deployment tree
  std::string tree_str =
      get_tree_str(world_star, "  ", host_actors, zone_hosts);
  XBT_INFO("Deployment Tree:\n%s", tree_str.c_str());

  // Export Topology
  XBT_INFO("Exporting topology...");
  std::string topo_path = (fs::path(ctx.output_dir) / "topology.json").string();
  ctx.serialize_topology(topo_path, host_actors);

  // Start Metric Monitor
  // Start Metric Monitor
  std::string metrics_path =
      (fs::path(ctx.output_dir) / "net_metrics.csv").string();

  e.add_actor("metric_monitor", e.get_all_hosts()[0], [metrics_path]() {
    sg4::Actor::self()->daemonize();
    MetricMonitor monitor(1.0, metrics_path);
    monitor();
  });

  auto start_time = std::chrono::steady_clock::now();
  /* Run the simulation */
  e.run();
  auto end_time = std::chrono::steady_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  XBT_INFO("Simulation wall time: %ld s", duration.count());

  // Write aggregated client metrics (T-Digest and throughput buckets)
  Client::write_aggregated_metrics();

  delete pgmap;
  return 0;
}