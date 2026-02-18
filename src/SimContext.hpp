#pragma once

#include "CLI11.hpp"
#include "OsdActor.hpp"
#include <simgrid/s4u.hpp>
#include <string>
#include <vector>

// --- Helper Functions ---
std::vector<std::string> split_str(const std::string &s, char delim);

double resolve_val(std::string spec, int self_id, int parent_id,
                   double fallback);

// --- Validator ---
struct ShapeValidator : public CLI::Validator {
  ShapeValidator();
};

// --- Context ---
struct SimContext {
  int pool_id = 0;
  int global_rack_id = 0;
  int global_host_id = 0;
  int global_osd_id = 0;

  long start_up_delay = 0;
  long shut_down_delay = 0;

  int disk_read_bandwidth = 0;
  int disk_write_bandwidth = 0;
  int iops = 0;
  SchedulerProfile profile = SchedulerProfile::BALANCED;

  std::vector<int> clients;
  int client_read_queue = 1;
  int client_write_queue = 1;
  int client_op_size = 4096;
  double client_read_bandwidth = 0;  // 0 = unlimited (use link default)
  double client_write_bandwidth = 0; // 0 = unlimited (use link default)

  std::vector<std::string> shapes;
  std::vector<std::string> speeds;

  int pg_objects = 1000;
  int object_size = 4 * 1024 * 1024;

  std::map<std::string, simgrid::s4u::NetZone *> host_zones;
  std::vector<std::string> ordered_host_names;
  std::string output_dir = ".";

  std::string get_hierarchy_spec(const std::vector<std::string> &vec,
                                 int dc_idx, int level);
  std::string to_string();

  // Topology Export
  void serialize_topology(
      const std::string &filename,
      const std::map<std::string, std::vector<std::string>> &host_actors);
};

// --- Top-Level Builder ---
void build_dc(SimContext &ctx, simgrid::s4u::NetZone *world, int dc_config_idx);
