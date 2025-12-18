#include "CephCommon.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex> // for std::scoped_lock
#include <sstream>
#include <string>

namespace fs = std::filesystem;

XBT_LOG_NEW_DEFAULT_CATEGORY(s4u_ceph_common,
                             "Messages specific for CephCommon");

int sender_str_to_int(std::string sender) {
  if (sender.rfind("osd.", 0) == 0) {
    return std::stoi(sender.substr(4));
  } else if (sender.rfind("client.", 0) == 0) {
    return -std::stoi(sender.substr(7));
  } else {
    xbt_die("Unknown sender: %s", sender.c_str());
  }
}

std::istream &operator>>(std::istream &is, std::vector<int> &v) {
  int osd_id;
  char ch;
  std::vector<int> members;

  // "(["
  if (!(is >> ch && ch == '(' && is >> ch && ch == '[')) {
    is.setstate(std::ios_base::failbit);
    return is;
  }

  // "1,20,30"
  is >> osd_id;
  members.push_back(osd_id);
  while (is.peek() == ',') {
    is >> ch;
    is >> osd_id;
    members.push_back(osd_id);
  }

  // "], p"
  if (!(is >> ch && ch == ']' && is >> ch && ch == ',' && is >> ch &&
        ch == 'p')) {
    is.setstate(std::ios_base::failbit);
    return is;
  }

  // "1"
  is >> osd_id;
  xbt_assert(osd_id == members.at(0),
             "Parsed primary %u does not match first member %u", osd_id,
             members.at(0));

  // "], p"
  if (!(is >> ch && ch == ')')) {
    is.setstate(std::ios_base::failbit);
    return is;
  }

  // only on success do we touch the set reference
  if (is)
    v = std::move(members);
  return is;
}

// PGShard
int PGShard::get_pg_id() const { return pg->get_id(); }

int PGShard::get_osd_id() const { return osd_id; }

bool PGShard::is_acting() {

  auto acting_members = pg->get_acting().members;
  return std::find(acting_members.begin(), acting_members.end(), this) !=
         acting_members.end();
}

// PG
PG::PG(std::string line, size_t object_size, size_t pg_objects)
    : pg_objects(pg_objects), object_size(object_size) {
  std::stringstream ss(line);
  // eg.
  // 3.18 raw ([123,73,98], p123) up ([123,73,98], p123) acting ([123,73,98],
  // p123)
  std::string temp;
  std::vector<int> members;

  // "3.18"
  char del = '.';
  getline(ss, temp, '.');
  ss >> temp;
  id = std::stoul(temp, nullptr, 16);

  // "raw"
  ss >> temp;
  // "([123,73,98], p123)""
  ss >> members; // do nothing with raw

  // "up ([123,73,98], p123)"
  ss >> temp;
  ss >> members;
  init_up(members);

  // "acting ([123,73,98], p123)""
  ss >> temp;
  ss >> members;
  init_acting(members);
}

void PG::init_up(const std::vector<int> &v) {
  xbt_assert(v.size() >= 1);

  PGShardSet new_up;
  // if acting[i] == v[i], we reuse the shard
  for (int i = 0; i < v.size(); ++i) {
    if (i < acting.members.size() && acting.members[i]->get_osd_id() == v[i]) {
      new_up.members.push_back(acting.members[i]);
    } else {
      auto shard = std::make_unique<PGShard>(this, v[i], i);
      new_up.members.push_back(shard.get());
      shards.insert(std::move(shard));
    }
  }
  up = new_up;
  prune_shards();
}

void PG::update_up(const std::vector<int> &v) {
  const std::scoped_lock lock(*mutex_);
  init_up(v);
}

void PG::init_acting(const std::vector<int> &v) {
  xbt_assert(v.size() >= 1);
  PGShardSet new_acting;
  // if up[i] == v[i], we reuse the shard
  for (int i = 0; i < v.size(); ++i) {
    if (i < up.members.size() && up.members[i]->get_osd_id() == v[i]) {
      new_acting.members.push_back(up.members[i]);
    } else {
      auto shard = std::make_unique<PGShard>(this, v[i], i);
      new_acting.members.push_back(shard.get());
      shards.insert(std::move(shard));
    }
  }
  acting = new_acting;
  prune_shards();
}

void PG::update_acting(const std::vector<int> &v) {
  const std::scoped_lock lock(*mutex_);
  init_acting(v);
}

void PG::prune_shards() {
  for (auto it = shards.begin(); it != shards.end();) {
    PGShard *shard = it->get();
    if (!up.contains_shard(shard) && !acting.contains_shard(shard)) {
      it = shards.erase(it);
    } else {
      ++it;
    }
  }
}

std::string PG::to_string() const {
  std::ostringstream ss;
  ss << "PG " << id << "(" << std::hex << id << ") ";
  ss << "up " << up.to_string() << " ";
  ss << "acting " << acting.to_string();
  return ss.str();
}

// PGMap
PGMap::PGMap(int pool_id, std::string path, size_t object_size,
             size_t pg_objects)
    : object_size(object_size), pg_objects(pg_objects) {
  // parse path + assert
  xbt_assert(fs::is_regular_file(path));

  // open
  std::ifstream file(path);
  std::string line;
  int current_pool = 0;
  int pg_num = 0;
  int pgs_to_parse = 0;
  std::string discard;

  // parse
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string first_word;
    ss >> first_word;
    if (first_word == "pool") {
      ss >> current_pool;
      if (current_pool == pool_id) {
        ss >> discard >> pg_num;
        pgs_to_parse = pg_num;
      }

    } else if (current_pool == pool_id && pgs_to_parse > 0) {
      pgs.push_back(std::make_unique<PG>(line, object_size, pg_objects));
      pgs_to_parse--;
    }
  }
}

const std::vector<PGShardSet> PGMap::get_up() const {
  const std::scoped_lock lock(*mutex_);
  std::vector<PGShardSet> up_set;
  for (const auto &pg_ptr : pgs) {
    up_set.push_back(pg_ptr.get()->get_up());
  }
  return up_set;
}

const std::vector<PGShardSet> PGMap::get_acting() const {
  const std::scoped_lock lock(*mutex_);
  std::vector<PGShardSet> acting_set;
  for (const auto &pg_ptr : pgs) {
    acting_set.push_back(pg_ptr->get_acting());
  }
  return acting_set;
}

void PGMap::_update_primary_osd_to_pg_index() {
  primary_osd_to_pg_index.clear();

  // find max osd id
  max_osd_id = 0;
  for (const auto &pg_ptr : pgs) {
    max_osd_id = std::max(max_osd_id, pg_ptr->get_max_osd_id());
  }

  // ensure every osd is in the index
  for (int o = 0; o <= max_osd_id; o++) {
    primary_osd_to_pg_index[o] = std::set<PG *>();
  }

  // populate index
  for (const auto &pg_ptr : pgs) {
    int primary_osd = pg_ptr->get_acting_ids().at(0); // first is primary
    primary_osd_to_pg_index[primary_osd].insert(pg_ptr.get());
  }
}

void PGMap::update_primary_osd_to_pg_index() {
  const std::scoped_lock lock(*mutex_);
  _update_primary_osd_to_pg_index();
}

std::vector<int> PGMap::get_osds() const {
  const std::scoped_lock lock(*mutex_);
  std::vector<int> osds;
  for (const auto &pair : primary_osd_to_pg_index) {
    osds.push_back(pair.first);
  }

  return osds;
}

bool PGMap::needs_backfill() const {
  const std::scoped_lock lock(*mutex_);
  for (const auto &pg_ptr : pgs) {
    if (pg_ptr->needs_backfill()) {
      return true;
    }
  }
  return false;
}

size_t PGMap::size() const {
  // just assume it doesn't change  const std::scoped_lock lock(*mutex_);
  return pgs.size();
}

size_t PGMap::get_object_size() const { return object_size; }

size_t PGMap::get_objects_per_pg() const { return pg_objects; }

bool PG::needs_backfill() const {
  const std::scoped_lock lock(*mutex_);
  return up.members != acting.members;
}

void PG::on_object_recovered() {
  const std::scoped_lock lock(*mutex_);
  objects_recovered++;
  // lazy simplified logic: if we recovered enough, we are done
  if (objects_recovered >= pg_objects) {
    acting = up;
    prune_shards();
    XBT_INFO("PG %d backfill complete (recovered %d objects)", id,
             objects_recovered);
  }
}

bool PG::schedule_recovery() {
  const std::scoped_lock lock(*mutex_);
  if (objects_recoveries_scheduled >= pg_objects) {
    return false;
  }
  objects_recoveries_scheduled++;
  return true;
}

std::string PGMap::primary_osds_to_pgs_string() const {
  const std::scoped_lock lock(*mutex_);
  std::ostringstream ss;
  for (const auto &pair : primary_osd_to_pg_index) {
    const char *sep = "";
    ss << "osd." << pair.first << " -> pg ";
    for (const auto &pg : pair.second) {
      ss << sep << pg->get_id() << "(" << std::hex << pg->get_id() << ")";
      sep = ", ";
    }
    ss << "\n";
  }
  return ss.str();
}

std::string PGMap::to_string() const {
  std::ostringstream ss;
  ss << "PGMap (Size " << pgs.size() << "):\n";
  for (const auto &pg : pgs) {
    ss << "  " << pg->to_string() << "\n";
  }
  return ss.str();
}