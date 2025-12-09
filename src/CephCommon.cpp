#include "CephCommon.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex> // for std::scoped_lock
#include <sstream>
#include <string>

namespace fs = std::filesystem;

std::istream &operator>>(std::istream &is, std::vector<unsigned int> &v) {
  unsigned int osd_id;
  char ch;
  std::vector<unsigned int> members;

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
unsigned int PGShard::get_pg_id() const { return pg->get_id(); }

unsigned int PGShard::get_osd_id() const { return osd_id; }

bool PGShard::is_acting() {

  auto acting_members = pg->get_acting().members;
  return std::find(acting_members.begin(), acting_members.end(), this) !=
         acting_members.end();
}

unsigned long long int PGShard::get_objects() const { return objects; }

// PG
PG::PG(std::string line) {
  std::stringstream ss(line);
  // eg.
  // 3.18 raw ([123,73,98], p123) up ([123,73,98], p123) acting ([123,73,98],
  // p123)
  std::string temp;
  std::vector<unsigned int> members;

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
  set_up(members); // data types don't make sense here

  // "acting ([123,73,98], p123)""
  ss >> temp;
  ss >> members;
  set_acting(members);
}

void PG::set_up(const std::vector<unsigned int> &v) {
  xbt_assert(v.size() >= 1);
  const std::scoped_lock lock(*mutex_);

  PGShardSet new_up;
  // if acting[i] == v[i], we reuse the shard
  for (int i = 0; i < v.size(); ++i) {
    if (i < acting.members.size() && acting.members[i]->get_osd_id() == v[i]) {
      new_up.members.push_back(acting.members[i]);
    } else {
      auto shard = std::make_unique<PGShard>(this, v[i], i, 0);
      new_up.members.push_back(shard.get());
      shards.insert(std::move(shard));
    }
  }
  up = new_up;
  prune_shards();
}

void PG::set_acting(const std::vector<unsigned int> &v) {
  xbt_assert(v.size() >= 1);
  const std::scoped_lock lock(*mutex_);

  PGShardSet new_acting;
  // if up[i] == v[i], we reuse the shard
  for (int i = 0; i < v.size(); ++i) {
    if (i < up.members.size() && up.members[i]->get_osd_id() == v[i]) {
      new_acting.members.push_back(up.members[i]);
    } else {
      auto shard = std::make_unique<PGShard>(this, v[i], i, 0);
      new_acting.members.push_back(shard.get());
      shards.insert(std::move(shard));
    }
  }
  acting = new_acting;
  prune_shards();
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
  const std::scoped_lock lock(*mutex_);
  std::ostringstream ss;
  ss << "PG " << id << "(" << std::hex << id << ") ";
  ss << "up " << up.to_string() << " ";
  ss << "acting " << acting.to_string();
  return ss.str();
}

// PGMap
PGMap::PGMap(std::string path, unsigned int pool) {
  // parse path + assert
  xbt_assert(fs::is_regular_file(path));

  // open
  std::ifstream file(path);
  std::string line;
  unsigned int current_pool;
  unsigned int pg_num;
  unsigned int pgs_to_parse;
  std::string discard;

  // parse
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string first_word;
    ss >> first_word;
    if (first_word == "pool") {
      ss >> current_pool;
      if (current_pool == pool) {
        ss >> discard >> pg_num;
        pgs_to_parse = pg_num;
      }

    } else if (current_pool == pool && pgs_to_parse > 0) {
      pgs.push_back(std::make_unique<PG>(line));
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

void PGMap::rebuild_primary_osd_to_pg_index() {
  const std::scoped_lock lock(*mutex_);
  primary_osd_to_pg_index.clear();
  for (const auto &pg_ptr : pgs) {
    unsigned int primary_osd =
        pg_ptr->get_acting_ids().at(0); // first is primary
    primary_osd_to_pg_index[primary_osd].insert(pg_ptr.get());
  }
}

size_t PGMap::size() const {
  const std::scoped_lock lock(*mutex_);
  return pgs.size();
}

bool PG::needs_backfill() const {
  const std::scoped_lock lock(*mutex_);
  return up.members != acting.members;
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
  const std::scoped_lock lock(*mutex_);
  std::ostringstream ss;
  ss << "PGMap (Size " << pgs.size() << "):\n";
  for (const auto &pg : pgs) {
    ss << "  " << pg->to_string() << "\n";
  }
  return ss.str();
}