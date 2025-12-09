#pragma once
#include <simgrid/s4u.hpp>
#include <variant>
#include <vector>

namespace sg4 = simgrid::s4u;

// PG forward declaration
class PG;

class PGShard {
  PG *pg;
  unsigned int osd_id;
  unsigned int index;
  unsigned long long int
      objects; // non-acting (up) pg shard might be incomplete

public:
  PGShard(PG *pg, unsigned int osd_id, unsigned int index, unsigned int objects)
      : pg(pg), osd_id(osd_id), index(index), objects(objects) {}
  unsigned int get_pg_id() const;
  unsigned int get_osd_id() const;
  unsigned int get_index() const { return index; }
  bool is_acting();
  unsigned long long int get_objects() const;
};

// PG Shard set
struct PGShardSet {
  std::vector<PGShard *> members;
  PGShard *primary() const { return members.at(0); };

  bool contains_shard(PGShard *shard) const {
    return members.size() > shard->get_index() &&
           members.at(shard->get_index()) == shard;
  }

  std::string to_string() const {
    std::ostringstream ss;
    ss << "[";
    for (size_t i = 0; i < members.size(); ++i) {
      ss << members[i]->get_osd_id() << (i < members.size() - 1 ? "," : "");
    }
    ss << "]";
    return ss.str();
  }
};

// PG class
class PG {
  // mutex is non-recursive, actor must ensure no deadlock
  sg4::MutexPtr mutex_ = sg4::Mutex::create();
  unsigned int id;
  // ownership
  std::set<std::unique_ptr<PGShard>> shards;
  // current state of those shards
  PGShardSet up;
  PGShardSet acting;
  void prune_shards();

public:
  explicit PG(std::string line);
  void set_up(const std::vector<unsigned int> &set);
  void set_acting(const std::vector<unsigned int> &set);
  const PGShardSet get_up() const { return up; }
  const std::vector<unsigned int> get_up_ids() const {
    std::vector<unsigned int> ids;
    for (const auto &shard : up.members) {
      ids.push_back(shard->get_osd_id());
    }
    return ids;
  }
  const std::vector<unsigned int> get_acting_ids() const {
    std::vector<unsigned int> ids;
    for (const auto &shard : acting.members) {
      ids.push_back(shard->get_osd_id());
    }
    return ids;
  }
  const PGShardSet get_acting() const { return acting; }
  unsigned int get_id() const { return id; }
  bool needs_backfill() const;
  std::string to_string() const;
};

// PGMap
class PGMap {
  // mutex is non-recursive, actor must ensure no deadlock
  sg4::MutexPtr mutex_ = sg4::Mutex::create();
  std::vector<std::unique_ptr<PG>> pgs;
  std::map<unsigned int, std::set<PG *>> primary_osd_to_pg_index;

public:
  PGMap(std::string path, unsigned int pool);

  size_t size() const;
  const std::vector<PGShardSet> get_up() const;
  const std::vector<PGShardSet> get_acting() const;

  const PG *get_pg(int pg) const { return pgs.at(pg).get(); };
  PG *get_pg(int pg) { return pgs.at(pg).get(); };

  void rebuild_primary_osd_to_pg_index();
  void refresh_primary_osd_to_pg_index_for_pg(unsigned int pg_id);

  std::string primary_osds_to_pgs_string() const;
  std::string to_string() const;
};

// message types
// c++17 variant visitor helper
template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

struct KillMsg {};
struct SubscribeToPGMapChangeMsg {};
struct SubscribeToPGChangeMsg {};

struct PGMapNotification {
  PGMap *pg_map;
};

struct PGNotification {
  int pg_id;
  int acting;
};

using MessagePayload =
    std::variant<KillMsg, SubscribeToPGMapChangeMsg, SubscribeToPGChangeMsg,
                 PGMapNotification, PGNotification>;

struct Message {
  std::string sender;
  MessagePayload payload;
};

template <typename T, typename... Args> Message *make_message(Args &&...args) {
  return new Message{.sender = simgrid::s4u::this_actor::get_host()->get_name(),
                     .payload = T{std::forward<Args>(args)...}};
}
