#pragma once
#include <simgrid/s4u.hpp>
#include <sstream>
#include <variant>
#include <vector>

namespace sg4 = simgrid::s4u;

// helper functions
int sender_str_to_int(std::string sender);

// PG forward declaration
class PG;

class PGShard {
  PG *pg;
  int osd_id;
  int index;

public:
  PGShard(PG *pg, int osd_id, int index)
      : pg(pg), osd_id(osd_id), index(index) {}
  int get_pg_id() const;
  int get_osd_id() const;
  int get_index() const { return index; }
  bool is_acting();
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
  int id;
  // ownership
  std::set<std::unique_ptr<PGShard>> shards;
  // current state of those shards
  PGShardSet up;
  PGShardSet acting;

  // Backfill state
  int pg_objects;
  int object_size;
  int objects_recovered = 0;
  int objects_recoveries_scheduled = 0;

  void prune_shards();

public:
  explicit PG(std::string line, size_t object_size, size_t pg_objects);
  void init_up(const std::vector<int> &set);
  void update_up(const std::vector<int> &set);
  void init_acting(const std::vector<int> &set);
  void update_acting(const std::vector<int> &set);

  // State transitions
  void on_object_recovered();
  bool schedule_recovery();

  int get_max_osd_id() const {
    int max = 0;
    for (const auto &shard : shards) {
      max = std::max(max, shard->get_osd_id());
    }
    return max;
  }
  size_t get_object_size() const { return object_size; }
  size_t get_objects_per_pg() const { return pg_objects; }
  const PGShardSet get_up() const { return up; }
  const std::vector<int> get_up_ids() const {
    std::vector<int> ids;
    for (const auto &shard : up.members) {
      ids.push_back(shard->get_osd_id());
    }
    return ids;
  }
  const std::vector<int> get_acting_ids() const {
    std::vector<int> ids;
    for (const auto &shard : acting.members) {
      ids.push_back(shard->get_osd_id());
    }
    return ids;
  }
  const PGShardSet get_acting() const { return acting; }
  int get_primary() const { return acting.primary()->get_osd_id(); }
  int get_id() const { return id; }
  bool needs_backfill() const;
  std::string to_string() const;
};

// PGMap
class PGMap {
  // mutex is non-recursive, actor must ensure no deadlock
  sg4::MutexPtr mutex_ = sg4::Mutex::create();
  std::vector<std::unique_ptr<PG>> pgs;
  std::map<int, std::set<PG *>> primary_osd_to_pg_index;
  int max_osd_id;
  size_t object_size;
  size_t pg_objects;

public:
  PGMap(int pool_id, std::string path, size_t object_size, size_t pg_objects);

  // getters
  size_t size() const;
  const std::vector<PGShardSet> get_up() const;
  const std::vector<PGShardSet> get_acting() const;
  const PG *get_pg(int pg) const { return pgs.at(pg).get(); };
  PG *get_pg(int pg) { return pgs.at(pg).get(); };

  std::vector<int> get_osds() const;

  size_t get_object_size() const;
  size_t get_objects_per_pg() const;

  // primary osd to pgs
  const std::set<PG *> primary_osd_get_pgs(int osd_id) const {
    const std::scoped_lock lock(*mutex_);
    std::set<PG *> result;
    auto pg_set = primary_osd_to_pg_index.at(osd_id);
    for (auto pg : pg_set) {
      result.insert(pg);
    }
    return result;
  }
  bool needs_backfill() const;

  // mailboxes
  std::vector<sg4::Mailbox *> osd_mailboxes;
  sg4::Mailbox *get_osd_mailbox(int osd_id) const {
    return osd_mailboxes.at(osd_id);
  }

  // refreshers
  void _update_primary_osd_to_pg_index();
  void update_primary_osd_to_pg_index();
  void refresh_primary_osd_to_pg_index_for_pg(int pg_id);

  std::string primary_osds_to_pgs_string() const;
  std::string to_string() const;
};

// OSD ops and op contexts
enum class OpType {
  BACKFILL,
  REPLICA_WRITE, // works for client write and replication
  CLIENT_READ,
  CLIENT_WRITE,
};

enum class OpState {
  OP_CREATED,
  OP_QUEUED, // not started yet
  OP_WAITING_DISK,
  OP_WAITING_PEER,
  OP_COMPLETED,
};

struct OpContext {
  int local_id;     // unique within this actor
  int client_op_id; // id from the sender (e.g. client or peer osd)
  OpType type;
  int pgid;
  int sender; // positive is osd id, negative is client id
  size_t size;
  double start_time = std::numeric_limits<double>::max();

  OpState state;
  std::set<int> pending_peers;
};

struct Op {
  OpType type;
  int id; // unique only within sender
  int recipient;
  int pgid;
  size_t size;
};

// message types
// c++17 variant visitor helper
template <class... Ts> struct overloaded : Ts... {
  using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

struct KillMsg {};
struct SubscribeToPGMapChangeMsg {};

struct PGMapNotification {
  PGMap *pgmap;
};

struct PGNotification {
  int pg_id;
  int acting;
};

struct OsdOpMsg {
  Op *op;
};

struct OsdOpAckMsg {
  int op_id;
};

using MessagePayload =
    std::variant<KillMsg, SubscribeToPGMapChangeMsg, PGMapNotification,
                 PGNotification, OsdOpMsg, OsdOpAckMsg>;

struct Message {
  std::string sender;
  MessagePayload payload;
};

template <typename T, typename... Args> Message *make_message(Args &&...args) {
  return new Message{.sender = simgrid::s4u::Actor::self()->get_name(),
                     .payload = T{std::forward<Args>(args)...}};
}
