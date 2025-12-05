#pragma once
#include <simgrid/s4u.hpp>
#include <variant>

struct OSDSet {
  std::vector<unsigned int> members;
  unsigned int primary() const { return members.at(0); };
  std::string to_string() const {
    std::ostringstream ss;
    ss << "[";
    for (size_t i = 0; i < members.size(); ++i) {
      ss << members[i] << (i < members.size() - 1 ? "," : "");
    }
    ss << "]";
    return ss.str();
  }
};

// PG struct
class PG {
  unsigned int id;
  OSDSet up;
  OSDSet acting;

public:
  explicit PG(std::string line);
  void set_up(const OSDSet &up_set) { up = up_set; }
  void set_acting(const OSDSet &acting_set) { acting = acting_set; }
  const OSDSet &get_up() const { return up; }
  const OSDSet &get_acting() const { return acting; }
  unsigned int get_id() const { return id; }
  bool needs_backfill() const;
  std::string to_string() const;
};

// pg map
class PGMap {
  std::vector<PG> pgs;

public:
  PGMap(std::string path, unsigned int pool);

  const std::vector<OSDSet> get_up() const;
  const std::vector<OSDSet> get_acting() const;
  void set_up(const std::vector<OSDSet> &up_list);
  void set_acting(const std::vector<OSDSet> &acting_list);

  size_t size() const;
  int find_pg_up(int pg) const;
  int find_pg_acting(int pg) const;
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
