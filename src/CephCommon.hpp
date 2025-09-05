#pragma once
#include <simgrid/s4u.hpp>
#include <variant>

struct OSDSet {
    unsigned int primary;
    std::vector<unsigned int> members;
};

// PG struct
class PG {
    unsigned int id;
    OSDSet acting;
    OSDSet up;

public:
    explicit PG(std::string line);
};

// pg map
class PGMap {
    std::vector<PG> pgs;

public:
    PGMap(std::string path, unsigned int pool);
    // methods for modifying pgmap which return a new instance of pgmap
        // eg. one PG's acting is changed helper. 
    
    size_t size() const;
    int find_pg_up(int pg) const;
    int find_pg_acting(int pg) const;
};

// message types
struct KillMsg {};
struct SubscribeToPGMapChangeMsg {};
struct SubscribeToPGChangeMsg {};

struct PGMapNotification {
    PGMap* pg_map;
};

struct PGNotification {
    int pg_id;
    int acting;
};


using MessagePayload = std::variant<
    KillMsg,
    SubscribeToPGMapChangeMsg,
    SubscribeToPGChangeMsg,
    PGMapNotification,
    PGNotification
>;

struct Message {
    std::string sender;
    MessagePayload payload;
};

template <typename T, typename... Args>
Message* make_message(Args&&... args)
{
    std::string my_name = simgrid::s4u::this_actor::get_host()->get_name();
    return new Message{
        .sender = my_name,
        .payload = T{std::forward<Args>(args)...}
    };
}
