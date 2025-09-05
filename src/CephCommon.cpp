#include "CephCommon.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

namespace fs = std::filesystem;


std::istream& operator>>(std::istream& is, OSDSet& set) {
    OSDSet temp;
    char ch;

    // "(["
    if(!(is >> ch && ch == '(' && is >> ch && ch == '[')) {
        is.setstate(std::ios_base::failbit);
        return is;
    }

    // "1,20,30"
    unsigned int member;
    is >> member;
    temp.members.push_back(member);
    while (is.peek() == ',') {
        is >> ch;
        is >> member;
        temp.members.push_back(member);
    }

    // "], p"
    if (!(is >> ch && ch == ']' && is >> ch && ch == ',' && is >> ch && ch == 'p')) {
        is.setstate(std::ios_base::failbit);
        return is;
    }

    // "1"
    is >> temp.primary;

    // "], p"
    if (!(is >> ch && ch == ')')) {
        is.setstate(std::ios_base::failbit);
        return is;
    }

    // only on success do we touch the set reference
    if (is)
        set = std::move(temp);
    return is;
}

PG::PG(std::string line) {
    std::stringstream ss(line);
    // eg.
    // 3.18 raw ([123,73,98], p123) up ([123,73,98], p123) acting ([123,73,98], p123)
    std::string temp;

    // "3.18"
    char del = '.';
    getline(ss, temp, '.');
    ss >> temp;
    id = std::stoul(temp, nullptr, 16);

    // "raw"
    ss >> temp;

    // "([123,73,98], p123)""
    ss >> acting;

    // "up ([123,73,98], p123)"
    ss >> temp;
    ss >> up;
    
    // "acting ([123,73,98], p123)""
    ss >> temp;
    ss >> acting;
}

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
    while(std::getline(file, line)) {
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
            PG pg(line);
            pgs.push_back(pg);
            pgs_to_parse--;
        }
    }
}

size_t PGMap::size() const {
    return pgs.size();
}

int PGMap::find_pg_up(int pg) const {
    return -1;
}

int PGMap::find_pg_acting(int pg) const {
    return -1;
}