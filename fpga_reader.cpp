#include "fpga_reader.h"

bool FPGAReader::readPlaceFile(const std::string& file_name) {
    std::ifstream ifs(file_name);
    if (!ifs) {
        std::cerr << "cannot open : " << file_name << "\n";
        return false;
    }

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        Point p;
        if (parsePlacement(line, p)) {
            p.die = (p.y < 120) ? 0 : 1;
            inst_map[p.name] = points.size();
            points.push_back(std::move(p));
        }
    }
    return true;
}

bool FPGAReader::readNetFile(const std::string& file_name) {
    std::ifstream ifs(file_name);
    if (!ifs) {
        std::cerr << "cannot open " << file_name << "\n";
        return false;
    }
    std::string content((std::istreambuf_iterator<char>(ifs)),
                        std::istreambuf_iterator<char>());
    parseNet(content);
    return true;
}

bool FPGAReader::parsePlacement(const std::string& line, Point& p) {
    std::istringstream ss(line);
    std::string name, x, z, last;
    int y;

    if (!(ss >> name >> x >> y >> z)) return false;
    p.name = name;
    p.y = y;
    p.is_fixed = false;
    if (ss >> last && last == "FIXED") p.is_fixed = true;

    return true;
}

void FPGAReader::parseNet(const std::string& content) {
    std::string compact = content;
    std::replace(compact.begin(), compact.end(), '\n', ' ');
    std::regex inst_re(R"([\w]+\s+([\w.]+)\s*\((.*?)\);)", std::regex::icase);
    auto it  = std::sregex_iterator(compact.begin(), compact.end(), inst_re);
    auto end = std::sregex_iterator();

    for (; it != end; ++it) {
        const std::smatch& m = *it;
        if (m.size() < 2)
            continue;

        std::string inst_name = m[1].str();
        std::string conn      = m[2].str();
        std::vector<std::string> nets;
        std::regex net_re(R"(\.\w+\(([\w.]+)\))");
        auto net_it  = std::sregex_iterator(conn.begin(), conn.end(), net_re);
        auto net_end = std::sregex_iterator();

        for (; net_it != net_end; ++net_it) {
            nets.push_back((*net_it)[1].str());
        }

        auto map_it = inst_map.find(inst_name);
        if (map_it != inst_map.end()) {
            size_t idx = map_it->second;
            points[idx].nets = std::move(nets);
        }
    }
}

void FPGAReader::buildIndexes() {
    net_map.clear();
    inst_map.clear();

    for (size_t i = 0; i < points.size(); ++i) {
        inst_map[points[i].name] = i;
        for (const auto& net : points[i].nets) {
            net_map[net].push_back(i);
        }
    }
}

void FPGAReader::print_Info() const {
    int die0_count = 0, die1_count = 0, fixed_count = 0;
    for (const auto& p : points) {
        if (p.is_fixed) fixed_count++;
        if (p.die == 0) die0_count++;
        else die1_count++;
    }
    for (int die = 0; die < 2; ++die) {
        std::string fname = "..\\die" + std::to_string(die) + ".txt";
        std::ofstream ofs(fname);
        if (!ofs) {
            return;
        }
        for (const auto& p : points) {
            if (p.die != die) continue;
            ofs << p.name << "\n";
        }
    }
}

void FPGAReader::FM() {
    int movable_total = 0;
    for (const auto& p : points)
        if (!p.is_fixed)
            ++movable_total;
    if (movable_total == 0) {
        return;
    }

    initial_Gain();
    int initial_cost = calculateCutSize();
    int iteration = 0;
    const int MAX_ITERATIONS = 100;
    const int MIN_IMPROVEMENT = 1;
    int best_overall_cost = initial_cost;
    std::vector<int> best_assignment(points.size());

    for (size_t i = 0; i < points.size(); ++i) {
        best_assignment[i] = points[i].die;
    }

    while (iteration < MAX_ITERATIONS) {
        iteration++;

        const size_t n = points.size();
        std::vector<char> locked(n, 0);
        std::vector<int>  orig_die(n);

        struct Move {
            size_t idx;
            int from;
            int to;
            int gain;
        };
        std::vector<Move> moves;
        moves.reserve(movable_total);

        for (size_t i = 0; i < n; ++i)
            orig_die[i] = points[i].die;

        auto rollback = [&](int inclusive) {
            for (auto t = static_cast<size_t>(inclusive + 1); t < moves.size(); ++t) {
                const auto& mv = moves[t];
                points[mv.idx].die = mv.from;
            }
            updateAllGains();
        };

        int moves_made = 0;
        while (moves_made < movable_total) {
            int best_gain = INT_MIN;
            size_t best_idx = SIZE_MAX;
            int best_from = -1, best_to = -1;

            for (size_t i = 0; i < n; ++i) {
                if (locked[i]) continue;
                const auto& p = points[i];
                if (p.is_fixed) continue;

                int g = cached_gains[i];
                if (g > best_gain) {
                    best_gain = g;
                    best_idx = i;
                    best_from = p.die;
                    best_to   = 1 - p.die;
                }
            }

            if (best_idx == SIZE_MAX) {
                break;
            }

            locked[best_idx] = 1;
            points[best_idx].die = best_to;
            moves.push_back({best_idx, best_from, best_to, best_gain});
            updateGain(best_idx);
            moves_made++;
        }

        if (moves.empty()) {
            break;
        }

        int cum = 0, best_cum = INT_MIN, best_k = -1;
        for (size_t k = 0; k < moves.size(); ++k) {
            cum += moves[k].gain;
            if (cum > best_cum) {
                best_cum = cum;
                best_k = static_cast<int>(k);
            }
        }

        if (best_k >= 0 && best_cum > MIN_IMPROVEMENT) {
            rollback(best_k);
            int current_cost = calculateCutSize();

            if (current_cost < best_overall_cost) {
                best_overall_cost = current_cost;
                for (size_t i = 0; i < points.size(); ++i) {
                    best_assignment[i] = points[i].die;
                }
            }
            std::cout << std::endl;
        } else {
            // 回退所有移动
            for (const auto& mv : moves) {
                points[mv.idx].die = mv.from;
            }
            updateAllGains();
            break;
        }
    }

    for (size_t i = 0; i < points.size(); ++i) {
        points[i].die = best_assignment[i];
    }
}

void FPGAReader::initial_Gain() {
    const size_t n = points.size();
    net_cache.clear();
    cached_gains.clear();
    cached_gains.resize(n);

    for (const auto& kv : net_map) {
        net_cache[kv.first] = NetStats(kv.first);
    }
    for (size_t i = 0; i < n; ++i) {
        cached_gains[i] = Point_Gain(i);
    }
    gains_initialized = true;
}

std::array<int, 2> FPGAReader::NetStats(const std::string& net_id) const {
    auto it = net_map.find(net_id);
    if (it == net_map.end()) {
        return {0, 0};
    }

    std::array<int, 2> counts = {0, 0};
    for (size_t point_idx : it->second) {
        int die = points[point_idx].die;
        if (die == 0 || die == 1)
            ++counts[die];
    }
    return counts;
}

int FPGAReader::Point_Gain(size_t point_idx) const {
    const Point& p = points[point_idx];
    if (p.is_fixed)
        return INT_MIN / 4;

    int gain = 0;
    const int current_die = p.die;
    const int target_die = 1 - current_die;

    for (const auto& net_id : p.nets) {
        std::array<int, 2> counts;
        auto cache_it = net_cache.find(net_id);
        if (cache_it != net_cache.end()) {
            counts = cache_it->second;
        } else {
            counts = NetStats(net_id);
        }

        int before_cut = (counts[0] > 0 && counts[1] > 0) ? 1 : 0;
        std::array<int, 2> after_counts = counts;
        after_counts[current_die]--;
        after_counts[target_die]++;

        int after_cut = (after_counts[0] > 0 && after_counts[1] > 0) ? 1 : 0;
        gain += (before_cut - after_cut);
    }
    return gain;
}

void FPGAReader::updateGain(size_t moved_idx) {
    if (!gains_initialized) {
        initial_Gain();
        return;
    }

    const Point& moved_point = points[moved_idx];
    std::unordered_set<size_t> affected_points;
    affected_points.insert(moved_idx);

    for (const auto& net_id : moved_point.nets) {
        net_cache[net_id] = NetStats(net_id);
        auto it = net_map.find(net_id);
        if (it != net_map.end()) {
            for (size_t idx : it->second) {
                affected_points.insert(idx);
            }
        }
    }
    for (size_t idx : affected_points) {
        cached_gains[idx] = Point_Gain(idx);
    }
}

void FPGAReader::updateAllGains() const {
    for (auto& kv : net_cache) {
        kv.second = NetStats(kv.first);
    }
    for (size_t i = 0; i < points.size(); ++i) {
        cached_gains[i] = Point_Gain(i);
    }
}

int FPGAReader::calculateCutSize() const {
    int cut_size = 0;
    for (const auto& kv : net_map) {
        const auto& net_points = kv.second;
        if (net_points.size() < 2)
            continue;

        bool has_die0 = false, has_die1 = false;
        for (size_t idx : net_points) {
            if (points[idx].die == 0) has_die0 = true;
            else if (points[idx].die == 1) has_die1 = true;

            if (has_die0 && has_die1) {
                cut_size++;
                break;
            }
        }
    }
    return cut_size;
}