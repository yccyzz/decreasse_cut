#include "fpga_reader.h"

bool FPGAReader::readPlaceFile(const std::string& file_name) {
    std::ifstream ifs(file_name);
    if (!ifs) {
        std::cerr << "cannot open : " << file_name << "\n";
        return false;
    }

    // 清空所有数据
    clear();

    std::string line;
    points.reserve(50000);
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        Point p;
        if (parsePlacement(line, p)) {
            if (p.y < 120) {
                p.die = 0;
            }else if(p.y < 240) {
                p.die = 1;
            }else if(p.y < 360) {
                p.die = 2;
            }else {
                p.die = 3;
            }
            points.insert({p.name, std::move(p)});
        }
    }
    detectDieMode();
    return true;
}

void FPGAReader::detectDieMode() {
    std::set<int> dies_used;
    for (const auto& [name, point] : points) {
        dies_used.insert(point.die);
    }

    // 如果只使用了die 0和1，则为2-die模式；否则为4-die模式
    if (dies_used.size() == 2 && dies_used.count(0) && dies_used.count(1)) {
        num_dies = 2;
        std::cout << "Detected 2-die mode" << std::endl;
    } else {
        num_dies = 4;
        std::cout << "Detected 4-die mode "<<std::endl;
    }
}

void FPGAReader::clear() {
    points.clear();
    net_map.clear();
    net_cache.clear();
    cached_gains.clear();
    sll_cache.clear();
    num_dies = 2; // 默认2-die模式
    max_gain = -INT32_MAX;
    iteration_counter = 0;
    cached_cut_size = -1;
    network_stats_cache_valid = false;
}

bool FPGAReader::readNetFile(std::string& file_name) {
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

void FPGAReader::parseNet(std::string& content) {
    std::replace(content.begin(), content.end(), '\n', ' ');
    std::regex inst_re(R"([\w]+\s+([\w.]+)\s*\((.*?)\);)", std::regex::icase);
    std::regex net_re(R"(\.\w+\((net_\d+)\))");

    // 使用string来临时存储，然后转为string
    std::unordered_map<std::string, std::vector<std::string>> temp_point_nets;
    std::unordered_map<std::string, std::unordered_set<std::string>> temp_net_to_points;
    temp_point_nets.reserve(points.size());
    temp_net_to_points.reserve(points.size());
    net_map.reserve(50000);
    // 第一步：解析每个实例的网络连接
    for (std::sregex_iterator it(content.begin(), content.end(), inst_re), end; it != end; ++it) {
        const auto& m = *it;
        std::string inst_name = m[1].str();
        std::string conn = m[2].str();

        if (!points.contains(inst_name)) {
            continue;
        }

        for (auto net_it = std::sregex_iterator(conn.begin(), conn.end(), net_re);
             net_it != std::sregex_iterator(); ++net_it) {
            std::string net_name = (*net_it)[1].str();
            temp_point_nets[inst_name].emplace_back(net_name);
            temp_net_to_points[net_name].insert(inst_name);
        }
    }

    // 第二步：构建最终的net_map（从网络到点的映射）
    for (const auto& [net_name, point_list] : temp_net_to_points) {
        auto& net_points = net_map[net_name];
        net_points.reserve(point_list.size());
        for (const auto& point_name_str : point_list) {
            if (points.contains(point_name_str)) {
                net_points.insert(point_name_str);
            }
        }
    }

    // 第三步：为每个Point填充nets字段
    for (auto& [point_name_sv, point] : points) {
        point.nets.clear();
        if (temp_point_nets.contains(point.name)) {
            const auto& net_list = temp_point_nets[point.name];
            point.nets.reserve(net_list.size());
            for (const auto& net_name : net_list) {
                point.nets.emplace_back(net_name);
            }
        }
    }

    // 解析完成后，使网络统计缓存失效
    invalidateNetworkStatsCache();
}

// ==================== 优化的网络统计计算 ====================

void FPGAReader::updateNetworkStatsCache() const {
    if (network_stats_cache_valid) return;

    sll_cache.clear();
    cached_cut_size = 0;

    // 初始化所有die对的SLL为0
    for (int i = 0; i < num_dies; ++i) {
        for (int j = i + 1; j < num_dies; ++j) {
            sll_cache[{i, j}] = 0;
        }
    }

    // 一次遍历计算cut size和所有SLL
    for (const auto& [net_name, point_set] : net_map) {
        std::set<int> dies_in_net;

        // 收集该网络涉及的所有die
        for (const std::string& point_name : point_set) {
            if (auto it = points.find(point_name); it != points.end()) {
                dies_in_net.insert(it->second.die);
            }
        }

        // 如果网络跨越多个die
        if (dies_in_net.size() > 1) {
            cached_cut_size++;  // cut size +1

            // 更新所有相关die对的SLL
            for (auto it1 = dies_in_net.begin(); it1 != dies_in_net.end(); ++it1) {
                for (auto it2 = std::next(it1); it2 != dies_in_net.end(); ++it2) {
                    int die1 = *it1, die2 = *it2;
                    if (die1 > die2) std::swap(die1, die2);
                    sll_cache[{die1, die2}]++;
                }
            }
        }
    }

    network_stats_cache_valid = true;
}

void FPGAReader::invalidateNetworkStatsCache() {
    network_stats_cache_valid = false;
    net_cache.clear(); // 同时清空网络缓存
}

int FPGAReader::calculateCutSize() const {
    updateNetworkStatsCache();
    return cached_cut_size;
}

int FPGAReader::calculateSLLBetweenDies(int die1, int die2) const {
    if (die1 == die2) return 0;

    // 确保die1 < die2，标准化键值
    if (die1 > die2) std::swap(die1, die2);

    updateNetworkStatsCache();

    auto it = sll_cache.find({die1, die2});
    return (it != sll_cache.end()) ? it->second : 0;
}

// ==================== SLL约束检查 ====================

int FPGAReader::calculateSLLChangeForDiePair(const std::string& point_name, int from_die, int to_die, int die1, int die2) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return 0;

    const Point& point = point_it->second;
    int sll_change = 0;

    // 遍历该点连接的所有网络
    for (const std::string& net_name : point.nets) {
        auto stats = NetStats(net_name);

        // 计算移动前该网络是否连接die1和die2
        bool before_connects_die1 = (stats[die1] > 0);
        bool before_connects_die2 = (stats[die2] > 0);
        bool before_sll = before_connects_die1 && before_connects_die2;

        // 模拟移动后的状态
        stats[from_die]--;
        stats[to_die]++;

        // 计算移动后该网络是否连接die1和die2
        bool after_connects_die1 = (stats[die1] > 0);
        bool after_connects_die2 = (stats[die2] > 0);
        bool after_sll = after_connects_die1 && after_connects_die2;

        // 计算SLL变化
        if (!before_sll && after_sll) {
            sll_change += 1;  // 新增一个SLL
        } else if (before_sll && !after_sll) {
            sll_change -= 1;  // 减少一个SLL
        }
        // 其他情况SLL不变
    }

    return sll_change;
}

bool FPGAReader::isValidMoveForSLL(const std::string& point_name, int from_die, int to_die) const {
    if (from_die == to_die) return true;

    // 对于2-die模式，只需要检查die0和die1之间
    if (num_dies == 2) {
        int current_sll = calculateSLLBetweenDies(0, 1);
        int sll_change = calculateSLLChangeForDiePair(point_name, from_die, to_die, 0, 1);
        return (current_sll + sll_change <= MAX_SLL_BETWEEN_DIES);
    }

    // 对于4-die模式，检查所有die对
    for (int i = 0; i < num_dies; ++i) {
        for (int j = i + 1; j < num_dies; ++j) {
            int current_sll = calculateSLLBetweenDies(i, j);
            int sll_change = calculateSLLChangeForDiePair(point_name, from_die, to_die, i, j);

            if (current_sll + sll_change > MAX_SLL_BETWEEN_DIES) {
                return false;  // 违反约束
            }
        }
    }

    return true;
}

// ==================== FM算法实现 ====================

void FPGAReader::FM() {
    if (points.empty()) return;
    auto cutsize_start = std::chrono::high_resolution_clock::now();
    int initial_cut_size = calculateCutSize();
    auto cutsize_end = std::chrono::high_resolution_clock::now();
    std::cout << "Initial cut size: " << initial_cut_size << std::endl;

    // 打印初始SLL统计
//    printSLLStats();

    auto init_gain_start = std::chrono::high_resolution_clock::now();
    // 初始化排序的增益表
    initSortedGains();
    auto init_gain_end = std::chrono::high_resolution_clock::now();

    // FM Pass - 严格只移动正增益点，并检查SLL约束
    while (!gain_map.empty() ) {

        if (gain_queue.empty()) {
            updateGainQueue();
            if (gain_queue.empty()) {
                break; // 没有更多可移动的点
            }
        }

        auto top_item = gain_queue.top();
        gain_queue.pop();

        std::string current_max_point_name = top_item.second;
        int current_max_gain = Point_Gain(current_max_point_name);

        // 检查该点是否已被固定（可能在队列更新间隔内被处理过）
        if (points[current_max_point_name].is_fixed) {
            continue;
        }

        // 验证该点在gain_map中仍然存在（确保数据一致性）
        if (gain_map[current_max_gain].find(current_max_point_name) == gain_map[current_max_gain].end()) {
            continue; // 该点可能已经被更新，跳过
        }

        // 确定目标die
        int target_die;
        if (num_dies == 2) {
            target_die = 1 - points[current_max_point_name].die;
        } else {
            // 4-die模式：获取最佳目标die
            auto [gain, best_target] = Point_Gain_Multi(current_max_point_name);
            target_die = best_target;
            // 重新验证gain值
            current_max_gain = gain;
        }

        // **新增：SLL约束检查**
        if (!isValidMoveForSLL(current_max_point_name,
                               points[current_max_point_name].die,
                               target_die)) {
            // 违反SLL约束，标记为固定，不允许移动
            points[current_max_point_name].is_fixed = true;
            // 从gain_map中移除该点
            gain_map[current_max_gain].erase(current_max_point_name);
            continue;
        }

        // 检查是否为正增益，如果不是则退出
        if (current_max_gain < 0) {
            break;
        }

        // 从gain_map中移除该点
        gain_map[current_max_gain].erase(current_max_point_name);

        // 执行移动
        points[current_max_point_name].die = target_die;
        points[current_max_point_name].is_fixed = true;

        // 更新受影响点的增益（这会更新gain_map）
        updateSortedGains(current_max_point_name);

        // 定期更新队列以保持gain_queue和gain_map同步
        iteration_counter++;
        if (iteration_counter >= UPDATE_FREQUENCY) {
            updateGainQueue();
            iteration_counter = 0;
        }
    }

    std::cout << "Final cut size: " << calculateCutSize() << std::endl;
    // 打印最终SLL统计
//    printSLLStats();
}

void FPGAReader::initSortedGains() {
    gain_map.clear();
    net_cache.clear();
    cached_gains.clear();
    cached_gains.reserve(points.size());

    // 计算所有点的增益并排序
    for (const auto& [name, point] : points) {
        int gain = Point_Gain(name);
        if (gain > max_gain) {
            max_gain = gain;
            max_gain_point_name = name;
        }
        cached_gains[name] = gain;
        gain_map[gain].insert(name);
    }
}

void FPGAReader::updateSortedGains(const std::string& moved_name) {
    const auto& moved_point = points[moved_name];

    std::unordered_set<std::string > affected_points;
    affected_points.reserve(500);
    for (const std::string& net_name : moved_point.nets) {
        net_cache.erase(net_name);

        if (auto net_it = net_map.find(net_name); net_it != net_map.end()) {
            for (const std::string& point_name : net_it->second) {
                if (point_name != moved_name && !points[point_name].is_fixed) {
                    affected_points.insert(point_name);
                }
            }
        }
    }

    for (const std::string& point_name : affected_points) {
        // 移除旧增益
        int old_gain = cached_gains[point_name];
        auto& range = gain_map[old_gain];
        if (range.contains(point_name)) {
            range.erase(point_name);
        }

        auto [new_gain, target_die] = Point_Gain_Multi(point_name);
        cached_gains[point_name] = new_gain;
        auto& new_gain_range = gain_map[new_gain];
        new_gain_range.insert(point_name);
    }
    cached_gains.erase(moved_name);

    // 移动后使网络统计缓存失效
    invalidateNetworkStatsCache();
}

std::vector<int> FPGAReader::NetStats(const std::string& net_id) const {
    if (auto it = net_cache.find(net_id); it != net_cache.end()) {
        return it->second;
    }

    std::vector<int> stats(num_dies, 0);
    if (auto net_it = net_map.find(net_id); net_it != net_map.end()) {
        for (const std::string& point_name : net_it->second) {
            if (auto point_it = points.find(point_name); point_it != points.end()) {
                int die = point_it->second.die;
                if (die >= 0 && die < num_dies) {
                    stats[die]++;
                }
            }
        }
    }

    net_cache[net_id] = stats;
    return stats;
}

int FPGAReader::Point_Gain(const std::string& point_name) const {
    // 2-die模式的原始实现，保持向后兼容
    auto [gain, target_die] = Point_Gain_Multi(point_name);
    return gain;
}

std::pair<int, int> FPGAReader::Point_Gain_Multi(const std::string& point_name) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return {0, -1};

    const Point& point = point_it->second;
    int current_die = point.die;
    int best_gain = -INT_MAX;
    int best_target_die = current_die;

    // 如果是2-die模式，只考虑另一个die
    if (num_dies == 2) {
        int target_die = 1 - current_die;
        int gain = calculateGainForMove(point_name, current_die, target_die);
        return {gain, target_die};
    }

    // 4-die模式：尝试移动到所有其他die，选择gain最大的
    for (int target_die = 0; target_die < num_dies; ++target_die) {
        if (target_die == current_die) continue;

        int gain = calculateGainForMove(point_name, current_die, target_die);
        if (gain > best_gain) {
            best_gain = gain;
            best_target_die = target_die;
        }
    }

    return {best_gain, best_target_die};
}

int FPGAReader::calculateGainForMove(const std::string& point_name, int from_die, int to_die) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return 0;

    const Point& point = point_it->second;
    int gain = 0;

    for (const std::string& net_name : point.nets) {
        auto stats = NetStats(net_name);

        // 计算移动前该网络是否跨die
        int dies_before = 0;
        for (int i = 0; i < num_dies; ++i) {
            if (stats[i] > 0) dies_before++;
        }
        bool before_cut = dies_before > 1;

        // 模拟移动后的情况
        stats[from_die]--;
        stats[to_die]++;

        // 计算移动后该网络是否跨die
        int dies_after = 0;
        for (int i = 0; i < num_dies; ++i) {
            if (stats[i] > 0) dies_after++;
        }
        bool after_cut = dies_after > 1;

        if (before_cut && !after_cut) {
            gain += 1;  // 消除割边
        } else if (!before_cut && after_cut) {
            gain -= 1;  // 增加割边
        }
    }

    return gain;
}

void FPGAReader::updateGainQueue() {
    // 清空现有队列
    std::priority_queue<std::pair<int, std::string>> empty_queue;
    std::swap(gain_queue, empty_queue);

    // 用当前增益值填充队列，按gain值从高到低排序
    for (const auto& [gain_value, point_set] : gain_map) {
        if (point_set.empty()) {
            continue;
        }
        for (const auto& point_name : point_set) {
            // 确保只添加未固定的点
            if (!points[point_name].is_fixed) {
                gain_queue.push(std::make_pair(gain_value, point_name));
            }
        }
    }
}

// ==================== 输出和调试函数 ====================

void FPGAReader::print_Info() const {
    int output_dies = (num_dies == 2) ? 2 : 4;

    for (int die = 0; die < output_dies; ++die) {
        std::string fname = "..\\die" + std::to_string(die) + ".txt";
        std::ofstream ofs(fname);
        if (!ofs) {
            return;
        }
        for (const auto& p : points) {
            if (p.second.die != die) continue;
            ofs << p.second.name << "\n";
        }
    }
}

void FPGAReader::printSLLStats() const {
    std::cout << "\n=== SLL Statistics ===" << std::endl;

    if (num_dies == 2) {
        int sll = calculateSLLBetweenDies(0, 1);
        std::cout << "Die 0 <-> Die 1: " << sll << " SLLs";
        if (sll > MAX_SLL_BETWEEN_DIES) {
            std::cout << " (VIOLATION! Limit: " << MAX_SLL_BETWEEN_DIES << ")";
        } else {
            std::cout << " (OK, Limit: " << MAX_SLL_BETWEEN_DIES << ")";
        }
        std::cout << std::endl;
    } else {
        // 4-die模式
        for (int i = 0; i < num_dies; ++i) {
            for (int j = i + 1; j < num_dies; ++j) {
                int sll = calculateSLLBetweenDies(i, j);
                std::cout << "Die " << i << " <-> Die " << j << ": " << sll << " SLLs";
                if (sll > MAX_SLL_BETWEEN_DIES) {
                    std::cout << " (VIOLATION! Limit: " << MAX_SLL_BETWEEN_DIES << ")";
                } else {
                    std::cout << " (OK, Limit: " << MAX_SLL_BETWEEN_DIES << ")";
                }
                std::cout << std::endl;
            }
        }
    }
    std::cout << "=====================\n" << std::endl;
}