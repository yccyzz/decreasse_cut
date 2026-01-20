#include "fpga_reader.h"

// ==================== 文件读取与初始化 ====================

bool FPGAReader::readPlaceFile(const std::string& file_name) {
    std::ifstream ifs(file_name);
    if (!ifs) {
        std::cerr << "cannot open : " << file_name << "\n";
        return false;
    }
    clear();
    std::string line;
    points.reserve(50000);
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        Point p;
        if(parsePlacement(line, p)) {
            if (!is_coarsing) {
                if (p.y < 120) {
                    p.die = 0;
                }else if(p.y < 240) {
                    p.die = 1;
                }else if(p.y < 360) {
                    p.die = 2;
                }else {
                    p.die = 3;
                }
            }
            points.insert({p.name, std::move(p)});
        }
    }
    detectDieMode();
    return true;
}

bool FPGAReader::parsePlacement(const std::string& line, Point& p) {
    std::istringstream ss(line);
    std::string name, last;
    double x_temp; // 占位
    int y, die;

    if(!is_coarsing){
        // 原始格式: name x y z [FIXED]
        // 这里假设第二列是x，但我们只用y
        std::string x_str, z_str;
        if (!(ss >> name >> x_str >> y >> z_str)) return false;
        p.name = name;
        p.y = y;
        p.is_fixed = false;
        p.is_originally_fixed = false;

        if (ss >> last && last == "FIXED") {
            p.is_fixed = true;
            p.is_originally_fixed = true;
        }
    } else{
        // 粗化格式: name y die [FIXED]
        // 注意：coarsing输出的可能是浮点y，这里最好用double读再转int，或者直接根据die算
        double y_double;
        if (!(ss >> name >> y_double >> die)) return false;
        p.name = name;
        p.die = die;
        p.y = (int)y_double;

        if (die < 0 || die > 3) {
            // std::cout << "Warning: Invalid die " << die << " for " << name << std::endl;
        }

        p.is_fixed = false;
        p.is_originally_fixed = false;

        if (ss >> last && last == "FIXED") {
            p.is_fixed = true;
            p.is_originally_fixed = true;
        }
    }

    return true;
}

void FPGAReader::detectDieMode() {
    std::set<int> dies_used;
    for (const auto& [name, point] : points) {
        dies_used.insert(point.die);
    }

    if (dies_used.size() <= 2 && dies_used.count(0) && dies_used.count(1) && !dies_used.count(2)) {
        num_dies = 2;
    } else {
        num_dies = 4;
    }
}

void FPGAReader::clear() {
    points.clear();
    net_map.clear();
    net_distribution.clear(); // 清除分布缓存
    cached_gains.clear();
    gain_map.clear();

    // 清空优先队列
    std::priority_queue<std::pair<double, std::string>> empty;
    std::swap(gain_queue, empty);

    die_type_count.clear();
    die_type_capacity.clear();
    tracked_resource_types.clear();

    num_dies = 2;
    max_gain = -1e9;
    iteration_counter = 0;
    cached_cut_size = -1;
    network_stats_cache_valid = false;

    base_penalty_coefficient = 1.0;

    initial_pairwise_cuts.clear();
    final_pairwise_cuts.clear();
    execution_time = 0.0;
}

// 【修复栈溢出】使用逐行读取替代正则全匹配
bool FPGAReader::readNetFile(std::string& file_name) {
    std::ifstream ifs(file_name);
    if (!ifs) {
        std::cerr << "cannot open " << file_name << "\n";
        return false;
    }

    std::cout << "Parsing netlist line by line..." << std::endl;

    std::unordered_map<std::string, std::vector<std::string>> temp_point_nets;
    std::unordered_map<std::string, std::unordered_set<std::string>> temp_net_to_points;

    temp_point_nets.reserve(points.size());
    net_map.reserve(points.size() * 2);

    std::string line;
    std::string current_inst_name;
    std::string current_inst_type;

    // 匹配实例头: "LUT6 inst_name ("
    std::regex header_re(R"(^([\w]+)\s+([\w\.]+)\s*\()", std::regex::icase);
    // 匹配引脚连接: ".PIN(net)"
    std::regex conn_re(R"(\.[\w\[\]\:]+\s*\(([^)]+)\))");

    while (std::getline(ifs, line)) {
        // 简单去除前导空白
        size_t first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos) continue; // 空行

        // 1. 检查实例定义头
        std::smatch header_match;
        if (std::regex_search(line, header_match, header_re)) {
            current_inst_type = header_match[1].str();
            current_inst_name = header_match[2].str();

            if (points.contains(current_inst_name)) {
                points[current_inst_name].type = current_inst_type;
            }
            continue;
        }

        // 2. 如果不在有效实例块内，跳过
        if (current_inst_name.empty() || !points.contains(current_inst_name)) {
            if (line.find(");") != std::string::npos) {
                current_inst_name.clear();
            }
            continue;
        }

        // 3. 解析连接
        std::string::const_iterator search_start(line.cbegin());
        std::smatch conn_match;
        while (std::regex_search(search_start, line.cend(), conn_match, conn_re)) {
            std::string net_name = conn_match[1].str();

            // Trim net name
            net_name.erase(0, net_name.find_first_not_of(" \t"));
            net_name.erase(net_name.find_last_not_of(" \t") + 1);

            if (!net_name.empty()) {
                temp_point_nets[current_inst_name].emplace_back(net_name);
                temp_net_to_points[net_name].insert(current_inst_name);
            }
            search_start = conn_match.suffix().first;
        }

        // 4. 检查块结束
        if (line.find(");") != std::string::npos) {
            current_inst_name.clear();
        }
    }

    // 回填数据
    for (const auto& [net_name, point_list] : temp_net_to_points) {
        auto& net_points = net_map[net_name];
        net_points.reserve(point_list.size());
        for (const auto& point_name_str : point_list) {
            if (points.contains(point_name_str)) {
                net_points.insert(point_name_str);
            }
        }
    }

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

    // network_stats_cache_valid = false; // 已废弃，不需要
    return true;
}

// ==================== 资源与工具函数 ====================

std::string FPGAReader::normalizeResourceType(const std::string& type_name) const {
    if (type_name.empty()) return "";
    std::string upper_type = type_name;
    std::transform(upper_type.begin(), upper_type.end(), upper_type.begin(), ::toupper);
    if (upper_type.find("DSP") == 0) return "DSP";
    if (upper_type.find("RAMB") == 0 || upper_type.find("RAM") == 0) return "RAMB";
    return "";
}

void FPGAReader::initResourceCapacity() {
    tracked_resource_types = {"RAMB", "DSP"};
    for (int die = 0; die < num_dies; ++die) {
        die_type_capacity[die]["RAMB"] = MAX_RESOURCES_RAM_DIE;
        die_type_capacity[die]["DSP"] = MAX_RESOURCES_DSP_DIE;
    }
}

void FPGAReader::initResourceCount() {
    die_type_count.clear();
    for (const auto& type : tracked_resource_types) {
        die_type_count[type] = std::vector<int>(num_dies, 0);
    }

    for (const auto& [name, point] : points) {
        std::string normalized_type = normalizeResourceType(point.type);
        if (!normalized_type.empty() &&
            std::find(tracked_resource_types.begin(), tracked_resource_types.end(), normalized_type) != tracked_resource_types.end()) {
            if (point.die >= 0 && point.die < num_dies) {
                die_type_count[normalized_type][point.die]++;
            }
        }
    }
}

void FPGAReader::updateResourceCount(const std::string& point_name, int from_die, int to_die) {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return;

    const Point& point = point_it->second;
    std::string normalized_type = normalizeResourceType(point.type);

    if (normalized_type.empty()) return;

    // 只有在追踪列表里的资源才更新
    bool tracked = false;
    for(const auto& t : tracked_resource_types) if(t == normalized_type) tracked = true;
    if(!tracked) return;

    if (from_die >= 0 && from_die < num_dies) die_type_count[normalized_type][from_die]--;
    if (to_die >= 0 && to_die < num_dies) die_type_count[normalized_type][to_die]++;
}

bool FPGAReader::isValidMoveForCapacity(const std::string& point_name, int from_die, int to_die) const {
    if (from_die == to_die) return true;

    auto point_it = points.find(point_name);
    if (point_it == points.end()) return true;

    std::string normalized_type = normalizeResourceType(point_it->second.type);
    if (normalized_type.empty()) return true;

    bool tracked = false;
    for(const auto& t : tracked_resource_types) if(t == normalized_type) tracked = true;
    if(!tracked) return true;

    auto cap_it = die_type_capacity.find(to_die);
    if (cap_it == die_type_capacity.end()) return true;

    auto type_cap_it = cap_it->second.find(normalized_type);
    if (type_cap_it == cap_it->second.end()) return true;

    int current_usage = die_type_count.at(normalized_type)[to_die];
    int max_capacity = type_cap_it->second;

    return (current_usage < max_capacity);
}

// ==================== 惩罚参数 ====================

void FPGAReader::initPenaltyParameters() {
    base_penalty_coefficient = 1.0;
}

double FPGAReader::calculateDynamicPenaltyCoefficient(double utilization) const {
    if (utilization <= UTILIZATION_SOFT_THRESHOLD) {
        return 0.0;
    } else if (utilization <= UTILIZATION_HARD_THRESHOLD) {
        double normalized = (utilization - UTILIZATION_SOFT_THRESHOLD) /
                            (UTILIZATION_HARD_THRESHOLD - UTILIZATION_SOFT_THRESHOLD);
        return base_penalty_coefficient * std::pow(normalized, 2);
    } else {
        double excess = (utilization - UTILIZATION_HARD_THRESHOLD) /
                        (1.0 - UTILIZATION_HARD_THRESHOLD);
        return base_penalty_coefficient * (1.0 + 10.0 * std::exp(5.0 * excess));
    }
}

double FPGAReader::calculateDieUtilization(int die, const std::string& resource_type) const {
    if (die < 0 || die >= num_dies) return 0.0;

    auto type_it = die_type_count.find(resource_type);
    if (type_it == die_type_count.end()) return 0.0;

    int usage = type_it->second[die];
    int capacity = die_type_capacity.at(die).at(resource_type);

    return (capacity > 0) ? (static_cast<double>(usage) / capacity) : 0.0;
}

double FPGAReader::calculateUtilizationPenalty(int die, int from_die, const std::string& point_type) const {
    if (die < 0 || die >= num_dies) return 0.0;
    std::string normalized_type = normalizeResourceType(point_type);
    if (normalized_type.empty()) return 0.0;

    double ram_util = calculateDieUtilization(die, "RAMB");
    double dsp_util = calculateDieUtilization(die, "DSP");

    // 模拟移动后的增加
    if (die != from_die) {
        int capacity = 0;
        if (normalized_type == "RAMB") capacity = die_type_capacity.at(die).at("RAMB");
        else if (normalized_type == "DSP") capacity = die_type_capacity.at(die).at("DSP");

        if (capacity > 0) {
            if (normalized_type == "RAMB") ram_util = static_cast<double>(die_type_count.at("RAMB")[die] + 1) / capacity;
            else if (normalized_type == "DSP") dsp_util = static_cast<double>(die_type_count.at("DSP")[die] + 1) / capacity;
        }
    }

    double ram_penalty_coeff = calculateDynamicPenaltyCoefficient(ram_util);
    double dsp_penalty_coeff = calculateDynamicPenaltyCoefficient(dsp_util);
    double ram_penalty = ram_penalty_coeff * std::pow(ram_util, 2);
    double dsp_penalty = dsp_penalty_coeff * std::pow(dsp_util, 2);

    return std::sqrt(std::pow(ram_penalty, 2) + std::pow(dsp_penalty, 2));
}

double FPGAReader::calculatePenalizedGain(const std::string& point_name, int base_gain, int from_die, int to_die) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return GAIN_WEIGHT * base_gain;

    double weighted_gain = GAIN_WEIGHT * base_gain;
    double penalty = calculateUtilizationPenalty(to_die, from_die, point_it->second.type);
    return weighted_gain - penalty;
}

// ==================== 网络分布管理 (核心优化) ====================

// 初始化网络分布表
void FPGAReader::initNetDistribution() {
    net_distribution.clear();
    net_distribution.reserve(net_map.size());

    for (const auto& [net_name, point_set] : net_map) {
        std::vector<int> stats(num_dies, 0);
        for (const std::string& point_name : point_set) {
            if (auto it = points.find(point_name); it != points.end()) {
                int die = it->second.die;
                if (die >= 0 && die < num_dies) {
                    stats[die]++;
                }
            }
        }
        net_distribution[net_name] = stats;
    }
}

// 增量更新网络分布 (O(1))
void FPGAReader::updateNetDistributionForMove(const std::string& point_name, int from_die, int to_die) {
    if (from_die == to_die) return;

    const auto& point = points.at(point_name);
    for (const std::string& net_name : point.nets) {
        auto it = net_distribution.find(net_name);
        if (it != net_distribution.end()) {
            if (from_die >= 0 && from_die < num_dies) it->second[from_die]--;
            if (to_die >= 0 && to_die < num_dies) it->second[to_die]++;
        }
    }
}

// 极速查表
std::vector<int> FPGAReader::NetStats(const std::string& net_id) const {
    auto it = net_distribution.find(net_id);
    if (it != net_distribution.end()) {
        return it->second;
    }
    static std::vector<int> empty_stats(num_dies, 0);
    return empty_stats;
}

// ==================== FM 逻辑 ====================

double FPGAReader::Point_Gain(const std::string& point_name) const {
    auto [gain, target_die] = Point_Gain_Multi(point_name);
    return gain;
}

std::pair<double, int> FPGAReader::Point_Gain_Multi(const std::string& point_name) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return {0.0, -1};

    const Point& point = point_it->second;
    int current_die = point.die;
    double best_gain = -1e9;
    int best_target_die = current_die;

    if (num_dies == 2) {
        int target_die = 1 - current_die;
        int base_gain = calculateGainForMove(point_name, current_die, target_die);
        double penalized_gain = calculatePenalizedGain(point_name, base_gain, current_die, target_die);
        return {penalized_gain, target_die};
    }

    for (int target_die = 0; target_die < num_dies; ++target_die) {
        if (target_die == current_die) continue;
        int base_gain = calculateGainForMove(point_name, current_die, target_die);
        double penalized_gain = calculatePenalizedGain(point_name, base_gain, current_die, target_die);

        if (penalized_gain > best_gain) {
            best_gain = penalized_gain;
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
        // 直接读取 current distribution
        // 注意：计算Gain时，NetStats返回的是移动前的状态
        std::vector<int> stats = NetStats(net_name);

        int dies_before = 0;
        for (int c : stats) if (c > 0) dies_before++;
        bool before_cut = dies_before > 1;

        // 模拟移动后的状态
        stats[from_die]--;
        stats[to_die]++;

        int dies_after = 0;
        for (int c : stats) if (c > 0) dies_after++;
        bool after_cut = dies_after > 1;

        if (before_cut && !after_cut) gain += 1;
        else if (!before_cut && after_cut) gain -= 1;
    }
    return gain;
}

void FPGAReader::initSortedGains() {
    gain_map.clear();
    cached_gains.clear();
    cached_gains.reserve(points.size());

    // 预估一个合理的 bucket 大小

    for (const auto& [name, point] : points) {
        if (point.is_fixed) continue;
        double gain = Point_Gain(name);
        cached_gains[name] = gain;
        gain_map[gain].insert(name);
        gain_queue.push({gain, name});
    }
}

// 【关键优化】忽略高扇出网络 + 移除无效缓存清除
void FPGAReader::updateSortedGains(const std::string& moved_name) {
    const auto& moved_point = points.at(moved_name);
    std::unordered_set<std::string> affected_points;
    affected_points.reserve(500);

    for (const std::string& net_name : moved_point.nets) {
        auto net_it = net_map.find(net_name);
        if (net_it == net_map.end()) continue;

        // 忽略超大网络，更新邻居太慢且无意义
        if (net_it->second.size() > NET_FANOUT_THRESHOLD) continue;

        for (const std::string& point_name : net_it->second) {
            if (point_name != moved_name && !points[point_name].is_fixed) {
                affected_points.insert(point_name);
            }
        }
    }

    for (const std::string& point_name : affected_points) {
        double old_gain = cached_gains[point_name];

        // 从 map 中移除
        auto& range = gain_map[old_gain];
        range.erase(point_name);
        if (range.empty()) gain_map.erase(old_gain); // 如果桶空了，可以清理key

        // 计算新 gain
        auto [new_gain, target_die] = Point_Gain_Multi(point_name);
        cached_gains[point_name] = new_gain;
        gain_map[new_gain].insert(point_name);

        // 注意：gain_queue 包含过期数据，这在 FM 循环中通过 check cached_gains 处理
        // 但为了更频繁的更新，可以不在这里 push，完全依赖 lazy delete
    }
    // 移动点本身已在 FM 主循环处理
}

void FPGAReader::updateGainQueue() {
    std::priority_queue<std::pair<double, std::string>> empty_queue;
    std::swap(gain_queue, empty_queue);

    for (const auto& [gain_value, point_set] : gain_map) {
        if (point_set.empty()) continue;
        for (const auto& point_name : point_set) {
            // 双重确认
            if (!points[point_name].is_fixed) {
                gain_queue.push(std::make_pair(gain_value, point_name));
            }
        }
    }
}

void FPGAReader::FM() {
    if (points.empty()) return;

    auto start_time = std::chrono::high_resolution_clock::now();

    initResourceCapacity();
    initResourceCount();
    initPenaltyParameters();

    // 1. 初始化网络状态
    initNetDistribution();

    recordInitialStats();
    initSortedGains(); // 初始化 gain_map 和 gain_queue

    std::cout << "Starting FM loop..." << std::endl;

    while (true) {
        if (gain_queue.empty()) {
            updateGainQueue();
            if (gain_queue.empty()) break;
        }

        auto top_item = gain_queue.top();
        gain_queue.pop();
        std::string current_max_point_name = top_item.second;
        double queue_gain = top_item.first;

        // Lazy Deletion 检查:
        // 如果 queue 里的 gain 和当前 cached 不一致，说明是旧数据
        if (cached_gains.count(current_max_point_name)) {
            if (std::abs(cached_gains[current_max_point_name] - queue_gain) > 1e-6) {
                continue;
            }
        } else {
            // 点可能已经被移除了或者fixed了
            continue;
        }

        if (points[current_max_point_name].is_fixed) continue;

        int from_die = points[current_max_point_name].die;
        // 重新计算最佳目标 (虽然 gain 应该是一样的，但我们需要 target_die)
        auto [final_gain, target_die] = Point_Gain_Multi(current_max_point_name);

        // 验证一下计算出的 gain 是否还有效（防止微小浮点误差）
        // 这里主要用于判断是否真的值得移动

        if (!isValidMoveForCapacity(current_max_point_name, from_die, target_die)) {
            points[current_max_point_name].is_fixed = true;
            // 从系统中移除
            double old_gain = cached_gains[current_max_point_name];
            gain_map[old_gain].erase(current_max_point_name);
            cached_gains.erase(current_max_point_name);
            continue;
        }

        // 增益判断阈值
        int base_gain = calculateGainForMove(current_max_point_name, from_die, target_die);
        if (base_gain <= 0 || final_gain <= 0.0) {
            points[current_max_point_name].is_fixed = true;
            double old_gain = cached_gains[current_max_point_name];
            gain_map[old_gain].erase(current_max_point_name);
            cached_gains.erase(current_max_point_name);
            continue;
        }

        // --- 执行移动 ---

        // 1. 从数据结构移除旧状态
        double old_gain = cached_gains[current_max_point_name];
        gain_map[old_gain].erase(current_max_point_name);
        cached_gains.erase(current_max_point_name); // 当前点移出，之后不再移动它（或者也可以不移出，允许反复移动，视策略而定）
        // 如果希望一个点只能移动一次(locking)，保持上面这行。如果允许反复移动，则需重新插入。
        // 标准 FM 是 locking 的。
        points[current_max_point_name].is_fixed = true; // Locking

        // 2. 更新物理状态
        updateResourceCount(current_max_point_name, from_die, target_die);
        updateNetDistributionForMove(current_max_point_name, from_die, target_die);
        points[current_max_point_name].die = target_die;

        // 3. 更新邻居
        updateSortedGains(current_max_point_name);

        iteration_counter++;
        if (iteration_counter % 1000 == 0) {
            // std::cout << "Iteration: " << iteration_counter << std::endl;
        }
        if (iteration_counter >= UPDATE_FREQUENCY) {
            updateGainQueue(); // 重建堆，清理垃圾
            iteration_counter = 0;
        }
    }

    recordFinalStats();

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    execution_time = duration.count();
    std::cout << "FM finished in " << execution_time << "s" << std::endl;
}

// ==================== 统计与输出 ====================

// 注意：此函数现在是 O(N) 扫描所有 Net，因为不再维护全局 cached_cut_size
// 但 FM 循环中不需要调用它，只在最后调用，所以没问题。
int FPGAReader::calculateCutSize() const {
    int cut = 0;
    // 使用 net_distribution 计算
    for (const auto& [net, stats] : net_distribution) {
        int dies_occupied = 0;
        for(int c : stats) if(c > 0) dies_occupied++;
        if (dies_occupied > 1) cut++;
    }
    return cut;
}

int FPGAReader::calculateSLLBetweenDies(int die1, int die2) const {
    if (die1 == die2) return 0;
    int sll = 0;
    for (const auto& [net, stats] : net_distribution) {
        if (stats[die1] > 0 && stats[die2] > 0) {
            sll++;
        }
    }
    return sll;
}

// 下面这几个函数暂时保留，虽然 FM 优化版可能不再频繁调用它们
void FPGAReader::updateNetworkStatsCache() const { /* 已废弃 */ }
void FPGAReader::recordInitialStats() {
    initial_pairwise_cuts.clear();
    for (int i = 0; i < num_dies; ++i) {
        for (int j = i + 1; j < num_dies; ++j) {
            initial_pairwise_cuts[{i, j}] = calculateSLLBetweenDies(i, j);
        }
    }
}
void FPGAReader::recordFinalStats() {
    final_pairwise_cuts.clear();
    for (int i = 0; i < num_dies; ++i) {
        for (int j = i + 1; j < num_dies; ++j) {
            final_pairwise_cuts[{i, j}] = calculateSLLBetweenDies(i, j);
        }
    }
}

// 输出函数保持不变
void FPGAReader::print_Info() const {
    std::string output_file_name = "../FPGA_FM_output.place";
    std::ofstream ofs(output_file_name);
    if (!ofs) {
        std::cerr << "ERROR: Cannot create output placement file: " << output_file_name << "\n";
        return;
    }
    for (const auto& pair : points) {
        const Point& p = pair.second;
        ofs << p.name << " " << (int)p.y << " " << p.die;
        if (p.is_originally_fixed) ofs << " FIXED";
        ofs << "\n";
    }
    std::cout << "Successfully saved FM output to: " << output_file_name << "\n";
}

void FPGAReader::saveResultsToExcel(const std::string& excel_file) {
    bool file_exists = std::ifstream(excel_file).good();
    if (!file_exists) {
        std::ofstream ofs(excel_file);
        if (!ofs.is_open()) return;
        ofs << "File,Mode,"
            << "Initial_Die0-Die1,Initial_Die0-Die2,Initial_Die0-Die3,"
            << "Initial_Die1-Die2,Initial_Die1-Die3,Initial_Die2-Die3,"
            << "Final_Die0-Die1,Final_Die0-Die2,Final_Die0-Die3,"
            << "Final_Die1-Die2,Final_Die1-Die3,Final_Die2-Die3,"
            << "Initial_Total,Final_Total,Improvement,";
        for (int die = 0; die < 4; ++die) ofs << "RAMB_Die" << die << "(%),";
        for (int die = 0; die < 4; ++die) ofs << "DSP_Die" << die << "(%),";
        ofs << "Time(s)\n";
        ofs.close();
    }
    std::ofstream ofs(excel_file, std::ios::app);
    if (!ofs.is_open()) return;

    ofs << filename << "," << num_dies << "-die,";

    // 简化输出逻辑，统一处理 2-die 和 4-die
    auto get_cut = [&](const std::map<std::pair<int,int>, int>& cuts, int d1, int d2) {
        if(d1 > d2) std::swap(d1, d2);
        auto it = cuts.find({d1, d2});
        return (it != cuts.end()) ? it->second : 0;
    };

    // Initial Pairs
    ofs << get_cut(initial_pairwise_cuts, 0, 1) << "," << get_cut(initial_pairwise_cuts, 0, 2) << "," << get_cut(initial_pairwise_cuts, 0, 3) << ","
        << get_cut(initial_pairwise_cuts, 1, 2) << "," << get_cut(initial_pairwise_cuts, 1, 3) << "," << get_cut(initial_pairwise_cuts, 2, 3) << ",";

    // Final Pairs
    ofs << get_cut(final_pairwise_cuts, 0, 1) << "," << get_cut(final_pairwise_cuts, 0, 2) << "," << get_cut(final_pairwise_cuts, 0, 3) << ","
        << get_cut(final_pairwise_cuts, 1, 2) << "," << get_cut(final_pairwise_cuts, 1, 3) << "," << get_cut(final_pairwise_cuts, 2, 3) << ",";

    int initial_total = 0, final_total = 0;
    for (const auto& [k, v] : initial_pairwise_cuts) initial_total += v;
    for (const auto& [k, v] : final_pairwise_cuts) final_total += v;

    ofs << initial_total << "," << final_total << "," << (initial_total - final_total) << ",";

    // Resources
    auto write_res = [&](const std::string& type) {
        if (die_type_count.count(type)) {
            for (int die = 0; die < 4; ++die) {
                if (die < num_dies) {
                    int usage = die_type_count.at(type)[die];
                    int capacity = die_type_capacity.at(die).at(type);
                    double pct = (capacity > 0) ? (100.0 * usage / capacity) : 0.0;
                    ofs << std::fixed << std::setprecision(1) << pct << ",";
                } else ofs << "0.0,";
            }
        } else for(int i=0;i<4;i++) ofs << "0.0,";
    };

    write_res("RAMB");
    write_res("DSP");

    ofs << std::fixed << std::setprecision(2) << execution_time << "\n";
    ofs.close();
}

// 其余辅助函数 SLL check 等...
bool FPGAReader::isValidMoveForSLL(const std::string& point_name, int from_die, int to_die) const {
    // 简化：如果只关心总Cut最小化，通常不需要严格的 max SLL 限制，或者可以在 Penalty 里加。
    // 如果必须限制，逻辑比较复杂，需要差分更新。鉴于性能优化，暂时略过或返回true。
    return true;
}
int FPGAReader::calculateSLLChangeForDiePair(const std::string&, int, int, int, int) const { return 0; }