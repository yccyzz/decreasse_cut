#include "fpga_reader.h"

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
        parsePlacement(line, p);
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
    detectDieMode();
    return true;
}

bool FPGAReader::parsePlacement(const std::string& line, Point& p) {
    std::istringstream ss(line);
    std::string name, x, z, last;
    int y,die;
    if(!is_coarsing){
        if (!(ss >> name >> x >> y >> z)) return false;
        p.name = name;
        p.y = y;
        p.is_fixed = false;
        p.is_originally_fixed = false;

        if (ss >> last && last == "FIXED") {
            p.is_fixed = true;
            p.is_originally_fixed = true;
        }
    } else{
        if (!(ss >> name >> y >> die)) return false;
        p.name = name;
        p.die = die;
        if (die < 0 || die > 3) {
            std::cout << die << std::endl;
        }
        p.y = die * 120.0 + 60.0;
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

    if (dies_used.size() == 2 && dies_used.count(0) && dies_used.count(1)) {
        num_dies = 2;
    } else {
        num_dies = 4;
    }
}

void FPGAReader::clear() {
    points.clear();
    net_map.clear();
    net_cache.clear();
    cached_gains.clear();

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


void FPGAReader::parseNet(std::string& content) {
    std::replace(content.begin(), content.end(), '\n', ' ');
    std::regex inst_re(R"(([\w]+)\s*([\w.]*)\s*\(([\s\S]*?)\);)", std::regex::icase);
    std::regex net_re(R"(\.[\w\.]+\s*\(([\w\.]*)\))");
    std::unordered_map<std::string, std::vector<std::string>> temp_point_nets;
    std::unordered_map<std::string, std::unordered_set<std::string>> temp_net_to_points;
    temp_point_nets.reserve(points.size());
    temp_net_to_points.reserve(points.size());
    net_map.reserve(50000);
    for (std::sregex_iterator it(content.begin(), content.end(), inst_re), end;
         it != end; ++it) {
        const auto& m = *it;
        std::string inst_type = m[1].str();
        std::string inst_name = m[2].str();
        std::string conn = m[3].str();

        if (!points.contains(inst_name)) {
            continue;
        }

        points[inst_name].type = inst_type;

        for (auto net_it = std::sregex_iterator(conn.begin(), conn.end(), net_re);
             net_it != std::sregex_iterator(); ++net_it) {
            std::string net_name = (*net_it)[1].str();
            temp_point_nets[inst_name].emplace_back(net_name);
            temp_net_to_points[net_name].insert(inst_name);
        }
    }
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
        }else {
            std::cout << point.name << " not found in temp_point_nets." << std::endl;
        }
    }
    invalidateNetworkStatsCache();
}

// ==================== 工具函数 ====================

std::string FPGAReader::normalizeResourceType(const std::string& type_name) const {
    if (type_name.empty()) return "";
    std::string upper_type = type_name;
    std::transform(upper_type.begin(), upper_type.end(), upper_type.begin(), ::toupper);
    if (upper_type.find("DSP") == 0) {
        return "DSP";
    } else if (upper_type.find("RAMB") == 0 || upper_type.find("RAM") == 0) {
        return "RAMB";
    }
    return "";
}

// ==================== 资源容量管理 ====================

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
            std::find(tracked_resource_types.begin(),
                      tracked_resource_types.end(),
                      normalized_type) != tracked_resource_types.end()) {

            if (point.die >= 0 && point.die < num_dies) {
                die_type_count[normalized_type][point.die]++;
            }
        }
    }
}

void FPGAReader::updateResourceCount(const std::string& point_name,
                                     int from_die, int to_die) {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return;

    const Point& point = point_it->second;
    std::string normalized_type = normalizeResourceType(point.type);

    if (normalized_type.empty() ||
        std::find(tracked_resource_types.begin(),
                  tracked_resource_types.end(),
                  normalized_type) == tracked_resource_types.end()) {
        return;
    }

    if (from_die >= 0 && from_die < num_dies) {
        die_type_count[normalized_type][from_die]--;
    }
    if (to_die >= 0 && to_die < num_dies) {
        die_type_count[normalized_type][to_die]++;
    }
}

bool FPGAReader::isValidMoveForCapacity(const std::string& point_name,
                                        int from_die, int to_die) const {
    if (from_die == to_die) return true;

    auto point_it = points.find(point_name);
    if (point_it == points.end()) return true;

    const Point& point = point_it->second;
    std::string normalized_type = normalizeResourceType(point.type);

    if (normalized_type.empty() ||
        std::find(tracked_resource_types.begin(),
                  tracked_resource_types.end(),
                  normalized_type) == tracked_resource_types.end()) {
        return true;
    }

    auto capacity_it = die_type_capacity.find(to_die);
    if (capacity_it == die_type_capacity.end()) return true;

    auto type_capacity_it = capacity_it->second.find(normalized_type);
    if (type_capacity_it == capacity_it->second.end()) return true;

    int current_usage = die_type_count.at(normalized_type)[to_die];
    int max_capacity = type_capacity_it->second;

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
        // 二次增长
        double normalized = (utilization - UTILIZATION_SOFT_THRESHOLD) /
                            (UTILIZATION_HARD_THRESHOLD - UTILIZATION_SOFT_THRESHOLD);
        return base_penalty_coefficient * std::pow(normalized, 2);
    } else {
        // 指数增长
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

    auto capacity_it = die_type_capacity.find(die);
    if (capacity_it == die_type_capacity.end()) return 0.0;

    auto res_capacity_it = capacity_it->second.find(resource_type);
    if (res_capacity_it == capacity_it->second.end()) return 0.0;

    int capacity = res_capacity_it->second;

    return (capacity > 0) ? (static_cast<double>(usage) / capacity) : 0.0;
}

double FPGAReader::calculateUtilizationPenalty(int die, int from_die,
                                               const std::string& point_type) const {
    if (die < 0 || die >= num_dies) return 0.0;

    std::string normalized_type = normalizeResourceType(point_type);
    if (normalized_type.empty()) return 0.0;

    double ram_util = calculateDieUtilization(die, "RAMB");
    double dsp_util = calculateDieUtilization(die, "DSP");

    // 如果是移动到目标die，模拟添加资源
    if (die != from_die) {
        auto capacity_it = die_type_capacity.find(die);
        if (capacity_it != die_type_capacity.end()) {
            auto res_capacity_it = capacity_it->second.find(normalized_type);
            if (res_capacity_it != capacity_it->second.end()) {
                int capacity = res_capacity_it->second;
                if (capacity > 0) {
                    if (normalized_type == "RAMB") {
                        ram_util = static_cast<double>(die_type_count.at("RAMB")[die] + 1) / capacity;
                    } else if (normalized_type == "DSP") {
                        dsp_util = static_cast<double>(die_type_count.at("DSP")[die] + 1) / capacity;
                    }
                }
            }
        }
    }

    double ram_penalty_coeff = calculateDynamicPenaltyCoefficient(ram_util);
    double dsp_penalty_coeff = calculateDynamicPenaltyCoefficient(dsp_util);

    double ram_penalty = ram_penalty_coeff * std::pow(ram_util, 2);
    double dsp_penalty = dsp_penalty_coeff * std::pow(dsp_util, 2);

    double total_penalty = std::sqrt(std::pow(ram_penalty, 2) + std::pow(dsp_penalty, 2));

    return total_penalty;
}

double FPGAReader::calculatePenalizedGain(const std::string& point_name,
                                          int base_gain, int from_die, int to_die) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return GAIN_WEIGHT * base_gain;

    const Point& point = point_it->second;
    double weighted_gain = GAIN_WEIGHT * base_gain;
    double penalty = calculateUtilizationPenalty(to_die, from_die, point.type);

    return weighted_gain - penalty;
}

// ==================== 网络统计计算 ====================

void FPGAReader::updateNetworkStatsCache() const {
    if (network_stats_cache_valid) return;

    sll_cache.clear();
    cached_cut_size = 0;

    for (int i = 0; i < num_dies; ++i) {
        for (int j = i + 1; j < num_dies; ++j) {
            sll_cache[{i, j}] = 0;
        }
    }

    for (const auto& [net_name, point_set] : net_map) {
        std::set<int> dies_in_net;

        for (const std::string& point_name : point_set) {
            if (auto it = points.find(point_name); it != points.end()) {
                dies_in_net.insert(it->second.die);
            }
        }

        if (dies_in_net.size() > 1) {
            cached_cut_size++;

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
    net_cache.clear();
}

int FPGAReader::calculateCutSize() const {
    updateNetworkStatsCache();
    return cached_cut_size;
}

int FPGAReader::calculateSLLBetweenDies(int die1, int die2) const {
    if (die1 == die2) return 0;

    if (die1 > die2) std::swap(die1, die2);

    updateNetworkStatsCache();

    auto it = sll_cache.find({die1, die2});
    return (it != sll_cache.end()) ? it->second : 0;
}

// ==================== SLL约束检查 ====================

int FPGAReader::calculateSLLChangeForDiePair(const std::string& point_name,
                                             int from_die, int to_die,
                                             int die1, int die2) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return 0;

    const Point& point = point_it->second;
    int sll_change = 0;

    for (const std::string& net_name : point.nets) {
        auto stats = NetStats(net_name);

        bool before_connects_die1 = (stats[die1] > 0);
        bool before_connects_die2 = (stats[die2] > 0);
        bool before_sll = before_connects_die1 && before_connects_die2;

        stats[from_die]--;
        stats[to_die]++;

        bool after_connects_die1 = (stats[die1] > 0);
        bool after_connects_die2 = (stats[die2] > 0);
        bool after_sll = after_connects_die1 && after_connects_die2;

        if (!before_sll && after_sll) {
            sll_change += 1;
        } else if (before_sll && !after_sll) {
            sll_change -= 1;
        }
    }

    return sll_change;
}

bool FPGAReader::isValidMoveForSLL(const std::string& point_name,
                                   int from_die, int to_die) const {
    if (from_die == to_die) return true;

    if (num_dies == 2) {
        int current_sll = calculateSLLBetweenDies(0, 1);
        int sll_change = calculateSLLChangeForDiePair(point_name, from_die, to_die, 0, 1);
        return (current_sll + sll_change <= MAX_SLL_BETWEEN_DIES);
    }

    for (int i = 0; i < num_dies; ++i) {
        for (int j = i + 1; j < num_dies; ++j) {
            int current_sll = calculateSLLBetweenDies(i, j);
            int sll_change = calculateSLLChangeForDiePair(point_name, from_die, to_die, i, j);

            if (current_sll + sll_change > MAX_SLL_BETWEEN_DIES) {
                return false;
            }
        }
    }

    return true;
}

// ==================== 统计辅助函数 ====================

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

// ==================== FM算法实现 ====================

void FPGAReader::FM() {
    if (points.empty()) return;

    auto start_time = std::chrono::high_resolution_clock::now();

    initResourceCapacity();
    initResourceCount();
    initPenaltyParameters();
    recordInitialStats();

    initSortedGains();

    while (!gain_map.empty()) {

        if (gain_queue.empty()) {
            updateGainQueue();
            if (gain_queue.empty()) {
                break;
            }
        }

        auto top_item = gain_queue.top();
        gain_queue.pop();

        std::string current_max_point_name = top_item.second;

        if (points[current_max_point_name].is_fixed) {
            continue;
        }

        // 重新计算当前点的base_gain和final_gain
        int target_die;
        int base_gain;
        double final_gain;
        int from_die = points[current_max_point_name].die;

        if (num_dies == 2) {
            target_die = 1 - from_die;
            base_gain = calculateGainForMove(current_max_point_name, from_die, target_die);
            final_gain = calculatePenalizedGain(current_max_point_name, base_gain,
                                                from_die, target_die);
        } else {
            auto [best_final_gain, best_target] = Point_Gain_Multi(current_max_point_name);
            target_die = best_target;
            base_gain = calculateGainForMove(current_max_point_name, from_die, target_die);
            final_gain = best_final_gain;
        }

        // 检查资源硬约束
        if (!isValidMoveForCapacity(current_max_point_name, from_die, target_die)) {
            points[current_max_point_name].is_fixed = true;
            if (cached_gains.count(current_max_point_name)) {
                double old_gain = cached_gains[current_max_point_name];
                gain_map[old_gain].erase(current_max_point_name);
            }
            continue;
        }

        // 结束条件判断：
        // 1. base_gain <= 0：移动无法减少cut size → fixed
        // 2. base_gain > 0 但 final_gain <= 0：需要移动但代价太大 → fixed
        if (base_gain <= 0 || final_gain <= 0) {
            points[current_max_point_name].is_fixed = true;
            if (cached_gains.count(current_max_point_name)) {
                double old_gain = cached_gains[current_max_point_name];
                gain_map[old_gain].erase(current_max_point_name);
            }
            continue;
        }

        // 只有当 base_gain > 0 且 final_gain > 0 时才执行移动
        // 从gain_map中移除该点
        if (cached_gains.count(current_max_point_name)) {
            double old_gain = cached_gains[current_max_point_name];
            gain_map[old_gain].erase(current_max_point_name);
        }

        // 主动更新资源计数
        updateResourceCount(current_max_point_name, from_die, target_die);

        // 执行移动（不标记为fixed，允许后续再移动）
        points[current_max_point_name].die = target_die;

        // 更新受影响点的增益
        updateSortedGains(current_max_point_name);

        iteration_counter++;
        if (iteration_counter >= UPDATE_FREQUENCY) {
            updateGainQueue();
            iteration_counter = 0;
        }
    }

    recordFinalStats();

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    execution_time = duration.count();
}

void FPGAReader::initSortedGains() {
    gain_map.clear();
    net_cache.clear();
    cached_gains.clear();
    cached_gains.reserve(points.size());

    for (const auto& [name, point] : points) {
        double gain = Point_Gain(name);
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

    std::unordered_set<std::string> affected_points;
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
        double old_gain = cached_gains[point_name];
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
        double penalized_gain = calculatePenalizedGain(point_name, base_gain,
                                                       current_die, target_die);
        return {penalized_gain, target_die};
    }

    for (int target_die = 0; target_die < num_dies; ++target_die) {
        if (target_die == current_die) continue;

        int base_gain = calculateGainForMove(point_name, current_die, target_die);
        double penalized_gain = calculatePenalizedGain(point_name, base_gain,
                                                       current_die, target_die);

        if (penalized_gain > best_gain) {
            best_gain = penalized_gain;
            best_target_die = target_die;
        }
    }

    return {best_gain, best_target_die};
}

int FPGAReader::calculateGainForMove(const std::string& point_name,
                                     int from_die, int to_die) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return 0;

    const Point& point = point_it->second;
    int gain = 0;

    for (const std::string& net_name : point.nets) {
        auto stats = NetStats(net_name);

        int dies_before = 0;
        for (int i = 0; i < num_dies; ++i) {
            if (stats[i] > 0) dies_before++;
        }
        bool before_cut = dies_before > 1;

        stats[from_die]--;
        stats[to_die]++;

        int dies_after = 0;
        for (int i = 0; i < num_dies; ++i) {
            if (stats[i] > 0) dies_after++;
        }
        bool after_cut = dies_after > 1;

        if (before_cut && !after_cut) {
            gain += 1;
        } else if (!before_cut && after_cut) {
            gain -= 1;
        }
    }

    return gain;
}

void FPGAReader::updateGainQueue() {
    std::priority_queue<std::pair<double, std::string>> empty_queue;
    std::swap(gain_queue, empty_queue);

    for (const auto& [gain_value, point_set] : gain_map) {
        if (point_set.empty()) {
            continue;
        }
        for (const auto& point_name : point_set) {
            if (!points[point_name].is_fixed) {
                gain_queue.push(std::make_pair(gain_value, point_name));
            }
        }
    }
}
void FPGAReader::print_Info() const {
    std::string output_file_name = "../FPGA_FM_output.place";
    std::ofstream ofs(output_file_name);

    if (!ofs) {
        std::cerr << "ERROR: Cannot create output placement file: " << output_file_name << "\n";
        return;
    }

    for (const auto& pair : points) {
        const Point& p = pair.second;

        // 【修改】新格式: name y die [FIXED]
        // x 坐标取整，die 直接使用 Point.die
        ofs << p.name << " " << (int)p.y << " " << p.die;

        // 仅当点在输入时就为 FIXED 时，才输出 FIXED
        if (p.is_originally_fixed) {
            ofs << " FIXED";
        }

        ofs << "\n";
    }

    std::cout << "Successfully saved FM output to: " << output_file_name << "\n";
}

void FPGAReader::saveResultsToExcel(const std::string& excel_file) {
    bool file_exists = std::ifstream(excel_file).good();

    if (!file_exists) {
        std::ofstream ofs(excel_file);
        if (!ofs.is_open()) {
            std::cerr << "ERROR: Cannot create file: " << excel_file << std::endl;
            return;
        }

        ofs << "File,Mode,";
        ofs << "Initial_Die0-Die1,Initial_Die0-Die2,Initial_Die0-Die3,"
            << "Initial_Die1-Die2,Initial_Die1-Die3,Initial_Die2-Die3,"
            << "Final_Die0-Die1,Final_Die0-Die2,Final_Die0-Die3,"
            << "Final_Die1-Die2,Final_Die1-Die3,Final_Die2-Die3,";
        ofs << "Initial_Total,Final_Total,Improvement,";

        for (int die = 0; die < 4; ++die) {
            ofs << "RAMB_Die" << die << "(%),";
        }

        for (int die = 0; die < 4; ++die) {
            ofs << "DSP_Die" << die << "(%),";
        }

        ofs << "Time(s)\n";

        ofs.close();
    }

    std::ofstream ofs(excel_file, std::ios::app);

    if (!ofs.is_open()) {
        std::cerr << "ERROR: Cannot open file for appending: " << excel_file << std::endl;
        return;
    }

    ofs << filename << "," << num_dies << "-die,";

    if (num_dies == 2) {
        auto initial_it = initial_pairwise_cuts.find({0, 1});
        auto final_it = final_pairwise_cuts.find({0, 1});

        if (initial_it == initial_pairwise_cuts.end() || final_it == final_pairwise_cuts.end()) {
            std::cerr << "ERROR: Missing pairwise cut data" << std::endl;
            ofs.close();
            return;
        }

        ofs << initial_it->second << ",0,0,0,0,0,";
        ofs << final_it->second << ",0,0,0,0,0,";
    } else {
        for (int i = 0; i < num_dies; ++i) {
            for (int j = i + 1; j < num_dies; ++j) {
                auto it = initial_pairwise_cuts.find({i, j});
                ofs << (it != initial_pairwise_cuts.end() ? it->second : 0) << ",";
            }
        }
        for (int i = 0; i < num_dies; ++i) {
            for (int j = i + 1; j < num_dies; ++j) {
                auto it = final_pairwise_cuts.find({i, j});
                ofs << (it != final_pairwise_cuts.end() ? it->second : 0) << ",";
            }
        }
    }

    int initial_total = 0, final_total = 0;
    for (const auto& [k, v] : initial_pairwise_cuts) initial_total += v;
    for (const auto& [k, v] : final_pairwise_cuts) final_total += v;

    ofs << initial_total << "," << final_total << ","
        << (initial_total - final_total) << ",";

    if (die_type_count.find("RAMB") != die_type_count.end()) {
        for (int die = 0; die < 4; ++die) {
            if (die < num_dies && die < die_type_count.at("RAMB").size()) {
                int usage = die_type_count.at("RAMB")[die];
                int capacity = die_type_capacity.at(die).at("RAMB");
                double percentage = (capacity > 0) ? (100.0 * usage / capacity) : 0.0;
                ofs << std::fixed << std::setprecision(1) << percentage << ",";
            } else {
                ofs << "0.0,";
            }
        }
    } else {
        for (int die = 0; die < 4; ++die) {
            ofs << "0.0,";
        }
    }

    if (die_type_count.find("DSP") != die_type_count.end()) {
        for (int die = 0; die < 4; ++die) {
            if (die < num_dies && die < die_type_count.at("DSP").size()) {
                int usage = die_type_count.at("DSP")[die];
                int capacity = die_type_capacity.at(die).at("DSP");
                double percentage = (capacity > 0) ? (100.0 * usage / capacity) : 0.0;
                ofs << std::fixed << std::setprecision(1) << percentage << ",";
            } else {
                ofs << "0.0,";
            }
        }
    } else {
        for (int die = 0; die < 4; ++die) {
            ofs << "0.0,";
        }
    }

    ofs << std::fixed << std::setprecision(2) << execution_time << "\n";

    ofs.close();
}