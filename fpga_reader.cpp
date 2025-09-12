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

    return true;
}

void FPGAReader::clear() {
    points.clear();
    net_map.clear();
    net_cache.clear();
    cached_gains.clear();
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
    std::regex net_re(R"(\.\w+\(([\w.]+)\))");

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

    std::cout << "Parsed " << net_map.size() << " nets connecting " << points.size() << " points" << std::endl;
}

int FPGAReader::calculateCutSize() const {
    int cut_size = 0;

    // 遍历所有网络，统计跨分区的网络数量
    for (const auto& [net_name, point_set] : net_map) {
        bool has_die0 = false, has_die1 = false;

        // 检查网络中的点分布在哪些分区
        for (const std::string& point_name : point_set) {
            if (auto it = points.find(point_name); it != points.end()) {
                if (it->second.die == 0) has_die0 = true;
                else if (it->second.die == 1) has_die1 = true;

                // 如果同时存在于两个分区，这是一个割边
                if (has_die0 && has_die1) {
                    cut_size++;
                    break;
                }
            }
        }
    }
    return cut_size;
}

void FPGAReader::FM() {
    if (points.empty()) return;
    auto cutsize_start = std::chrono::high_resolution_clock::now();
    int initial_cut_size = calculateCutSize();
    auto cutsize_end = std::chrono::high_resolution_clock::now();
    std::cout << "CutSize time: " << std::chrono::duration<double>(cutsize_end - cutsize_start).count() << std::endl;
    std::cout << "Initial cut size: " << initial_cut_size << std::endl;

    auto init_gain_start = std::chrono::high_resolution_clock::now();
    // 初始化排序的增益表
    initSortedGains();
    auto init_gain_end = std::chrono::high_resolution_clock::now();
    std::cout << "Init gain time: " << std::chrono::duration<double>(init_gain_end - init_gain_start).count() << std::endl;


    // FM Pass - 严格只移动正增益点
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

        // 检查是否为正增益，如果不是则退出
        if (current_max_gain < 0) {
            break;
        }

        // 检查该点是否已被固定（可能在队列更新间隔内被处理过）
        if (points[current_max_point_name].is_fixed) {
            continue;
        }

        // 验证该点在gain_map中仍然存在（确保数据一致性）
        if (gain_map[current_max_gain].find(current_max_point_name) == gain_map[current_max_gain].end()) {
            continue; // 该点可能已经被更新，跳过
        }

        // 从gain_map中移除该点
        gain_map[current_max_gain].erase(current_max_point_name);


        // 执行移动
        points[current_max_point_name].die = 1 - points[current_max_point_name].die;
        points[current_max_point_name].is_fixed = true;


        // 更新受影响点的增益（这会更新gain_map）
        updateSortedGains(current_max_point_name);

        // 定期更新队列以保持gain_queue和gain_map同步
        iteration_counter++;
        if (iteration_counter >= UPDATE_FREQUENCY) {
            updateGainQueue();
            iteration_counter = 0;
        }


        // 更新成员变量（用于其他函数可能的访问）
//        max_gain = current_max_gain;
//        max_gain_point_name = current_max_point_name;

//        // 获取最大增益的点（严格 > 0）
//        gain_map[max_gain].erase(max_gain_point_name);
//        max_gain = -INT32_MAX;
//        // 执行移动
//        points[max_gain_point_name].die = 1 - points[max_gain_point_name].die;
//        points[max_gain_point_name].is_fixed = true;


//        auto update_gain_start = std::chrono::high_resolution_clock::now();
//        max_gain = -INT32_MAX;
//        for (const auto& gain : gain_map) {
//            if (gain.second.empty()) {
//                continue;
//            }
//            if (gain.first > max_gain) {
//                max_gain = gain.first;
//                max_gain_point_name = *gain.second.begin();
//            }
//        }
//        auto update_gain_end = std::chrono::high_resolution_clock::now();
//        std::cout << "Update gain time: " << std::chrono::duration<double>(update_gain_end - update_gain_start).count() << std::endl;
    }

    std::cout << "Final cut size: " << calculateCutSize() << std::endl;
    std::cout << "Improvement: " << (initial_cut_size - calculateCutSize()) << std::endl;
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
            break;
        }
        int new_gain = Point_Gain(point_name);
        cached_gains[point_name] = new_gain;
        auto& new_gain_range = gain_map[new_gain];
        new_gain_range.insert(point_name);
    }
    cached_gains.erase(moved_name);
}

std::array<int, 2> FPGAReader::NetStats(const std::string& net_id) const {
    if (auto it = net_cache.find(net_id); it != net_cache.end()) {
        return it->second;
    }

    std::array<int, 2> stats = {0, 0};
    if (auto net_it = net_map.find(net_id); net_it != net_map.end()) {
        for (const std::string& point_name : net_it->second) {
            if (auto point_it = points.find(point_name); point_it != points.end()) {
                int die = point_it->second.die;
                if (die >= 0 && die <= 1) {
                    stats[die]++;
                }
            }
        }
    }

    net_cache[net_id] = stats;
    return stats;
}

int FPGAReader::Point_Gain(const std::string& point_name) const {
    auto point_it = points.find(point_name);
    if (point_it == points.end()) return 0;

    const Point& point = point_it->second;
    int gain = 0;
    int current_die = point.die;
    int other_die = 1 - current_die;

    for (const std::string& net_name : point.nets) {
        auto stats = NetStats(net_name);
        bool before_cut = (stats[0] > 0) && (stats[1] > 0);

        stats[current_die]--;
        stats[other_die]++;
        bool after_cut = (stats[0] > 0) && (stats[1] > 0);

        if (before_cut && !after_cut) {
            gain += 1;  // 消除割边
        } else if (!before_cut && after_cut) {
            gain -= 1;  // 增加割边
        }
    }
    return gain;
}

void FPGAReader::print_Info() const {
    for (int die = 0; die < 2; ++die) {
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