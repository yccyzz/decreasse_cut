#ifndef FPGA_READER_H
#define FPGA_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <regex>
#include <queue>
#include <cmath>
#include <filesystem>
#include <map>
#include <iomanip>

struct Point {
    std::string name;
    int y;
    int die; // 0, 1, 2, 3
    std::string type; // LUT, RAM, DSP, etc.
    bool is_fixed;
    bool is_originally_fixed; // 记录是否在输入文件中就是FIXED
    std::vector<std::string> nets;
};

class FPGAReader {
private:
    std::unordered_map<std::string, Point> points;
    std::unordered_map<std::string, std::unordered_set<std::string>> net_map;

    // --- 优化部分：永久状态维护，替代反复清除的缓存 ---
    // key: net_name, value: vector<int> (每个die上的计数)
    std::unordered_map<std::string, std::vector<int>> net_distribution;

    // 增益缓存
    std::unordered_map<std::string, double> cached_gains;
    // 增益桶：key=gain, value=set of point names
    std::map<double, std::unordered_set<std::string>> gain_map;
    // 优先队列，用于快速获取最大增益
    std::priority_queue<std::pair<double, std::string>> gain_queue;

    // 资源统计
    std::unordered_map<std::string, std::vector<int>> die_type_count;
    std::unordered_map<int, std::unordered_map<std::string, int>> die_type_capacity;
    std::vector<std::string> tracked_resource_types;

    // SLL 缓存 (die_pair -> count)
    mutable std::map<std::pair<int, int>, int> sll_cache;
    mutable int cached_cut_size = -1;
    mutable bool network_stats_cache_valid = false;

    // 统计数据
    std::map<std::pair<int, int>, int> initial_pairwise_cuts;
    std::map<std::pair<int, int>, int> final_pairwise_cuts;
    double execution_time = 0.0;

    // 参数
    int num_dies = 2;
    double max_gain = -1e9;
    std::string max_gain_point_name;
    int iteration_counter = 0;

    // 常量参数
    const int UPDATE_FREQUENCY = 5000;
    const int MAX_RESOURCES_RAM_DIE = 1500;
    const int MAX_RESOURCES_DSP_DIE = 1000;
    const int MAX_SLL_BETWEEN_DIES = 15000;
    const double GAIN_WEIGHT = 10.0;
    const double UTILIZATION_SOFT_THRESHOLD = 0.6;
    const double UTILIZATION_HARD_THRESHOLD = 0.85;

    // 【核心优化参数】忽略高扇出网络的阈值
    const int NET_FANOUT_THRESHOLD = 100;

    double base_penalty_coefficient = 1.0;

    // 私有辅助函数
    bool parsePlacement(const std::string& line, Point& p);
    void detectDieMode();

    // 资源相关
    std::string normalizeResourceType(const std::string& type_name) const;
    void initResourceCapacity();
    void initResourceCount();
    void updateResourceCount(const std::string& point_name, int from_die, int to_die);
    bool isValidMoveForCapacity(const std::string& point_name, int from_die, int to_die) const;

    // 惩罚计算
    void initPenaltyParameters();
    double calculateDynamicPenaltyCoefficient(double utilization) const;
    double calculateDieUtilization(int die, const std::string& resource_type) const;
    double calculateUtilizationPenalty(int die, int from_die, const std::string& point_type) const;

    // 增益计算核心
    // 【修改】直接查表，不再计算
    std::vector<int> NetStats(const std::string& net_id) const;

    double Point_Gain(const std::string& point_name) const;
    std::pair<double, int> Point_Gain_Multi(const std::string& point_name) const;
    int calculateGainForMove(const std::string& point_name, int from_die, int to_die) const;
    double calculatePenalizedGain(const std::string& point_name, int base_gain, int from_die, int to_die) const;

    // 优化后的状态更新函数
    void initNetDistribution();
    void updateNetDistributionForMove(const std::string& point_name, int from_die, int to_die);

    void initSortedGains();
    void updateSortedGains(const std::string& moved_name);
    void updateGainQueue();

    // SLL 相关
    void updateNetworkStatsCache() const;
    int calculateSLLBetweenDies(int die1, int die2) const;
    int calculateSLLChangeForDiePair(const std::string& point_name, int from_die, int to_die, int die1, int die2) const;
    bool isValidMoveForSLL(const std::string& point_name, int from_die, int to_die) const;

    // 统计记录
    void recordInitialStats();
    void recordFinalStats();

public:
    std::string filename;
    bool is_coarsing = false;

    FPGAReader() {}

    bool readPlaceFile(const std::string& file_name);
    bool readNetFile(std::string& file_name); // 修改了实现，去掉了 parseNet

    void FM();
    int calculateCutSize() const;
    void clear();
    void print_Info() const;
    void saveResultsToExcel(const std::string& excel_file);
};

#endif // FPGA_READER_H