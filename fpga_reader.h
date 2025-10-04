#ifndef FPGA_READER_H
#define FPGA_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <regex>
#include <algorithm>
#include <queue>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <filesystem>
struct Point {
    std::string name;
    int y = 0;
    int die = 0;
    bool is_fixed = false;
    std::string type;  // 实例类型：RAMB36E2, RAMB18E2, DSP48E2 等
    std::vector<std::string> nets;
};

class FPGAReader {
public:
    // 常量定义
    static constexpr int MAX_SLL_BETWEEN_DIES = 800;
    static constexpr int MAX_RESOURCES_RAM_DIE = 432;  // RAM硬约束
    static constexpr int MAX_RESOURCES_DSP_DIE = 192;  // DSP硬约束
    static constexpr int UPDATE_FREQUENCY = 100;
    static constexpr double GAIN_WEIGHT = 100.0;
    static constexpr double UTILIZATION_SOFT_THRESHOLD = 0.6;
    static constexpr double UTILIZATION_HARD_THRESHOLD = 0.9;

    // 基础文件读取
    bool readPlaceFile(const std::string& file_name);
    bool readNetFile(std::string& file_name);

    // FM算法
    void FM();

    // 输出
    void print_Info() const;
    void saveResultsToExcel(const std::string& excel_file);

    std::string filename;

private:
    // 数据存储
    std::unordered_map<std::string, Point> points;
    std::unordered_map<std::string, std::unordered_set<std::string>> net_map;

    // 缓存
    mutable std::unordered_map<std::string, std::vector<int>> net_cache;
    mutable std::unordered_map<std::string, double> cached_gains;
    mutable std::map<std::pair<int, int>, int> sll_cache;
    mutable int cached_cut_size = -1;
    mutable bool network_stats_cache_valid = false;

    // 资源缓存
    std::map<std::string, std::vector<int>> die_type_count;  // "RAMB" or "DSP" -> [die0_count, die1_count, ...]
    std::map<int, std::map<std::string, int>> die_type_capacity;  // die -> {"RAMB": 432, "DSP": 192}
    std::vector<std::string> tracked_resource_types;

    // FM算法状态
    std::map<double, std::unordered_set<std::string>> gain_map;
    std::priority_queue<std::pair<double, std::string>> gain_queue;
    int iteration_counter = 0;
    double max_gain = -1e9;
    std::string max_gain_point_name;
    double base_penalty_coefficient = 1.0;

    // 配置
    int num_dies = 2;

    // 统计数据
    std::map<std::pair<int, int>, int> initial_pairwise_cuts;
    std::map<std::pair<int, int>, int> final_pairwise_cuts;
  double execution_time = 0.0;

    // 解析函数
    bool parsePlacement(const std::string& line, Point& p);
    void parseNet(std::string& content);

    // Die模式检测
    void detectDieMode();

    // 清理
    void clear();

    // 资源统计
    void initResourceCapacity();
    void initResourceCount();
    void updateResourceCount(const std::string& point_name, int from_die, int to_die);
    bool isValidMoveForCapacity(const std::string& point_name, int from_die, int to_die) const;
    std::string normalizeResourceType(const std::string& type_name) const;

    // 惩罚计算
    void initPenaltyParameters();
    double calculateDynamicPenaltyCoefficient(double utilization) const;
    double calculateDieUtilization(int die, const std::string& resource_type) const;
    double calculateUtilizationPenalty(int die, int from_die, const std::string& point_type) const;
    double calculatePenalizedGain(const std::string& point_name, int base_gain,
                                  int from_die, int to_die) const;

    // 网络统计
    void updateNetworkStatsCache() const;
    void invalidateNetworkStatsCache();
    int calculateCutSize() const;
    int calculateSLLBetweenDies(int die1, int die2) const;
    std::vector<int> NetStats(const std::string& net_id) const;

    // SLL约束
    int calculateSLLChangeForDiePair(const std::string& point_name,
                                     int from_die, int to_die,
                                     int die1, int die2) const;
    bool isValidMoveForSLL(const std::string& point_name,
                           int from_die, int to_die) const;

    // Gain计算
    double Point_Gain(const std::string& point_name) const;
    std::pair<double, int> Point_Gain_Multi(const std::string& point_name) const;
    int calculateGainForMove(const std::string& point_name,
                             int from_die, int to_die) const;

    // FM算法辅助
    void initSortedGains();
    void updateSortedGains(const std::string& moved_name);
    void updateGainQueue();

    // 统计辅助
    void recordInitialStats();
    void recordFinalStats();
    void initializeCapacity();
};

#endif // FPGA_READER_H