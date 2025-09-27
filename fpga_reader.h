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
#include <map>
#include <regex>
#include <algorithm>
#include <chrono>
#include <queue>
#include <array>
#include <algorithm>
#include <filesystem>

struct Point {
    std::string name;
    int y;
    int die;
    bool is_fixed;
    std::vector<std::string> nets;
};

class FPGAReader {
private:
    static const int UPDATE_FREQUENCY = 100;
    static const int MAX_SLL_BETWEEN_DIES = 14400;  // SLL硬约束

    std::unordered_map<std::string, Point> points;
    std::unordered_map<std::string, std::unordered_set<std::string>> net_map;

    // 缓存系统
    mutable std::unordered_map<std::string, std::vector<int>> net_cache;
    std::unordered_map<std::string, int> cached_gains;

    // SLL缓存系统
    mutable std::map<std::pair<int,int>, int> sll_cache;
    mutable int cached_cut_size = -1;
    mutable bool network_stats_cache_valid = false;

    // Gain管理
    std::map<int, std::unordered_set<std::string>> gain_map;
    std::priority_queue<std::pair<int, std::string>> gain_queue;

    int num_dies = 2;
    int max_gain = -INT32_MAX;
    std::string max_gain_point_name;
    int iteration_counter = 0;

public:
    bool readPlaceFile(const std::string& file_name);
    bool readNetFile(std::string& file_name);
    void detectDieMode();
    void clear();

    // 解析函数
    bool parsePlacement(const std::string& line, Point& p);
    void parseNet(std::string& content);

    // 核心算法
    void FM();
    void initSortedGains();
    void updateSortedGains(const std::string& moved_name);
    void updateGainQueue();

    // Cut和SLL计算（优化版本）
    void updateNetworkStatsCache() const;
    void invalidateNetworkStatsCache();
    int calculateCutSize() const;
    int calculateSLLBetweenDies(int die1, int die2) const;

    // SLL约束检查
    bool isValidMoveForSLL(const std::string& point_name, int from_die, int to_die) const;
    int calculateSLLChangeForDiePair(const std::string& point_name, int from_die, int to_die, int die1, int die2) const;

    // Gain计算
    std::vector<int> NetStats(const std::string& net_id) const;
    int Point_Gain(const std::string& point_name) const;
    std::pair<int, int> Point_Gain_Multi(const std::string& point_name) const;
    int calculateGainForMove(const std::string& point_name, int from_die, int to_die) const;

    // 输出和调试
    void print_Info() const;
    void printSLLStats() const;
};

#endif // FPGA_READER_H