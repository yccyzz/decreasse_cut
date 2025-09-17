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
#include <regex>
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
    std::unordered_map<std::string, Point> points;
    std::unordered_map<std::string, std::unordered_set<std::string>> net_map;

    // 增益相关数据结构
    std::map<int, std::unordered_set<std::string>> gain_map;
    std::priority_queue<std::pair<int, std::string>> gain_queue;
    std::unordered_map<std::string, int> cached_gains;

    // 缓存和优化  
    mutable std::unordered_map<std::string, std::vector<int>> net_cache;

    // FM算法状态
    int max_gain = -INT32_MAX;
    std::string max_gain_point_name;
    int iteration_counter = 0;
    static const int UPDATE_FREQUENCY = 1000;  // 队列更新频率

    // Die模式支持
    int num_dies = 2;  // 默认2-die模式，可以是2或4

public:
    bool readPlaceFile(const std::string& file_name);
    bool readNetFile(std::string& file_name);
    void clear();
    void FM();
    void print_Info() const;
    int calculateCutSize() const;

private:
    bool parsePlacement(const std::string& line, Point& p);
    void parseNet(std::string& content);
    void detectDieMode();  // 检测当前是2-die还是4-die模式

    // 增益计算函数
    std::vector<int> NetStats(const std::string& net_id) const;
    int Point_Gain(const std::string& point_name) const;  // 保持向后兼容的2-die版本
    std::pair<int, int> Point_Gain_Multi(const std::string& point_name) const;  // 多die版本，返回(gain, target_die)
    int calculateGainForMove(const std::string& point_name, int from_die, int to_die) const;

    // 增益管理函数
    void initSortedGains();
    void updateSortedGains(const std::string& moved_name);
    void updateGainQueue();
};

#endif // FPGA_READER_H