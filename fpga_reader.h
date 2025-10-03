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
#include <iomanip>
#include <cctype>
#include <chrono>
#include <queue>
#include <array>
#include <filesystem>
#include <cmath>

constexpr int UPDATE_FREQUENCY = 100;
constexpr double UTILIZATION_SOFT_THRESHOLD = 0.60;
constexpr double UTILIZATION_HARD_THRESHOLD = 0.90;
constexpr double GAIN_WEIGHT = 100.0;

struct Point {
    std::string name;
    int y;
    int die;
    bool is_fixed;
    std::vector<std::string> nets;
    std::string type;
};

class FPGAReader {
public:
    bool readPlaceFile(const std::string& file_name);
    bool readNetFile(std::string& file_name);
    void FM();
    void print_Info() const;
    void setInputFileName(const std::string& filename);
    void saveResultsToExcel(const std::string& excel_file);
private:
    std::unordered_map<std::string, Point> points;
    std::unordered_map<std::string, std::unordered_set<std::string>> net_map;
    std::string excel_output_path = "C:\\Users\\Lenovo\\Desktop\\2.4.xlsx";
    mutable std::unordered_map<std::string, std::vector<int>> net_cache;
    mutable std::unordered_map<std::string, double> cached_gains;
    std::map<double, std::unordered_set<std::string>> gain_map;
    std::priority_queue<std::pair<double, std::string>> gain_queue;

    mutable int cached_cut_size;
    mutable bool network_stats_cache_valid;

    std::map<int, std::map<std::string, int>> die_type_capacity;
    std::map<std::string, std::vector<int>> die_type_count;
    std::vector<std::string> tracked_resource_types;

    double base_penalty_coefficient;

    int num_dies;
    double max_gain;
    std::string max_gain_point_name;
    int iteration_counter;
    std::string current_input_file;
    int initial_cut_size;
    int final_cut_size;
    std::map<std::pair<int, int>, int> initial_pairwise_cuts;
    std::map<std::pair<int, int>, int> final_pairwise_cuts;
    double execution_time;

    void clear();
    bool parsePlacement(const std::string& line, Point& p);
    void parseNet(std::string& content);
    void detectDieMode();

    void initSortedGains();
    void updateSortedGains(const std::string& moved_name);
    double Point_Gain(const std::string& point_name) const;
    std::pair<double, int> Point_Gain_Multi(const std::string& point_name) const;
    int calculateGainForMove(const std::string& point_name, int from_die, int to_die) const;
    void updateGainQueue();

    std::vector<int> NetStats(const std::string& net_id) const;
    void invalidateNetworkStatsCache();

    void initResourceCapacity();
    void initResourceCount();
    void updateResourceCount(const std::string& point_name, int from_die, int to_die);
    bool isValidMoveForCapacity(const std::string& point_name, int from_die, int to_die) const;
    std::string normalizeResourceType(const std::string& type_name) const;

    void initPenaltyParameters();
    double calculateDynamicPenaltyCoefficient(double utilization) const;
    double calculateUtilizationPenalty(int die, int from_die, const std::string& point_type) const;
    double calculateDieUtilization(int die, const std::string& resource_type) const;
    double calculatePenalizedGain(const std::string& point_name, int base_gain, int from_die, int to_die) const;
    void printResourceUtilization() const;// 在头文件中添加函数声明（如果你有fpga_reader.h）
    std::map<std::pair<int, int>, int> calculatePairwiseCutSizes() const;
    void printPairwiseCutSizes() const;
};

#endif // FPGA_READER_H