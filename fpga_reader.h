#ifndef FPGA_READER_H
#define FPGA_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <regex>
#include <array>
#include <algorithm>
#include <climits>
#include <chrono>
#include <iomanip>
#include <cctype>
#include<queue>
#include <numeric>

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <array>
#include <map>
struct Point {
    std::string name;
    int y = 0;
    int die = 0;
    bool is_fixed = false;
    std::vector<std::string> nets;

    Point() = default;
    Point(const Point&) = default;
    Point(Point&&) = default;
    Point& operator=(const Point&) = default;
    Point& operator=(Point&&) = default;
    Point(const std::string& n) : name(n) {}
};

class FPGAReader {
public:
    bool readPlaceFile(const std::string& file_name);
    bool readNetFile(std::string& file_name);
    bool parsePlacement(const std::string& line, Point& p);
    void parseNet(std::string& content);
    void print_Info() const;
    void FM();
    void clear();
    const std::unordered_map<std::string, Point>& getPoints() const { return points; }
    size_t getPointCount() const { return points.size(); }
    int calculateCutSize() const;

private:
    std::unordered_map<std::string, Point> points;
    std::unordered_map<std::string, std::unordered_set<std::string>> net_map;
    std::unordered_map<int, std::unordered_set<std::string>> gain_map;
    std::unordered_map<std::string, int> cached_gains;
    int max_gain = -INT32_MAX;
    std::string max_gain_point_name;
    mutable std::unordered_map<std::string, std::array<int, 2>> net_cache;

    void updateGainQueue();
    int iteration_counter = 0;
    const int UPDATE_FREQUENCY = 100; // 或其他合适的值

    std::priority_queue<std::pair<int, std::string>> gain_queue;

    void initSortedGains();
    void updateSortedGains(const std::string& moved_name);
    std::array<int, 2> NetStats(const std::string& net_id) const;
    int Point_Gain(const std::string& point_name) const;
};
#endif // FPGA_READER_H