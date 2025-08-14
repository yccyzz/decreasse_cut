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
};

class FPGAReader {
public:
    bool readPlaceFile(const std::string& file_name);
    bool readNetFile(const std::string& file_name);
    bool parsePlacement(const std::string& line, Point& p);
    void parseNet(const std::string& content);
    void buildIndexes();
    void print_Info() const;
    void FM();
    const std::vector<Point>& getPoints() const { return points; }
    size_t getPointCount() const { return points.size(); }
    int calculateCutSize() const;

private:
    std::vector<Point> points;
    std::unordered_map<std::string, size_t> inst_map;
    std::unordered_map<std::string, std::vector<size_t>> net_map;

    mutable std::unordered_map<std::string, std::array<int, 2>> net_cache;
    mutable std::vector<int> cached_gains;
    mutable bool gains_initialized = false;

    void initial_Gain();
    std::array<int, 2> NetStats(const std::string& net_id) const;
    int Point_Gain(size_t point_idx) const;
    void updateGain(size_t moved_idx);
    void updateAllGains() const;
};

#endif // FPGA_READER_H