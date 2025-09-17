#include "fpga_reader.h"

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    FPGAReader reader;
    std::string placement_file = "C:\\Users\\Lenovo\\Desktop\\plresults\\FPGA02.pl";
    std::string netlist_file = "C:\\Users\\Lenovo\\Desktop\\ispd2016_flexshelf\\FPGA02\\designwithout.osv";
    std::filesystem::path p(placement_file);
    std::string filename = p.filename().string();

    std::cout << filename << std::endl;
    if (!reader.readPlaceFile(placement_file) || !reader.readNetFile(netlist_file)) {
        std::cerr << "File read error\n";
        return 1;
    }
//    auto fm_start = std::chrono::high_resolution_clock::now();
    reader.FM();
//    auto fm_end = std::chrono::high_resolution_clock::now();
//    std::cout << "FM Time" << std::chrono::duration<double>(fm_end - fm_start).count() << std::endl;
//    reader.print_Info();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout<< "using time : "<<std::chrono::duration<double>(end - start).count() << "s"<<std::endl;
    return 0;
}

