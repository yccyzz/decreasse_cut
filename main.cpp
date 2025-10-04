#include "fpga_reader.h"

int main() {
    FPGAReader reader;
    std::string placement_file = "C:\\Users\\Lenovo\\Desktop\\plresults\\FPGA12.pl";
    std::string netlist_file = "C:\\Users\\Lenovo\\Desktop\\ispd2016_flexshelf\\FPGA12\\design.osv";
    std::filesystem::path p(placement_file);
    reader.filename = p.filename().string();

    std::cout << reader.filename << std::endl;
    if (!reader.readPlaceFile(placement_file) || !reader.readNetFile(netlist_file)) {
        std::cerr << "File read error\n";
        return 1;
    }
    reader.FM();
    reader.saveResultsToExcel("C:\\Users\\Lenovo\\Desktop\\2.4.csv");
//    reader.print_Info();
   return 0;
}

