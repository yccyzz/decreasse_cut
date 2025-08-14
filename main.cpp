#include "fpga_reader.h"

int main() {
    FPGAReader reader;
    std::string placement_file = "C:\\Users\\Lenovo\\Downloads\\FPGA01.pl";
    std::string netlist_file = "C:\\Users\\Lenovo\\Downloads\\design.osv";

    if (!reader.readPlaceFile(placement_file) || !reader.readNetFile(netlist_file)) {
        std::cerr << "File read error\n";
        return 1;
    }

    reader.buildIndexes();
    reader.FM();
    reader.print_Info();

    return 0;
}

