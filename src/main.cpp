#include <TISCC/pipeline.hpp>

#include <iostream>

int main(int argc, const char* argv[]) {
    return TISCC::run_tiscc(argc, argv, std::cin, std::cout, std::cerr);

    return 0;
}