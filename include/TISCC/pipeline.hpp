#ifndef TISCC_PIPELINE_HPP
#define TISCC_PIPELINE_HPP

#include <TISCC/instruction.hpp>

#include <iostream>
#include <vector>
#include <set>

namespace TISCC {
    int run_tiscc(
        int argc, const char* argv[],
        std::istream& in_stream,
        std::ostream& out_stream,
        std::ostream& err_stream);
}

#endif //TISCC_PIPELINE_HPP