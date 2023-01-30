#include <TISCC/pipeline.hpp>

#include <argparse/argparse.h>
#include <iostream>

namespace TISCC 
{
    int run_tiscc(
        int argc, const char* argv[],
        std::istream& in_stream,
        std::ostream& out_stream,
        std::ostream& err_stream)
    {
        std::string prog_name{argv[0]};
        argparse::ArgumentParser parser(prog_name, "Trapped-Ion Surface Code Compiler");
        parser.add_argument()
                .names({"-d", "--dummy"})
                .description("Dummy argument to test argparse")
                .required(false);
        parser.enable_help();
        auto err = parser.parse(argc, argv);
        if (err)
        {
            err_stream << err << std::endl;
            parser.print_help();
            return -1;
        }
        if (parser.exists("help"))
        {
            parser.print_help();
            return 0;
        }
        if (parser.exists("d"))
        {
            std::cout << "Hello world!" << std::endl;
        }

        return 0;
    }
}