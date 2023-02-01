#include <TISCC/pipeline.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/logicalqubit.hpp>

#include <argparse/argparse.h>
#include <iostream>

namespace TISCC 
{
    // Parses command line arguments and runs through the code pipeline
    int run_tiscc(
        int argc, const char* argv[],
        std::istream& in_stream,
        std::ostream& out_stream,
        std::ostream& err_stream)
    {
        // Setting up command line argument parser
        std::string prog_name{argv[0]};
        argparse::ArgumentParser parser(prog_name, "Trapped-Ion Surface Code Compiler");
        parser.add_argument()
                .names({"-x", "--dx"})
                .description("Code distance for X errors (Pauli weight of the minimum-weight logical X operator)")
                .required(true);
        parser.add_argument()
                .names({"-z", "--dz"})
                .description("Code distance for Z errors (Pauli weight of the minimum-weight logical Z operator)")
                .required(true);
        parser.add_argument()
                .names({"-t", "--dt"})
                .description("Code distance along the time dimension (number of surface code cycles)")
                .required(true);
        parser.add_argument()
                .names({"-o", "--operation"})
                .description("Surface code operation to be compiled")
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

        // Extracting required arguments from command line input
        unsigned int dx = parser.get<unsigned int>("x");
        unsigned int dz = parser.get<unsigned int>("z");
        unsigned int cycles = parser.get<unsigned int>("t");

        // Initializing grid on the heap using GridManager object
        unsigned int nrows = dz+1; 
        unsigned int ncols = dx+1;
        GridManager a(nrows, ncols);
        a.print_grid();
        // unsigned int b = 45;
        // std::cout << a[b] << " " << a.grid(b) << " " << *(a.grid(b)) << std::endl;
        // std::cout << a.get_row(b) << " " << a.get_col(b) << " " << a.get_idx(b) << " " << a.val_from_coords(a.get_row(b), a.get_col(b), a.get_idx(b)) << std::endl;

        // Operation-dependent logic
        if (parser.exists("o"))
        {
            std::string s = parser.get<std::string>("o");
            if (s == "idle") {
                // Initialize logical qubit using the grid
                LogicalQubit lq(dx, dz, a);

                // Perform 'idle' operation
                lq.idle(cycles);
            }
        }
        else {std::cout << "Grid constructed and checks performed, but no valid operation selected. Quitting." << std::endl;}

        return 0;
    }
}