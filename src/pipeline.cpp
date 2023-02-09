#include <TISCC/pipeline.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/logicalqubit.hpp>
#include <TISCC/plaquette.hpp>
#include <TISCC/hardwaremodel.hpp>
#include <TISCC/instruction.hpp>

#include <argparse/argparse.h>
#include <iostream>
#include <cassert>

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
                .names({"-i", "--info"})
                .description("Information to be queried (no operation will be performed). Options: {``instructions'', ``plaquettes'', ``grid''}")
                .required(false);
        parser.add_argument()
                .names({"-o", "--operation"})
                .description("Surface code operation to be compiled")
                .required(false);
        parser.add_argument()
                .names({"-d", "--debug"})
                .description("Provide extra output for the purpose of debugging")
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

        // Initializing grid using GridManager object
        unsigned int nrows = dz+1; 
        unsigned int ncols = dx+1;
        GridManager grid(nrows, ncols);

        // Query information (no operation will be performed)
        if (parser.exists("i")) 
        {
            std::string s = parser.get<std::string>("i");
            if (s == "instructions") {
                HardwareModel TI_model;
                TI_model.print_TI_ops();
            }

            else if (s == "plaquettes") {
                LogicalQubit lq(dx, dz, grid);
                lq.print_stabilizers();
            }

            else if (s == "grid") {
                grid.print_grid();
            }

            return 0;
        }

        // If debugging output is requested, create a bool to be passed into the functions below
        bool debug = false;
        if (parser.exists("d")) {
            debug = true;
        }

        // Operation-dependent logic
        if (parser.exists("o"))
        {
            std::string s = parser.get<std::string>("o");
            if (s == "idle") {
                // Initialize logical qubit using the grid
                LogicalQubit lq(dx, dz, grid);

                // Perform 'idle' operation
                lq.idle(cycles, grid, debug);

            }
        }
        else {std::cout << "Grid constructed but no valid operation selected. Quitting." << std::endl;}

        return 0;
    }
}