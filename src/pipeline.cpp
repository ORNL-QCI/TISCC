#include <TISCC/pipeline.hpp>
#include <TISCC/gridmanager.hpp>
#include <TISCC/logicalqubit.hpp>
#include <TISCC/plaquette.hpp>
#include <TISCC/hardwaremodel.hpp>
#include <TISCC/instruction.hpp>

#include <argparse/argparse.h>
#include <iostream>
#include <cassert>
#include <limits>

namespace TISCC 
{
    // Once operations are compiled into hardware instructions, this function is used to output them
    void print_hw_master(const std::vector<HW_Instruction>& hw_master, const std::set<unsigned int>& occupied_sites, bool debug) 
    {
        // I/O settings
        int W = 15;
        std::cout << std::setprecision(1);
        std::cout << std::setiosflags(std::ios::fixed);

        // Dump qsites 
        for (unsigned int site : occupied_sites) {
            std::cout << std::setw(W) << -1.0;
            std::cout << std::setw(W) << "Qubit_at";
            std::cout << std::setw(W) << site;
            std::cout << std::endl;
        }

        // Output HW instructions to file
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        for (const HW_Instruction& instruction : hw_master) {
            std::cout << std::setw(W) << instruction.get_time();
            std::cout << std::setw(W) << instruction.get_name();
            if (instruction.get_site2() != uint_max) {
                std::cout << std::setw(W) << ((std::to_string(instruction.get_site1()) += ",") += std::to_string(instruction.get_site2()));
            }
            else {
                std::cout << std::setw(W) << instruction.get_site1();
            }
            if (debug) {
                std::cout << std::setw(W) << instruction.get_step();
                std::cout << std::setw(W) << instruction.get_q1();
                std::cout << std::setw(W) << instruction.get_q2();
                std::cout << std::setw(W) << instruction.get_shape();
                std::cout << std::setw(W) << instruction.get_type();
            }
            std::cout << std::endl;
        }
    }

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
        // TODO: Implement resource counting and make available via this flag
        // parser.add_argument()
        //         .names({"-r", "--resources"})
        //         .description("Provide resource estimates")
        //         .required(false);
        
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

            else {
                std::cerr << "No valid query selected. Options: {instructions, plaquettes, grid}" << std::endl;
            }
            
            return 0;
        }

        // If debugging output is requested, create a bool pass into functions below
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

                // Grab all of the initially occupied sites
                std::set<unsigned int> occupied_sites = lq.occupied_sites();

                // Initialize vector of hardware instructions
                std::vector<HW_Instruction> hw_master;

                // Perform 'idle' operation
                lq.idle(cycles, grid, hw_master);

                // Make sure all of the same qsites are occupied
                assert(occupied_sites == lq.occupied_sites());

                // Test validity of final instruction list
                grid.check_hw_master_validity(hw_master);

                // Print hardware instructions
                print_hw_master(hw_master, occupied_sites, debug);
            }
            else {std::cerr << "No valid operation selected. Options: {idle}" << std::endl;}
        }

        return 0;
    }
}