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
#include <algorithm>
#include <set>
#include <iterator>

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
                .description("Single-tile code distance for X errors (Pauli weight of the minimum-weight logical X operator)")
                .required(true);
        parser.add_argument()
                .names({"-z", "--dz"})
                .description("Single-tile code distance for Z errors (Pauli weight of the minimum-weight logical Z operator)")
                .required(true);
        parser.add_argument()
                .names({"-t", "--dt"})
                .description("Code distance along the time dimension (number of surface code cycles for a round of error correction)")
                .required(true);
        parser.add_argument()
                .names({"-i", "--info"})
                .description("Information to be queried (no operation will be performed). Options: {instructions, plaquettes, grid}")
                .required(false);
        parser.add_argument()
                .names({"-o", "--operation"})
                .description("Surface code operation to be compiled. Options: {idle, prepz, prepx, measz, measx, extendx, extendz}")
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
        if ((parser.exists("i")) && (parser.exists("o"))) {
            std::cerr << "The -i and -o flags are mutually exclusive." << std::endl;
        }
        else if (parser.exists("i"))
        {
            std::string s = parser.get<std::string>("i");
            
            if (s == "instructions") {
                HardwareModel TI_model;
                TI_model.print_TI_ops();
            }

            else if (s == "plaquettes") {
                LogicalQubit lq(dx, dz, 0, 0, grid);
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

            // Single-qubit operations
            if ((s == "idle") || (s == "prepz") || (s == "prepx") || (s == "measz") || (s == "measx") || (s == "test")) {

                // Initialize logical qubit object using the grid
                LogicalQubit lq(dx, dz, 0, 0, grid);

                // Grab all of the initially occupied sites (to be used in printing)
                std::set<unsigned int> occupied_sites = lq.occupied_sites();

                // Initialize vector of hardware instructions
                std::vector<HW_Instruction> hw_master;

                // Initialize time tracker
                float time = 0;

                // Perform associated transversal operation
                if (s != "idle") {
                    lq.transversal_op(s, grid, hw_master, time);
                }

                // Append an idle operation if applicable 
                if ((s == "idle") || (s == "prepz") || (s == "prepx")) {
                    time = lq.idle(cycles, grid, hw_master, time);
                }

                // Placeholder input to help implement little test circuits
                if (s == "test") {
                    time = lq.test_circuits(grid, hw_master, time);
                }

                // Enforce validity of final instruction list 
                grid.enforce_hw_master_validity(hw_master);

                // Print hardware instructions
                print_hw_master(hw_master, occupied_sites, debug);

            }

            else if (s == "extendx") {

                // Two cases: odd vs. even dx determine whether we need one or two strips of qubits between 
                unsigned int extra_strip = 0;
                if (dx%2==0) {
                    extra_strip = 1;
                }

                // Construct grid with appropriate dimensions to hold two tiles side-by-side with a vertical strip of qubits in between
                unsigned int nrows = dz+1; 
                unsigned int ncols = 2*(dx+1) + extra_strip;
                GridManager grid_2(nrows, ncols);

                // Initialize logical qubit object using the grid
                LogicalQubit lq1(dx, dz, 0, 0, grid_2);

                // Initialize second logical qubit object to the right of the first
                LogicalQubit lq2(dx, dz, 0, dx+1+extra_strip, grid_2);

                // Initialize vector of hardware instructions
                std::vector<HW_Instruction> hw_master;

                // Prepare the physical qubits on lq2 in the X basis
                float time = 0;
                lq2.transversal_op("prepx", grid_2, hw_master, time);             

                // Create a merged qubit
                LogicalQubit lq = merge(lq1, lq2, grid_2);

                // std::cout << "Logical Qubit 1:" << std::endl;
                // lq1.print_stabilizers();

                // std::cout << "Logical Qubit 2:" << std::endl;
                // lq2.print_stabilizers();

                // std::cout << "Logical Qubit (merged):" << std::endl;
                // lq.print_stabilizers();

                // Grab all of the larger patch's occupied sites (to be used in printing)
                std::set<unsigned int> all_qsites = lq.occupied_sites();

                // Grab all of the qsites on the `strip' between lq1 and lq2
                std::set<unsigned int> strip = lq.get_strip(lq1, lq2);

                // Prepare qsites on the strip in the X basis
                float time_tmp = 0;
                HardwareModel TI_model;
                for (unsigned int site : strip) {
                    time_tmp = TI_model.add_init(site, time, 0, grid_2, hw_master);
                    time_tmp = TI_model.add_H(site, time_tmp, 1, grid_2, hw_master);
                }

                // Perform 'idle' operation on the merged qubit
                time = lq.idle(cycles, grid_2, hw_master, time);

                // Enforce validity of final instruction list 
                grid_2.enforce_hw_master_validity(hw_master);

                // Print hardware instructions
                print_hw_master(hw_master, all_qsites, debug);
            }

            else if (s == "extendz") {

                // Two cases: odd vs. even dz determine whether we need one or two strips of qubits between 
                unsigned int extra_strip = 0;
                if (dz%2==0) {
                    extra_strip = 1;
                }

                // Construct grid with appropriate dimensions to hold two tiles top and bottom with a horizontal strip of qubits in between
                unsigned int nrows = 2*(dz+1) + extra_strip; 
                unsigned int ncols = dx+1;
                GridManager grid_2(nrows, ncols);

                // Initialize logical qubit object using the grid
                LogicalQubit lq1(dx, dz, 0, 0, grid_2);

                // Initialize second logical qubit object to the bottom of the first
                LogicalQubit lq2(dx, dz, dz+1+extra_strip, 0, grid_2);

                // Initialize vector of hardware instructions
                std::vector<HW_Instruction> hw_master;

                // Prepare the physical qubits on lq2 in the Z basis
                float time = 0;
                lq2.transversal_op("prepz", grid_2, hw_master, time);

                // Create a merged qubit
                LogicalQubit lq = merge(lq1, lq2, grid_2);

                // std::cout << "Logical Qubit 1:" << std::endl;
                // lq1.print_stabilizers();

                // std::cout << "Logical Qubit 2:" << std::endl;
                // lq2.print_stabilizers();

                // std::cout << "Logical Qubit (merged):" << std::endl;
                // lq.print_stabilizers();

                // Grab all of the larger patch's occupied sites (to be used in printing)
                std::set<unsigned int> all_qsites = lq.occupied_sites();

                // Grab all of the qsites on the `strip' between lq1 and lq2
                std::set<unsigned int> strip = lq.get_strip(lq1, lq2);

                // Prepare qsites on the strip in the Z basis
                HardwareModel TI_model;
                for (unsigned int site : strip) {
                    TI_model.add_init(site, time, 0, grid_2, hw_master);
                }

                // Perform 'idle' operation on the merged qubit
                time = lq2.idle(cycles, grid_2, hw_master, time);

                // Enforce validity of final instruction list 
                grid_2.enforce_hw_master_validity(hw_master);

                // Print hardware instructions
                print_hw_master(hw_master, all_qsites, debug);
            }

            else {std::cerr << "No valid operation selected. Options: {idle, prepz, prepx, measz, measx, extendx, extendz}" << std::endl;}
        }

        return 0;
    }
}