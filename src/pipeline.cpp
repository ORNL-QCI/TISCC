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
#include <string>

namespace TISCC 
{
    // Once operations are compiled into hardware instructions, this function is used to output them
    void print_hw_master(std::ostream& output, const std::vector<HW_Instruction>& hw_master, const std::set<unsigned int>& occupied_sites, bool debug) 
    {
        // I/O settings
        int W = 15;
        output << std::setprecision(1);
        output << std::setiosflags(std::ios::fixed);

        // Dump qsites 
        for (unsigned int site : occupied_sites) {
            output << std::setw(W) << -1.0;
            output << std::setw(W) << "Qubit_at";
            output << std::setw(W) << site;
            output << std::endl;
        }

        // Output HW instructions to file
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        for (const HW_Instruction& instruction : hw_master) {
            output << std::setw(W) << instruction.get_time();
            output << std::setw(W) << instruction.get_name();
            if (instruction.get_site2() != uint_max) {
                output << std::setw(W) << ((std::to_string(instruction.get_site1()) += ",") += std::to_string(instruction.get_site2()));
            }
            else {
                output << std::setw(W) << instruction.get_site1();
            }
            if (debug) {
                output << std::setw(W) << instruction.get_step();
                output << std::setw(W) << instruction.get_q1();
                output << std::setw(W) << instruction.get_q2();
                output << std::setw(W) << instruction.get_shape();
                output << std::setw(W) << instruction.get_type();
            }
            output << std::endl;
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
                .description("Information to be queried. Options: {instructions, plaquettes, grid, parity}")
                .required(false);
        parser.add_argument()
                .names({"-o", "--operation"})
                .description("Surface code operation to be compiled. Options: {idle, prepz, prepx, measz, measx, hadamard, inject_y, inject_t, flip_patch, swap_left, rotation, extension, contraction, merge, split, bellprep, bellmeas}")
                .required(false);
        parser.add_argument()
                .names({"-s", "--tile_spec"})
                .description("Size of grid. Options: {single (default), double-vert, double-horiz}")
                .required(false);      
        parser.add_argument()
                .names({"-d", "--debug"})
                .description("Provide extra output for the purpose of debugging")
                .required(false);
        parser.add_argument()
                .names({"-p", "--printgates"})
                .description("Print the time-resolved gate sequence for the operation given")
                .required(false);
        parser.add_argument()
                .names({"-r", "--resources"})
                .description("Provide resource analysis for the operation given")
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

        // If debugging output is requested, create a bool to pass into functions below
        bool debug = false;
        if (parser.exists("d")) {
            debug = true;
        }

        // Extracting required arguments from command line input
        unsigned int dx = parser.get<unsigned int>("x");
        unsigned int dz = parser.get<unsigned int>("z");
        unsigned int cycles = parser.get<unsigned int>("t");

        // Designate the number rows and columns for a single tile
        // Note: We add one or two buffer strips of qubits depending on whether the code distance is even or odd
        unsigned int nrows = dz+1+!(dz%2); 
        unsigned int ncols = dx+1+!(dx%2);

        // Declare needed variables for one- and two-tile operations
        TISCC::LogicalQubit* lq;
        TISCC::LogicalQubit* lq1;
        TISCC::LogicalQubit* lq2;
        TISCC::GridManager* grid;
        std::set<unsigned int> strip;

        // Extract tile_spec
        std::string tile_spec = "single";
        if (parser.exists("s")) {
            tile_spec = parser.get<std::string>("s");
        }

        // Case-dependent grid & lq initialization
        if (tile_spec == "single") {

            // Initialize grid
            grid = new TISCC::GridManager(nrows, ncols);

            // Initialize logical qubit on  grid
            lq = new TISCC::LogicalQubit(dx, dz, 0, 0, *grid);
        }

        else if (tile_spec == "double-vert") {

            // Construct grid with appropriate dimensions to hold two tiles top and bottom with a horizontal strip of qubits in between
            grid = new TISCC::GridManager(2*nrows, ncols);

            // Initialize logical qubit object using the grid
            lq1 = new TISCC::LogicalQubit(dx, dz, 0, 0, *grid);

            // Initialize second logical qubit object to the bottom of the first
            lq2 = new TISCC::LogicalQubit(dx, dz, nrows, 0, *grid);

            // Create a merged qubit
            lq = merge(*lq1, *lq2, *grid);

            // Grab all of the qsites on the `strip' between lq1 and lq2
            strip = lq->get_strip(*lq1, *lq2);
        }

        else if (tile_spec == "double-horiz") {

            // Construct grid with appropriate dimensions to hold two tiles top and bottom with a horizontal strip of qubits in between
            grid = new TISCC::GridManager(nrows, 2*ncols);

            // Initialize logical qubit object using the grid
            lq1 = new TISCC::LogicalQubit(dx, dz, 0, 0, *grid);

            // Initialize second logical qubit object to the bottom of the first
            lq2 = new TISCC::LogicalQubit(dx, dz, 0, ncols, *grid);

            // Create a merged qubit
            lq = merge(*lq1, *lq2, *grid);

            // Grab all of the qsites on the `strip' between lq1 and lq2
            strip = lq->get_strip(*lq1, *lq2);
        }

        else {
            std::cerr << "Invalid tile_spec." << std::endl;
            abort();
        }

        if (parser.exists("i")) {
            std::string s = parser.get<std::string>("i");
            
            if (s == "instructions") {
                HardwareModel TI_model;
                TI_model.print_TI_ops();
            }

            else if (s == "plaquettes") {
                lq->print_stabilizers();
            }
          
            else if (s == "grid") {
                grid->print_qsite_mapping();
                std::cout << std::endl;
                std::vector<std::string> ascii_grid = grid->ascii_grid(false);
                grid->print_grid(ascii_grid);
                std::cout << std::endl;
                ascii_grid = grid->ascii_grid(true);
                grid->print_grid(ascii_grid);
            }

            else if (s == "parity") {
                lq->print_parity_check_matrix();
            }

            else {
                std::cerr << "No valid query selected. Options: {instructions, plaquettes, grid, parity}" << std::endl;
            }
            
        }

        // Operation-dependent logic
        if (parser.exists("o")) {

            // Initialize Hardware Model
            HardwareModel TI_model;

            // Initialize vector of hardware instructions
            std::vector<HW_Instruction> hw_master;

            // Initialize time tracker
            double time = 0;

            std::string s = parser.get<std::string>("o");

            // Single-patch operations
            if ((s == "idle") || (s == "prepz") || (s == "prepx") || (s == "measz") || (s == "measx") || (s == "inject_y") || (s == "inject_t") || (s == "flip_patch") || (s == "hadamard") || (s == "swap_left") || (s == "rotation")) {

                // Perform associated transversal operation
                if ((s == "prepz") || (s == "prepx") || (s == "measz") || (s == "measx") || (s == "hadamard")) {
                    lq->transversal_op(s, *grid, hw_master, time);
                }

                else if ((s == "inject_y") || (s == "inject_t")) {
                    char state_label = s.substr(7)[0];
                    lq->inject_state(state_label, *grid, hw_master, time);
                }

                else if ((s == "flip_patch") || (s == "rotation")) {
                    time = lq->flip_patch(*grid, hw_master, time, false);
                }

                else if ((s == "rotation")) {
                    // Will need to allocate a grid with one extra column to the right, actually
                    // In this block, will need to extend to the right and then contract
                }

                else if ((s == "swap_left") || (s == "rotation")) {
                    time = lq->swap_left(*grid, hw_master, time);
                }

                // Append an idle operation if applicable 
                if ((s == "idle") || (s == "prepz") || (s == "prepx") || (s == "inject_y") || (s == "inject_t") || (s == "flip_patch") || (s == "swap_left") || (s == "rotation") || (s == "hadamard")) {
                    time = lq->idle(cycles, *grid, hw_master, time);
                }

            }

            // Patch extensions
            else if (s == "extension") {

                if (tile_spec == "single") {std::cerr << "extension: invalid tile_spec given." << std::endl; abort();}

                // Prepare lq2 in approp. basis depending on direction of extension
                else if (tile_spec == "double-horiz") {
                    lq2->transversal_op("prepx", *grid, hw_master, time);
                }

                else if (tile_spec == "double-vert") {
                    lq2->transversal_op("prepz", *grid, hw_master, time);
                }

                // Prepare strip qubits
                double time_tmp = 0;
                for (unsigned int site : strip) {
                    time_tmp = TI_model.add_init(site, time, 0, *grid, hw_master);

                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time_tmp, 1, *grid, hw_master);
                }

                // Perform idle operation
                time = lq->idle(cycles, *grid, hw_master, time);

            }

            // Patch contractions
            else if (s == "contraction") {

                if (tile_spec == "single") {std::cerr << "contraction: invalid tile_spec given." << std::endl; abort();}

                // Measure lq2 in approp. basis depending on direction of contraction
                else if (tile_spec == "double-vert") {
                    lq2->transversal_op("measz", *grid, hw_master, time);
                }

                // Note we have to Hadamard-transform back after measurement in order to obtain correct expectation values
                else if (tile_spec == "double-horiz") { 
                    double time_tmp = lq2->transversal_op("measx", *grid, hw_master, time);
                    lq2->transversal_op("hadamard", *grid, hw_master, time_tmp);
                }

                // Measure strip qubits
                double time_tmp;
                for (unsigned int site : strip) {
                    time_tmp = time;
                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time, 0, *grid, hw_master);

                    time_tmp = TI_model.add_meas(site, time_tmp, 0, *grid, hw_master);

                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time_tmp, 0, *grid, hw_master);
                    
                }

                // Perform idle operation on lq1
                time = lq1->idle(cycles, *grid, hw_master, time);

            }

            // Merge two patches
            else if (s == "merge") {

                if (tile_spec == "single") {std::cerr << "merge: invalid tile_spec given." << std::endl; abort();}

                // Prepare strip qubits
                double time_tmp = 0;
                for (unsigned int site : strip) {
                    time_tmp = TI_model.add_init(site, time, 0, *grid, hw_master);

                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time_tmp, 1, *grid, hw_master);
                }

                // Perform 'idle' operation on the merged qubit
                time = lq->idle(cycles, *grid, hw_master, time);

            }

            // Split a two-tile patch
            else if (s == "split") {

                if (tile_spec == "single") {std::cerr << "split: invalid tile_spec given." << std::endl; abort();}

                // Measure strip qubits
                double time_tmp;
                for (unsigned int site : strip) {
                    time_tmp = time;
                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time, 0, *grid, hw_master);

                    time_tmp = TI_model.add_meas(site, time_tmp, 1, *grid, hw_master);

                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time_tmp, 2, *grid, hw_master);
                    
                }

                // Perform 'idle' operation on each separate qubit
                lq1->idle(cycles, *grid, hw_master, time);
                time = lq2->idle(cycles, *grid, hw_master, time);

            }

            else if (s == "bellmeas") {

                if (tile_spec == "single") {std::cerr << "bellmeas: invalid tile_spec given." << std::endl; abort();}

                // Prepare strip qubits
                double time_tmp = 0;
                for (unsigned int site : strip) {
                    time_tmp = TI_model.add_init(site, time, 0, *grid, hw_master);

                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time_tmp, 1, *grid, hw_master);
                }

                // Perform 'idle' operation on the merged qubit
                time = lq->idle(cycles, *grid, hw_master, time);

                // Measure out the whole merged patch
                if (tile_spec == "double-horiz") {
                    time = lq->transversal_op("measx", *grid, hw_master, time);
                }

                else if (tile_spec == "double-vert") {
                    time = lq->transversal_op("measz", *grid, hw_master, time);
                }

            }

            else if (s == "bellprep") {

                if (tile_spec == "single") {std::cerr << "bellprep: invalid tile_spec given." << std::endl; abort();}

                // Prepare state and do an idle on the merged patch
                else if (tile_spec == "double-horiz") {
                    lq->transversal_op("prepx", *grid, hw_master, time);
                }

                else if (tile_spec == "double-vert") {
                    lq->transversal_op("prepx", *grid, hw_master, time);
                }

                time = lq->idle(cycles, *grid, hw_master, time);

                // Measure strip qubits
                double time_tmp;
                for (unsigned int site : strip) {
                    time_tmp = time;
                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time, 0, *grid, hw_master);

                    time_tmp = TI_model.add_meas(site, time_tmp, 1, *grid, hw_master);

                    if (tile_spec == "double-horiz")
                        time_tmp = TI_model.add_H(site, time_tmp, 2, *grid, hw_master);
                    
                }

                // Perform 'idle' operation on each separate qubit
                lq1->idle(cycles, *grid, hw_master, time);
                time = lq2->idle(cycles, *grid, hw_master, time);

                // Post-processing: Pauli Z correction depending on measurement outcome (assumed to be tracked in software for now) (see notes)

            }

            else {std::cerr << "No valid operation selected. Options: {idle, prepz, prepx, measz, measx, hadamard, inject_y, inject_t, flip_patch, swap_left, rotation, extension, contraction, merge, split, bellprep, bellmeas}" << std::endl;}

            // Grab all of the occupied sites (to be used in printing)
            std::set<unsigned int> all_qsites = grid->get_occ_sites();

            // Enforce validity of final instruction list 
            grid->enforce_hw_master_validity(hw_master);

            // Print hardware instructions
            if (parser.exists("p")) {
                print_hw_master(std::cout, hw_master, all_qsites, debug);
            }

            // Count resources
            if (parser.exists("r")) {
                grid->resource_counter(hw_master);
            }

            // Free resources
            delete lq;
            if (tile_spec == "double-vert" || tile_spec == "double-horiz") {
                delete lq1;
                delete lq2;
            }
            delete grid;

        }

        return 0;
    }
}