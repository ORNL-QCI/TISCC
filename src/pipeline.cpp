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
                .description("Surface code operation to be compiled. Options: {idle, prepz, prepx, measz, measx, pauli_x, pauli_y, pauli_z, hadamard, inject_y, inject_t, flip_patch, move_right, swap_left, extension, contraction, move, merge, split, jointmeas, mergecontract, extendsplit, bellprep, bellmeas}")
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
        LogicalQubit* lq = nullptr;
        LogicalQubit* lq1 = nullptr;
        LogicalQubit* lq2 = nullptr;
        GridManager* grid = nullptr;

        // Extract tile_spec
        std::string tile_spec = "single";
        if (parser.exists("s")) {
            tile_spec = parser.get<std::string>("s");
        }

        // Some operations require the usage of an additional column to the right
        bool extra_col = 0;
        if (parser.exists("o")) {
            std::string s = parser.get<std::string>("o");
            if ((s == "move_right") || (s == "swap_left") || (s == "single_tile_rotation") || (s == "hadamard")) {
                extra_col = 1;
            }
        }

        // Case-dependent grid & lq initialization
        if (tile_spec == "single") {

            // Initialize grid
            grid = new GridManager(nrows, ncols + extra_col);

            // Initialize logical qubit on  grid
            lq = new LogicalQubit(dx, dz, 0, 0, *grid);
        }

        else if (tile_spec == "double-vert") {

            // Construct grid with appropriate dimensions to hold two tiles top and bottom with a horizontal strip of qubits in between
            grid = new GridManager(2*nrows, ncols + extra_col);

            // Initialize logical qubit object using the grid
            lq1 = new LogicalQubit(dx, dz, 0, 0, *grid);

            // Initialize second logical qubit object to the bottom of the first
            lq2 = new LogicalQubit(dx, dz, nrows, 0, *grid);

            // Create a merged qubit
            lq = lq1->get_merged_lq(*lq2, *grid);
        }

        else if (tile_spec == "double-horiz") {

            // Construct grid with appropriate dimensions to hold two tiles top and bottom with a horizontal strip of qubits in between
            grid = new GridManager(nrows, 2*ncols + extra_col);

            // Initialize logical qubit object using the grid
            lq1 = new LogicalQubit(dx, dz, 0, 0, *grid);

            // Initialize second logical qubit object to the right of the first
            lq2 = new LogicalQubit(dx, dz, 0, ncols, *grid);

            // Create a merged qubit
            lq = lq1->get_merged_lq(*lq2, *grid);

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
                std::cout << std::endl;
                ascii_grid = grid->ascii_grid_with_operator(lq->syndrome_measurement_qsites(), true);
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
            if ((s == "idle") || (s == "prepz") || (s == "prepx") || (s == "measz") || (s == "measx") || (s == "pauli_x") || (s == "pauli_y") || (s == "pauli_z") 
            || (s == "inject_y") || (s == "inject_t") || (s == "flip_patch")
            || (s == "hadamard") || (s == "move_right") || (s == "swap_left")) {

                // Perform associated transversal operation
                if ((s == "prepz") || (s == "prepx") || (s == "measz") || (s == "measx") || (s == "hadamard")) {
                    lq->transversal_op(s, *grid, hw_master, time);
                }

                else if ((s == "inject_y") || (s == "inject_t")) {
                    char state_label = s.substr(7)[0];
                    lq->inject_state(state_label, *grid, hw_master, time);
                }

                else if ((s == "pauli_x") || (s == "pauli_y") || (s == "pauli_z")) {
                    char pauli_label = s.substr(6)[0];
                    lq->apply_pauli(pauli_label, *grid, hw_master, time);
                }

                else if (s == "flip_patch") {
                    time = lq->flip_patch(*grid, hw_master, time, true, debug);
                }

                else if (s == "move_right") {

                    // move_right takes in pointers to lq_extended and lq_contracted that can be used later if desired
                    LogicalQubit* lq_extended = nullptr;
                    LogicalQubit* lq_contracted = nullptr;
                    time = lq->move_right(cycles, lq_extended, lq_contracted, *grid, hw_master, time);

                    // Visualize if debug flag is on
                    std::vector<std::string> ascii_grid;
                    if (debug) {
                        std::cout << "Configuration before move_right:" << std::endl;
                        ascii_grid = grid->ascii_grid_with_operator(lq->syndrome_measurement_qsites(), true);
                        grid->print_grid(ascii_grid);
                        std::cout << "Configuration after extension:" << std::endl;
                        ascii_grid = grid->ascii_grid_with_operator(lq_extended->syndrome_measurement_qsites(), true);
                        grid->print_grid(ascii_grid);
                        std::cout << "Configuration after move_right:" << std::endl;
                        ascii_grid = grid->ascii_grid_with_operator(lq_contracted->syndrome_measurement_qsites(), true);
                        grid->print_grid(ascii_grid);
                    }

                    // Transfer resources appropriately for later processing
                    if (lq1 != nullptr) {
                        delete lq1;
                        lq1 = nullptr;
                    }
                    if (lq2 != nullptr) {
                        delete lq2;
                        lq2 = nullptr;
                    }

                    delete lq;
                    lq = nullptr;
                    lq = lq_contracted; // Contracted patch, being final state, defines logical operators for later processing
                    lq2 = lq_extended; // Retain extended patch for calculating operator deformations

                }

                else if (s == "swap_left") {
                    
                    // re-allocate lq and swap x<->z stabilizers
                    LogicalQubit* tmp_lq = lq;
                    lq = new LogicalQubit(tmp_lq->get_dx_init(), tmp_lq->get_dz_init(), 0, 1, *grid);
                    delete tmp_lq;

                    std::vector<std::string> ascii_grid;
                    if (debug) {
                        std::cout << "Configuration before swap_left:" << std::endl;
                        ascii_grid = grid->ascii_grid_with_operator(lq->syndrome_measurement_qsites(), true);
                        grid->print_grid(ascii_grid);
                    }
                    
                    time = lq->swap_left(*grid, hw_master, time);

                    if (debug) {
                        std::cout << std::endl << "Configuration after swap_left:" << std::endl;
                        ascii_grid = grid->ascii_grid_with_operator(lq->syndrome_measurement_qsites(), true);
                        grid->print_grid(ascii_grid);
                    }

                }

                // Append an idle operation if applicable 
                if ((s == "idle") || (s == "prepz") || (s == "prepx") || (s == "inject_y") || (s == "inject_t") || (s == "flip_patch")) {
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

                // Perform merge operation
                time = lq->merge(cycles, *grid, hw_master, time);

            }

            // Patch contractions
            else if (s == "contraction") {

                if (tile_spec == "single") {std::cerr << "contraction: invalid tile_spec given." << std::endl; abort();}

                // Measure lq2 in approp. basis depending on direction of contraction
                else if (tile_spec == "double-vert") {
                    lq2->transversal_op("measz", *grid, hw_master, time);
                }

                else if (tile_spec == "double-horiz") { 
                    lq2->transversal_op("measx", *grid, hw_master, time);
                }

                // Perform split operation
                time = lq->split(*grid, hw_master, time);
            }

            // Move (extension followed by contraction)
            else if (s == "move") {

                if (tile_spec == "single") {std::cerr << "move: invalid tile_spec given." << std::endl; abort();}

                // Prepare lq2 in approp. basis depending on direction of extension
                else if (tile_spec == "double-horiz") {
                    lq2->transversal_op("prepx", *grid, hw_master, time);
                }

                else if (tile_spec == "double-vert") {
                    lq2->transversal_op("prepz", *grid, hw_master, time);
                }

                // Perform merge operation
                time = lq->merge(cycles, *grid, hw_master, time);

                // Measure lq1 in approp. basis depending on direction of contraction
                if (tile_spec == "double-vert") {
                    lq1->transversal_op("measz", *grid, hw_master, time);
                }

                else if (tile_spec == "double-horiz") { 
                    lq1->transversal_op("measx", *grid, hw_master, time);
                }

                // Perform split operation
                time = lq->split(*grid, hw_master, time);

            }

            // Merge two patches
            else if (s == "merge") {

                if (tile_spec == "single") {std::cerr << "merge: invalid tile_spec given." << std::endl; abort();}

                // Perform merge operation
                time = lq->merge(cycles, *grid, hw_master, time);

            }

            // Split a two-tile patch
            else if (s == "split") {

                if (tile_spec == "single") {std::cerr << "split: invalid tile_spec given." << std::endl; abort();}

                // Perform split operation
                time = lq->split(*grid, hw_master, time);

            }

            // Combine into joint XX or ZZ measurement
            else if (s == "jointmeas") {

                if (tile_spec == "single") {std::cerr << "jointmeas: invalid tile_spec given." << std::endl; abort();}

                // Perform merge operation
                time = lq->merge(cycles, *grid, hw_master, time);

                // Perform split operation
                time = lq->split(*grid, hw_master, time);

            }

            // Merge-Contract
            else if (s == "mergecontract") {

                if (tile_spec == "single") {std::cerr << "mergecontract: invalid tile_spec given." << std::endl; abort();}

                // Perform merge operation
                time = lq->merge(cycles, *grid, hw_master, time);

                // Measure lq2 in approp. basis depending on direction of contraction
                if (tile_spec == "double-vert") {
                    lq2->transversal_op("measz", *grid, hw_master, time);
                }

                else if (tile_spec == "double-horiz") { 
                    lq2->transversal_op("measx", *grid, hw_master, time);
                }

                // Perform split operation
                time = lq->split(*grid, hw_master, time);

            }

            // Extend-Split
            else if (s == "extendsplit") {

                if (tile_spec == "single") {std::cerr << "extension: invalid tile_spec given." << std::endl; abort();}

                // Prepare lq2 in approp. basis depending on direction of extension
                else if (tile_spec == "double-horiz") {
                    lq2->transversal_op("prepx", *grid, hw_master, time);
                }

                else if (tile_spec == "double-vert") {
                    lq2->transversal_op("prepz", *grid, hw_master, time);
                }

                // Perform merge operation
                time = lq->merge(cycles, *grid, hw_master, time);

                // Perform split operation
                time = lq->split(*grid, hw_master, time);

            }

            else if (s == "bellmeas") {

                if (tile_spec == "single") {std::cerr << "bellmeas: invalid tile_spec given." << std::endl; abort();}

                // Perform merge operation
                time = lq->merge(cycles, *grid, hw_master, time);

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

                // Prepare lq2 in approp. basis depending on direction of extension
                else if (tile_spec == "double-horiz") {
                    lq1->transversal_op("prepx", *grid, hw_master, time);
                    lq2->transversal_op("prepx", *grid, hw_master, time);
                }

                else if (tile_spec == "double-vert") {
                    lq1->transversal_op("prepz", *grid, hw_master, time);
                    lq2->transversal_op("prepz", *grid, hw_master, time);
                }

                // Perform merge operation
                time = lq->merge(cycles, *grid, hw_master, time);

                // Perform split operation
                time = lq->split(*grid, hw_master, time);

                // Post-processing: Pauli Z correction depending on measurement outcome

            }

            else {std::cerr << "No valid operation selected. Options: {idle, prepz, prepx, measz, measx, pauli_x, pauli_y, pauli_z, hadamard, inject_y, inject_t, flip_patch, move_right, swap_left, extension, contraction, move, merge, split, jointmeas, mergecontract, extendsplit, bellprep, bellmeas}" << std::endl;}

            // Grab all of the occupied sites (to be used in printing)
            // **Note: This being after circuit generation assumes that the occupancy of the grid post-circuit is equivalent to the occupancy pre-circuit
            std::set<unsigned int> all_qsites = grid->get_occ_sites();

            // Enforce validity of final instruction list 
            std::stable_sort(hw_master.begin(), hw_master.end());
            grid->enforce_hw_master_validity(hw_master);

            // Print hardware instructions
            if (parser.exists("p")) {
                HW_Instruction::print_hw_master(std::cout, hw_master, all_qsites, debug);
            }

            // Count resources
            if (parser.exists("r")) {
                grid->resource_counter(hw_master);
            }

            // Free resources
            delete lq;
            if (lq1 != nullptr) delete lq1;
            if (lq2 != nullptr) delete lq2;
            delete grid;

        }

        return 0;
    }
}