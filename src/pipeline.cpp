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
                .description("Information to be queried (no operation will be performed). Options: {instructions, plaquettes, grid, parity}")
                .required(false);
        parser.add_argument()
                .names({"-o", "--operation"})
                .description("Surface code operation to be compiled. Options: {idle, prepz, prepx, measz, measx, extendx, extendz, mergex, mergez, splitx, splitz, bellprepx, bellprepz, bellmeasx, bellmeasz}")
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

        // Query information (no operation will be performed)
        if ((parser.exists("i")) && (parser.exists("o"))) {
            std::cerr << "The -i and -o flags are mutually exclusive." << std::endl;
        }
        else if (parser.exists("i"))
        {
            std::string s = parser.get<std::string>("i");
            
            /* TODO: This option does not depend on any grid allocation, so it should be decoupled from 
                command line options related to that. */
            if (s == "instructions") {
                HardwareModel TI_model;
                TI_model.print_TI_ops();
            }

            else if (s == "plaquettes") {
                GridManager grid(nrows, ncols);
                LogicalQubit lq(dx, dz, 0, 0, grid);
                lq.print_stabilizers();
            }
          
            else if (s == "grid") {
                GridManager grid(nrows, ncols);
                grid.print_qsite_mapping();
                std::cout << std::endl;
                std::vector<std::string> ascii_grid = grid.ascii_grid(false);
                grid.print_grid(ascii_grid);
                std::cout << std::endl;
                LogicalQubit lq(dx, dz, 0, 0, grid);
                ascii_grid = grid.ascii_grid(true);
                grid.print_grid(ascii_grid);
            }

            else if (s == "parity") {
                GridManager grid(nrows, ncols);
                LogicalQubit lq(dx, dz, 0, 0, grid);
                lq.print_parity_check_matrix();
                std::vector<HW_Instruction> hw_master;
            }

            else {
                std::cerr << "No valid query selected. Options: {instructions, plaquettes, grid, parity}" << std::endl;
            }
            
            return 0;
        }

        // Operation-dependent logic
        if (parser.exists("o"))
        {
            // Initialize vector of hardware instructions
            std::vector<HW_Instruction> hw_master;

            // Initialize hardware model
            HardwareModel TI_model;

            // Initialize time tracker
            double time = 0;

            std::string s = parser.get<std::string>("o");

            // Single-tile operations
            if ((s == "idle") || (s == "prepz") || (s == "prepx") || (s == "measz") || (s == "measx") || (s == "test")) {
                
                // Initializing grid using GridManager object. 
                GridManager grid(nrows, ncols);

                // Initialize logical qubit object using the grid
                LogicalQubit lq(dx, dz, 0, 0, grid);

                // Perform associated transversal operation
                if (s != "idle") {
                    lq.transversal_op(s, grid, hw_master, time);
                }

                // Append an idle operation if applicable 
                if ((s == "idle") || (s == "prepz") || (s == "prepx")) {
                    time = lq.idle(cycles, grid, hw_master, time);
                }

                // Grab all of the occupied sites (to be used in printing)
                std::set<unsigned int> occupied_sites = grid.get_occ_sites();

                // Enforce validity of final instruction list
                /* TODO: 
                    - Instead of 'enforcing' hardware validity in this fashion, maybe the circuits itself 
                    should contain Idles in locations where we know waiting might need to occur.
                    - However, what we have might make more sense since then a user does not need to account 
                    for validity up front, but rather the compiler 'makes it work' given hardware constraints. */ 
                grid.enforce_hw_master_validity(hw_master);

                // Print hardware instructions
                if (parser.exists("p")) {
                    print_hw_master(std::cout, hw_master, occupied_sites, debug);
                }

                // Count resources
                if (parser.exists("r")) {
                    grid.resource_counter(hw_master);
                }

            }

            // Horizontal two-tile operations 
            else if ((s == "contractx") || (s == "mergex") || (s == "bellmeasx") || (s == "extendx") ||
                (s == "splitx") || (s == "bellprepx") || (s == "hadamardx")) {
                
                // Construct grid with room for two tiles arranged horizontally
                GridManager grid(nrows, 2*ncols);

                // Initialize logical qubit object using the grid
                LogicalQubit lq1(dx, dz, 0, 0, grid);

                // Initialize second logical qubit object to the right of the first
                LogicalQubit lq2(dx, dz, 0, ncols, grid);

                // Create a merged qubit
                LogicalQubit* lq = merge(lq1, lq2, grid);  

                // Grab all of the qsites on the `strip' between lq1 and lq2
                std::set<unsigned int> strip = lq->get_strip(lq1, lq2);        

                // Debugging output
                if (debug) {
                    std::cout << "Logical Qubit 1:" << std::endl;
                    lq1.print_stabilizers();

                    std::cout << "Logical Qubit 2:" << std::endl;
                    lq2.print_stabilizers();

                    std::cout << "Logical Qubit (merged):" << std::endl;
                    lq->print_stabilizers();

                    std::cout << "Data qubits on strip:" << std::endl;
                    for (unsigned int site : strip) {
                        std::cout << site << std::endl;
                    }
                }

                // Operation-specific instructions
                if (s == "mergex") {

                    // Prepare qsites on the strip in the X basis
                    double time_tmp = 0;
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_init(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_H(site, time_tmp, 1, grid, hw_master);
                    }

                    // Perform 'idle' operation on the merged qubit
                    time = lq->idle(cycles, grid, hw_master, time);

                }

                else if (s == "extendx") {

                    // Prepare the physical qubits on lq2 in the X basis
                    lq2.transversal_op("prepx", grid, hw_master, time);

                    // Prepare qsites on the strip in the X basis
                    double time_tmp = 0;
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_init(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_H(site, time_tmp, 1, grid, hw_master);
                    }

                    // Perform 'idle' operation on the merged qubit
                    time = lq->idle(cycles, grid, hw_master, time);

                    // Pauli (X) correction depending on measurement outcome (X^m on final patch) (?) (see notes. i am not convinced this is needed as long as the correct mapping from two to one-qubit states is used.)

                }

                else if (s == "contractx") {

                    // Perform measure x on the half to be cropped
                    lq2.transversal_op("measx", grid, hw_master, time);

                    // Perform measure x on the strip
                    double time_tmp = 0;
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_H(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_meas(site, time_tmp, 1, grid, hw_master);
                    }

                    // Perform idle on remaining patch
                    time = lq1.idle(cycles, grid, hw_master, time);

                    // Pauli (Z) correction depending on measurement outcome (Z^m on final patch) (assumed to be tracked in software for now)

                }

                else if (s == "bellmeasx") {

                    // Prepare qsites on the strip in the X basis
                    double time_tmp = 0;
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_init(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_H(site, time_tmp, 1, grid, hw_master);
                    }

                    // Perform 'idle' operation on the merged qubit
                    time = lq->idle(cycles, grid, hw_master, time);

                    // Measure out the whole merged patch
                    time = lq->transversal_op("measx", grid, hw_master, time);

                }

                else if (s == "splitx") {

                    // Measure qsites on the strip in the X basis
                    double time_tmp = 0;
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_H(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_meas(site, time_tmp, 1, grid, hw_master);
                    }

                    // Perform 'idle' operation on each separate qubit
                    lq1.idle(cycles, grid, hw_master, time);
                    time = lq2.idle(cycles, grid, hw_master, time);

                }

                else if (s == "bellprepx") {

                    // Prepare X and do an idle on the merged patch
                    lq->transversal_op("prepx", grid, hw_master, time);
                    time = lq->idle(cycles, grid, hw_master, time);

                    // Measure qsites on the strip in the X basis
                    double time_tmp = 0;
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_H(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_meas(site, time_tmp, 1, grid, hw_master);
                    }

                    // Perform 'idle' operation on each separate qubit
                    lq1.idle(cycles, grid, hw_master, time);
                    time = lq2.idle(cycles, grid, hw_master, time);

                    // Pauli Z correction depending on measurement outcome (assumed to be tracked in software for now) (see notes)

                }

                else if (s == "hadamardx") {

                    // First apply transversal Hadamard to qubit 1
                    lq1.transversal_op("hadamard", grid, hw_master, time);

                    // Prepare the physical qubits on lq2 in the Z basis
                    lq2.transversal_op("prepz", grid, hw_master, time);

                    // Prepare qsites on the strip in the Z basis
                    for (unsigned int site : strip) {
                        TI_model.add_init(site, time, 0, grid, hw_master);
                    }

                    // Swap roles of X and Z for merged patch
                    lq->xz_swap(grid);

                    // Next extend patch rightward by running 'idle' on the merged patch
                    time = lq->idle(cycles, grid, hw_master, time);

                    // Corner movements: all measurements commute so I think they can be done at once

                    /* Note: extend_logical_operator_clockwise is still experimental and needs to be used with care. */
                    if (((lq->get_dx_init() % 2) == 0) || ((lq->get_dz_init() % 2) == 0)) {
                        std::cerr << "Patch rotation only currently implemented for odd code distances." << std::endl;
                        abort();
                    }
                    lq->extend_logical_operator_clockwise('X', "opposite", lq->get_dz_init()-1, grid, hw_master, time, debug); 
                    lq->extend_logical_operator_clockwise('Z', "default", lq->get_dz_init()-1, grid, hw_master, time, debug); 
                    lq->extend_logical_operator_clockwise('X', "default", lq->get_dx_init()-1, grid, hw_master, time, debug); 

                    std::cout << "Made it out." << std::endl;

                    // Reset stabilizer circuit patterns to default values. This should be improved later.
                    lq->reset_stabilizer_circuit_patterns();

                    // Measure the new stabilizers fault-tolerantly
                    time = lq->idle(cycles, grid, hw_master, time);

                    /* Contraction */
                    // Perform measure X on the half to be cropped
                    lq1.transversal_op("measx", grid, hw_master, time);

                    // Perform measure x on the strip
                    double time_tmp = 0;
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_H(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_meas(site, time_tmp, 1, grid, hw_master);
                    }

                    // Set time to start next operation
                    time = time_tmp;

                    /* Extension */
                    // Prepare the physical qubits on lq2 in the X basis
                    lq1.transversal_op("prepx", grid, hw_master, time);

                    // Prepare qsites on the strip in the X basis
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_init(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_H(site, time_tmp, 1, grid, hw_master);
                    }

                    time = lq->idle(cycles, grid, hw_master, time);

                    /* Final Contraction */
                    lq1.transversal_op("measx", grid, hw_master, time);

                    // Perform measure x on the strip
                    time_tmp = 0;
                    for (unsigned int site : strip) {
                        time_tmp = TI_model.add_H(site, time, 0, grid, hw_master);
                        time_tmp = TI_model.add_meas(site, time_tmp, 1, grid, hw_master);
                    }
                }

                // Grab all of the merged patch's occupied sites (to be used in printing)
                std::set<unsigned int> all_qsites = grid.get_occ_sites();

                // Enforce validity of final instruction list 
                grid.enforce_hw_master_validity(hw_master);

                // Print hardware instructions
                if (parser.exists("p")) {
                    print_hw_master(std::cout, hw_master, all_qsites, debug);
                }

                // Count resources
                if (parser.exists("r")) {
                    grid.resource_counter(hw_master);
                }

                // Free resources
                delete lq;
 
            }

            // Vertical two-patch operations
            else if ((s == "contractz") || (s == "mergez") || (s == "bellmeasz") || (s == "extendz") ||
                (s == "splitz") || (s == "bellprepz")) {

                // Construct grid with appropriate dimensions to hold two tiles top and bottom with a horizontal strip of qubits in between
                GridManager grid(2*nrows, ncols);

                // Initialize logical qubit object using the grid
                LogicalQubit lq1(dx, dz, 0, 0, grid);

                // Initialize second logical qubit object to the bottom of the first
                LogicalQubit lq2(dx, dz, nrows, 0, grid);

                // Create a merged qubit
                LogicalQubit* lq = merge(lq1, lq2, grid);

                // Grab all of the qsites on the `strip' between lq1 and lq2
                std::set<unsigned int> strip = lq->get_strip(lq1, lq2);

                // Debugging output
                if (debug) {
                    std::cout << "Logical Qubit 1:" << std::endl;
                    lq1.print_stabilizers();

                    std::cout << "Logical Qubit 2:" << std::endl;
                    lq2.print_stabilizers();

                    std::cout << "Logical Qubit (merged):" << std::endl;
                    lq->print_stabilizers();

                    std::cout << "Data qubits on strip:" << std::endl;
                    for (unsigned int site : strip) {
                        std::cout << site << std::endl;
                    }
                }

                // Operation-specific instructions
                if (s == "mergez") {

                    // Prepare qsites on the strip in the Z basis
                    for (unsigned int site : strip) {
                        TI_model.add_init(site, time, 0, grid, hw_master);
                    }

                    // Perform 'idle' operation on the merged qubit
                    time = lq->idle(cycles, grid, hw_master, time);

                }

                else if (s == "extendz") {

                    // Prepare the physical qubits on lq2 in the Z basis
                    lq2.transversal_op("prepz", grid, hw_master, time);

                    // Prepare qsites on the strip in the Z basis
                    for (unsigned int site : strip) {
                        TI_model.add_init(site, time, 0, grid, hw_master);
                    }

                    // Perform 'idle' operation on the merged qubit
                    time = lq->idle(cycles, grid, hw_master, time);

                    // Pauli (Z) correction depending on measurement outcome (Z^m on final patch) (?) (see notes. i am not convinced this is needed as long as the correct mapping from two to one-qubit states is used.)

                }

                else if (s == "contractz") {

                    // Perform measure z on the half to be cropped
                    lq2.transversal_op("measz", grid, hw_master, time);

                    // Perform measure z on the strip
                    for (unsigned int site : strip) {
                        TI_model.add_meas(site, time, 0, grid, hw_master);
                    }

                    // Perform idle on remaining patch
                    time = lq1.idle(cycles, grid, hw_master, time);

                    // Pauli (X) correction depending on measurement outcome (X^m) (assumed to be tracked in software for now)

                }

                else if (s == "bellmeasz") {

                    // Prepare qsites on the strip in the Z basis
                    for (unsigned int site : strip) {
                        TI_model.add_init(site, time, 0, grid, hw_master);
                    }

                    // Perform 'idle' operation on the merged qubit
                    time = lq->idle(cycles, grid, hw_master, time);

                    // Measure out the whole merged patch
                    time = lq->transversal_op("measz", grid, hw_master, time);

                }

                else if (s == "splitz") {

                    // Measure qsites on the strip in the Z basis
                    for (unsigned int site : strip) {
                        TI_model.add_meas(site, time, 0, grid, hw_master);
                    }

                    // Perform 'idle' operation on each separate qubit
                    lq1.idle(cycles, grid, hw_master, time);
                    time = lq2.idle(cycles, grid, hw_master, time);

                }

                else if (s == "bellprepz") {

                    // Prepare Z and idle on the merged patch
                    lq->transversal_op("prepz", grid, hw_master, time);
                    time = lq->idle(cycles, grid, hw_master, time);

                    // Measure qsites on the strip in the Z basis
                    for (unsigned int site : strip) {
                        TI_model.add_meas(site, time, 0, grid, hw_master);
                    }

                    // Perform 'idle' operation on each separate qubit
                    lq1.idle(cycles, grid, hw_master, time);
                    time = lq2.idle(cycles, grid, hw_master, time);

                    // Pauli X correction depending on measurement outcome (see notes) (assumed to be tracked in software for now)

                    /* TODO: consider whether this op can be done in a single `time step'. */

                }

                // Grab all of the larger patch's occupied sites (to be used in printing)
                std::set<unsigned int> all_qsites = grid.get_occ_sites();

                // Enforce validity of final instruction list 
                grid.enforce_hw_master_validity(hw_master);

                // Print hardware instructions
                if (parser.exists("p")) {
                    print_hw_master(std::cout, hw_master, all_qsites, debug);
                }

                // Count resources
                if (parser.exists("r")) {
                    grid.resource_counter(hw_master);
                }

                // Free resources
                delete lq;

            }

            else {std::cerr << "No valid operation selected. Options: {idle, prepz, prepx, measz, measx, extendx, extendz, mergex, mergez, splitx, splitz, bellprepx, bellprepz, bellmeasx, bellmeasz}" << std::endl;}
        }

        return 0;
    }
}