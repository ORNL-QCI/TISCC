#include <TISCC/logicalqubit.hpp>

#include <iomanip>
#include <cassert>
#include <algorithm>
#include <limits>

namespace TISCC 
{
    // Construct stabilizers and update set of occupied sites in grid
    void LogicalQubit::init_stabilizers(unsigned int dx, unsigned int dz, unsigned int row, unsigned int col, GridManager& grid) {
        
        // Construct the 4-qubit plaquettes of the surface code
        for (unsigned int i=row+1; i<row+dz; i++) {
            for (unsigned int j=col+1; j<col+dx; j++) {
                if (((i-row)+(j-col))%2 == 0) {
                    z_plaquettes.push_back(grid.get_plaquette(i, j, 'f', 'Z'));
                }
                else if (((i-row)+(j-col))%2 == 1) {
                    x_plaquettes.push_back(grid.get_plaquette(i, j, 'f', 'X'));
                }
            }
        }

        /* Construct the 2-qubit plaquettes of the surface code */
        int gap_tmp = 0;
        // Left edge of surface code
        for (unsigned int i=row+1; i<row+dz; i+=2) {
            x_plaquettes.push_back(grid.get_plaquette(i, col, 'w', 'X'));
            if (i==row+(dz-1)) {gap_tmp = 1;}
        }
        int gap = gap_tmp; gap_tmp = 0;
        for (unsigned int i=col+1; i<col+dx; i+=2) {
            // Top edge of surface code 
            if (i!=col+(dx-1)) {
                z_plaquettes.push_back(grid.get_plaquette(row, i+1, 'n', 'Z'));
                if (i+1==col+(dx-1)) {gap_tmp = 1;}
            }
            // Bottom edge of surface code
            if ((gap==1) && (i==col+(dx-1))) {continue;}
            z_plaquettes.push_back(grid.get_plaquette(row+dz, i+gap, 's', 'Z'));
        }
        gap = gap_tmp;
        // Right edge of surface code
        for (unsigned int i=row+1+gap; i<row+dz; i+=2) {
            x_plaquettes.push_back(grid.get_plaquette(i, col+dx, 'e', 'X'));
        }
    }

    // Set up the circuits that we intend to use
    void LogicalQubit::init_circuits() {

        // Add instructions to the Z plaquette measurement circuit
        Z_Circuit_Z_Type.push_back(Instruction("Prepare_Z", 'm', ' ')); 
        Z_Circuit_Z_Type.push_back(Instruction("Idle", ' ', ' '));   
        Z_Circuit_Z_Type.push_back(Instruction("CNOT", 'a', 'm'));
        Z_Circuit_Z_Type.push_back(Instruction("CNOT", 'b', 'm'));
        Z_Circuit_Z_Type.push_back(Instruction("CNOT", 'c', 'm'));
        Z_Circuit_Z_Type.push_back(Instruction("CNOT", 'd', 'm'));
        Z_Circuit_Z_Type.push_back(Instruction("Idle", ' ', ' ')); 
        Z_Circuit_Z_Type.push_back(Instruction("Measure_Z", 'm', ' '));

        // Add instructions to the X plaquette measurement circuit
        X_Circuit_N_Type.push_back(Instruction("Prepare_Z", 'm', ' '));  
        X_Circuit_N_Type.push_back(Instruction("Hadamard", 'm', ' '));   
        X_Circuit_N_Type.push_back(Instruction("CNOT", 'm', 'a'));
        X_Circuit_N_Type.push_back(Instruction("CNOT", 'm', 'c'));
        X_Circuit_N_Type.push_back(Instruction("CNOT", 'm', 'b'));
        X_Circuit_N_Type.push_back(Instruction("CNOT", 'm', 'd'));
        X_Circuit_N_Type.push_back(Instruction("Hadamard", 'm', ' '));   
        X_Circuit_N_Type.push_back(Instruction("Measure_Z", 'm', ' '));

        // Add instructions to the Z plaquette measurement circuit
        Z_Circuit_N_Type.push_back(Instruction("Prepare_Z", 'm', ' ')); 
        Z_Circuit_N_Type.push_back(Instruction("Idle", ' ', ' '));   
        Z_Circuit_N_Type.push_back(Instruction("CNOT", 'a', 'm'));
        Z_Circuit_N_Type.push_back(Instruction("CNOT", 'c', 'm'));
        Z_Circuit_N_Type.push_back(Instruction("CNOT", 'b', 'm'));
        Z_Circuit_N_Type.push_back(Instruction("CNOT", 'd', 'm'));
        Z_Circuit_N_Type.push_back(Instruction("Idle", ' ', ' ')); 
        Z_Circuit_N_Type.push_back(Instruction("Measure_Z", 'm', ' '));

        // Add instructions to the X plaquette measurement circuit
        X_Circuit_Z_Type.push_back(Instruction("Prepare_Z", 'm', ' '));  
        X_Circuit_Z_Type.push_back(Instruction("Hadamard", 'm', ' '));   
        X_Circuit_Z_Type.push_back(Instruction("CNOT", 'm', 'a'));
        X_Circuit_Z_Type.push_back(Instruction("CNOT", 'm', 'b'));
        X_Circuit_Z_Type.push_back(Instruction("CNOT", 'm', 'c'));
        X_Circuit_Z_Type.push_back(Instruction("CNOT", 'm', 'd'));
        X_Circuit_Z_Type.push_back(Instruction("Hadamard", 'm', ' '));   
        X_Circuit_Z_Type.push_back(Instruction("Measure_Z", 'm', ' '));
    }

    void LogicalQubit::print_stabilizers() {

        // Print all z stabilizers
        for (const Plaquette& p : z_plaquettes) {
            std::cout << p.get_type() << " " << p.get_row() << " " << p.get_col() << " " << p.get_shape()
                 << " " << (p.grid()).index_from_coords(p.get_row(), p.get_col(), 1) << std::endl;
        }

        // Print all x stabilizers
        for (const Plaquette& p : x_plaquettes) {
            std::cout << p.get_type() << " " << p.get_row() << " " << p.get_col() << " " << p.get_shape()
                 << " " << (p.grid()).index_from_coords(p.get_row(), p.get_col(), 1) << std::endl;
        }
    }

    // Construct parity check matrix. Rows are plaquettes and columns refer to qsites (repeated twice)
    void LogicalQubit::construct_parity_check_matrix(const GridManager& grid) {

        // Construct qsites to indices map
        std::map<unsigned int, unsigned int> qsite_to_index_;
        std::set<unsigned int> data = data_qsites();
        unsigned int counter = 0;
        for (unsigned int site: data) {
            qsite_to_index_[site] = counter;
            counter++;
        }
        qsite_to_index = std::move(qsite_to_index_);

        // Construct parity check matrix to represent both stabilizers and observables
        std::vector<std::vector<bool> > parity_check_(x_plaquettes.size() + z_plaquettes.size() + 2, std::vector<bool>(2*data.size(), 0));

        // Iterate over z stabilizers
        counter = 0;
        for (const Plaquette& p : z_plaquettes) {
            for (char qubit : {'a', 'b', 'c', 'd'}) {
                if (qsite_to_index->count(p.get_qsite(qubit))) 
                    parity_check_[counter][qsite_to_index.value()[p.get_qsite(qubit)]] = 1;
            }
            counter++;
        }

       // Iterate over x stabilizers
        for (const Plaquette& p : x_plaquettes) {
            for (char qubit : {'a', 'b', 'c', 'd'}) {
                if (qsite_to_index->count(p.get_qsite(qubit))) 
                    parity_check_[counter][data.size() + qsite_to_index.value()[p.get_qsite(qubit)]] = 1;
            }
            counter++;
        }

        // Obtain observables
        for (unsigned int qsite : data) {
            // The col_'th column of the grid contains the left-most column of data qubits contained in the surface code patch
            if (grid.get_col(qsite) == col_) {
                parity_check_[counter][qsite_to_index.value()[qsite]] = 1;
            }
            // The (row_+1)'th row of the grid contains the top-most row of data qubits contained in the surface code patch
            if (grid.get_row(qsite) == row_+1) {
                parity_check_[counter+1][data.size() + qsite_to_index.value()[qsite]] = 1;
            }
        }

        parity_check_matrix = std::move(parity_check_);
    }

    // Construct and print parity check matrix. rows are plaquettes and columns refer to qsites (repeated twice)
    void LogicalQubit::print_parity_check_matrix(const GridManager& grid) {

        if (!parity_check_matrix) {
            construct_parity_check_matrix(grid);
        }

        // Print map
        std::cout << "Map from column indices to qsites:" << std::endl;
        for (std::pair<unsigned int, unsigned int> pair : qsite_to_index.value()) {
            std::cout << pair.second << " " << pair.first << std::endl;
        }
        for (std::pair<unsigned int, unsigned int> pair : qsite_to_index.value()) {
            std::cout << qsite_to_index->size() + pair.second << " " << pair.first << std::endl;
        }

        // Print matrix
        std::cout << std::endl;
        std::cout << "Parity check matrix (rows in same order as given by print_stabilizers function):" << std::endl;
        for (unsigned int i = 0; i<parity_check_matrix->size(); i++) {
            for (unsigned int j = 0; j<parity_check_matrix.value()[i].size(); j++) {
                std::cout << parity_check_matrix.value()[i][j];
            }
            std::cout << std::endl;
        }
    }

    // Test stabilizers (not fully implemented)
    void LogicalQubit::test_stabilizers(unsigned int dx, unsigned int dz) {
        unsigned int num_data_qubits = dx*dz;
        unsigned int num_measure_qubits = num_data_qubits-1;

        // Check that number of stabilizers is equal to the number expected
        unsigned int num_plaquettes = z_plaquettes.size() + x_plaquettes.size();
        if (num_plaquettes != num_measure_qubits) {
            std::cerr << "Incorrect number of plaquettes: " << num_plaquettes << " " << num_measure_qubits << std::endl; 
            abort();
        }

        // Check that measure qubits are all unique
        std::vector<unsigned int> measure_qubits;
        for (Plaquette p : z_plaquettes) {
            for (unsigned int q : measure_qubits) {
                if (p.get_qsite('m') == q) {std::cerr << "Found duplicate measure qubit among plaquettes."; abort();}
            }
            measure_qubits.push_back(p.get_qsite('m'));
        }
        for (Plaquette p : x_plaquettes) {
            for (unsigned int q : measure_qubits) {
                if (p.get_qsite('m') == q) {std::cerr << "Found duplicate measure qubit among plaquettes."; abort();}
            }
            measure_qubits.push_back(p.get_qsite('m'));
        }

        // TODO: Create a way to find all plaquettes containing a particular qubit
    }

    LogicalQubit::LogicalQubit(unsigned int dx, unsigned int dz, unsigned int row, unsigned int col, GridManager& grid) : 
        dx_(dx), dz_(dz), row_(row), col_(col), TI_model() { 
        init_stabilizers(dx, dz, row, col, grid);
        init_circuits();
        test_stabilizers(dx, dz);
    }

    // Function to return all qsites occupied by the surface code
    std::set<unsigned int> LogicalQubit::occupied_sites() {
        std::set<unsigned int> sites;
        for (char qubit : {'a', 'b', 'c', 'd', 'm'}) {
            for (const Plaquette& p : z_plaquettes) {
                sites.insert(p.get_qsite(qubit));
            }
            for (const Plaquette& p : x_plaquettes) {
                sites.insert(p.get_qsite(qubit));
            }
        }
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        sites.erase(uint_max);
        return sites;
    }

    // Function to return all qsites occupied by data qubits on the surface code
    /* TODO: Combine this somehow with occupied_sites by passing in a string with default value abcdm. */
    std::set<unsigned int> LogicalQubit::data_qsites() {
        std::set<unsigned int> sites;
        for (char qubit : {'a', 'b', 'c', 'd'}) {
            for (const Plaquette& p : z_plaquettes) {
                sites.insert(p.get_qsite(qubit));
            }
            for (const Plaquette& p : x_plaquettes) {
                sites.insert(p.get_qsite(qubit));
            }
        }
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        sites.erase(uint_max);
        return sites;
    }

    // Apply a given ``qubit-level'' instruction to all plaquettes in a given vector and add corresponding HW_Instructions to hw_master
    double LogicalQubit::apply_instruction(const Instruction& instr, std::vector<Plaquette>& plaquettes, double time, unsigned int step, 
        const GridManager& grid, std::vector<HW_Instruction>& hw_master) {

        // Create tmp variable for time
        double time_tmp = 0;

        // Don't explicitly apply an "Idle" operation
        if (instr.get_name() != "Idle") {

            // Loop over plaquettes
            for (Plaquette& p : plaquettes) {

                // Certain stabilizer shapes do not accomodate certain qubit labels
                if (p.is_instr_valid(instr)) 
                {
                    
                    // Add single-qubit instructions
                    if ((instr.get_name() == "Prepare_Z") && (instr.get_q2() == ' ')) {
                        time_tmp = TI_model.add_init(p, instr.get_q1(), time, step, hw_master);
                    }

                    else if ((instr.get_name() == "Hadamard") && (instr.get_q2() == ' ')) {
                        time_tmp = TI_model.add_H(p, instr.get_q1(), time, step, hw_master);
                    }

                    else if ((instr.get_name() == "Measure_Z") && (instr.get_q2() == ' ')) {
                        time_tmp = TI_model.add_meas(p, instr.get_q1(), time, step, hw_master);
                    }

                    else if ((instr.get_name() == "Test_Gate") && (instr.get_q2() == ' ')) {
                        time_tmp = TI_model.add_test(p, instr.get_q1(), time, step, hw_master);
                    }

                    // Add CNOT gate
                    else if (instr.get_name() == "CNOT") {
                        time_tmp = TI_model.add_CNOT(p, instr.get_q1(), instr.get_q2(), time, step, grid, hw_master);
                    }

                    else {std::cerr << "LogicalQubit::apply_instruction: Invalid instruction given." << std::endl; abort();}
                    
                }
            }
        }

        return time_tmp;
    }

    double LogicalQubit::idle(unsigned int cycles, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {
        
        // Loop over surface code cycles
        for (unsigned int cycle=0; cycle < cycles; cycle++) {

            // Enforce that the two circuits contain the same number of instructions
            assert(Z_Circuit_Z_Type.size() == X_Circuit_N_Type.size());
            unsigned int num_instructions = Z_Circuit_Z_Type.size();

            // Loop over ``qubit-level'' instructions and apply them to plaquettes while adding HW_Instructions to hw_master
            for (unsigned int i=0; i<num_instructions; i++) {
                double t1 = apply_instruction(Z_Circuit_Z_Type[i], z_plaquettes, time, i, grid, hw_master);
                double t2 = apply_instruction(X_Circuit_N_Type[i], x_plaquettes, time, i, grid, hw_master);

                // Increment time counter
                if (t1 == 0) {time = t2;}
                else if (t2 == 0) {time = t1;}
                else {assert(t1==t2); time = t1;}  

            }
            
        }

        // Sort master list of hardware instructions according to overloaded operator<
        std::stable_sort(hw_master.begin(), hw_master.end());

        return time;

    }

    // Swap roles of x and z for this patch (used during Hadamard and patch rotation)
    void LogicalQubit::xz_swap() {

        // Every X operator becomes a Z operator and vice versa

        // Thus, Z are horizontal and X are vertical

        // Thus, Z stabilizer circuits need to be N type and X stabilizer circuits need to be Z type

        // x_plaquettes and z_plaquettes vectors should be swapped and dx and dz should be swapped

        // Z_Circuit_N_Type and X_Circuit_Z_Type should be created and used for any subsequent operation

        // A patch rotation needs to be required after this, because other operations (such as merge) will no longer be compatible with this patch

    }

    double LogicalQubit::transversal_op(const std::string& op, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {
        
        // Get all occupied sites
        std::set<unsigned int> sites = data_qsites();

        // Create tmp variable for time
        double time_tmp = 0;

        // TODO: Check whether all data qubits are at their 'home' positions

        // Loop over all occupied sites
        for (unsigned int site : sites) {

            // Add instructions by case
            if (op == "prepz") {
                time_tmp = TI_model.add_init(site, time, 0, grid, hw_master);
            }

            else if (op == "prepx") {
                time_tmp = TI_model.add_init(site, time, 0, grid, hw_master);
                time_tmp = TI_model.add_H(site, time_tmp, 1, grid, hw_master);
            }

            else if (op == "measz") {
                time_tmp = TI_model.add_meas(site, time, 0, grid, hw_master);
            }

            else if (op == "measx") {
                time_tmp = TI_model.add_H(site, time, 0, grid, hw_master);
                time_tmp = TI_model.add_meas(site, time_tmp, 1, grid, hw_master);
            }

            else if (op == "hadamard") {
                time_tmp = TI_model.add_H(site, time, 0, grid, hw_master);
            }
        }

        // Sort master list of hardware instructions according to overloaded operator<
        std::stable_sort(hw_master.begin(), hw_master.end());

        return time_tmp;

    }

    // Helper function to return the data qubits from this patch that are NOT occupied by two others
    std::set<unsigned int> LogicalQubit::get_strip(LogicalQubit& lq1, LogicalQubit& lq2) {
        std::set<unsigned int> data = data_qsites();
        std::set<unsigned int> lq1_data = lq1.data_qsites();
        std::set<unsigned int> lq2_data = lq2.data_qsites();
        std::set<unsigned int> strip;
        for (unsigned int lq_site : data) {
            if ((lq1_data.find(lq_site) == lq1_data.end()) && (lq2_data.find(lq_site) == lq2_data.end())) {
                strip.insert(lq_site);
            }
        }
        return strip;
    }

    // Placeholder function to help implement little test circuits
    double LogicalQubit::test_circuits(const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {
        // Bell state preparation
        apply_instruction(Instruction("Prepare_Z", 'a', ' '), z_plaquettes, time, 0, grid, hw_master);
        time = apply_instruction(Instruction("Prepare_Z", 'm', ' '), z_plaquettes, time, 0, grid, hw_master);
        time = apply_instruction(Instruction("Hadamard", 'a', ' '), z_plaquettes, time, 1, grid, hw_master);
        time = apply_instruction(Instruction("CNOT", 'a', 'm'), z_plaquettes, time, 2, grid, hw_master);
        apply_instruction(Instruction("Measure_Z", 'a', ' '), z_plaquettes, time, 3, grid, hw_master);
        time = apply_instruction(Instruction("Measure_Z", 'm', ' '), z_plaquettes, time, 3, grid, hw_master);

        // SWAP gate
        // apply_instruction(Instruction("Prepare_Z", 'a', ' '), z_plaquettes, time, 0, grid, hw_master);
        // time = apply_instruction(Instruction("Prepare_Z", 'm', ' '), z_plaquettes, time, 0, grid, hw_master);
        // time = apply_instruction(Instruction("Hadamard", 'a', ' '), z_plaquettes, time, 1, grid, hw_master);
        // time = apply_instruction(Instruction("CNOT", 'a', 'm'), z_plaquettes, time, 2, grid, hw_master);
        // time = apply_instruction(Instruction("CNOT", 'm', 'a'), z_plaquettes, time, 3, grid, hw_master);
        // time = apply_instruction(Instruction("CNOT", 'a', 'm'), z_plaquettes, time, 4, grid, hw_master);
        // time = apply_instruction(Instruction("Hadamard", 'm', ' '), z_plaquettes, time, 5, grid, hw_master);
        // time = apply_instruction(Instruction("Measure_Z", 'a', ' '), z_plaquettes, time, 6, grid, hw_master);
        // time = apply_instruction(Instruction("Measure_Z", 'm', ' '), z_plaquettes, time, 6, grid, hw_master);

        // time = apply_instruction(Instruction("Prepare_Z", 'a', ' '), z_plaquettes, time, 0, grid, hw_master);
        // time = apply_instruction(Instruction("Hadamard", 'a', ' '), z_plaquettes, time, 1, grid, hw_master);
        // time = apply_instruction(Instruction("Test_Gate", 'a', ' '), z_plaquettes, time, 2, grid, hw_master);
        // time = apply_instruction(Instruction("Hadamard", 'a', ' '), z_plaquettes, time, 3, grid, hw_master);
        // time = apply_instruction(Instruction("Measure_Z", 'a', ' '), z_plaquettes, time, 4, grid, hw_master);

        return time;
    }

// Construct and return a logical qubit that represents the merged product of two input logical qubits
LogicalQubit merge(LogicalQubit& lq1, LogicalQubit& lq2, GridManager& grid) {

    // Determine whether to merge horizontally or vertically and set parameters

    // If they are horizontally displaced,
    if (lq1.get_row() == lq2.get_row()) {

        // We require them to have the same code distance
        assert(lq1.get_dz() == lq2.get_dz());

        // Merge horizontally
        unsigned int extra_strip = 0;
        if (lq1.get_dx()%2 == 0) {extra_strip = 1;}
        assert(lq2.get_col() == lq1.get_col() + lq1.get_dx() + 1 + extra_strip);
        return LogicalQubit(lq1.get_dx() + lq2.get_dx() + 1 + extra_strip, lq1.get_dz(), lq1.get_row(), lq1.get_col(), grid);

    }

    // Otherwise, if they are vertically displaced,
    else if (lq1.get_col() == lq2.get_col()) {

        // We require them to have the same x code distance
        assert(lq1.get_dx() == lq2.get_dx());

        // Merge vertically
        unsigned int extra_strip = 0;
        if (lq1.get_dz()%2 == 0) {extra_strip = 1;} 
        assert(lq2.get_row() == lq1.get_row() + lq1.get_dz() + 1 + extra_strip);
        return LogicalQubit(lq1.get_dx(), lq1.get_dz() + lq2.get_dz() + 1 + extra_strip, lq1.get_row(), lq1.get_col(), grid); 

    }

    else {
        std::cerr << "merge: this operation must take place between logical qubits either vertically or horizontally separated, but not both." << std::endl;
        abort();
    }



}
}