#include <TISCC/logicalqubit.hpp>

#include <iomanip>
#include <cassert>
#include <algorithm>
#include <limits>
#include <iterator>
#include <utility>

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

    // Transform operators from binary representation to pair<qsite unsigned int, Pauli char>
    std::vector<std::pair<unsigned int, char>> LogicalQubit::binary_operator_to_qsites(const std::vector<bool>& binary_rep) {
        std::vector<std::pair<unsigned int, char>> qsite_rep;
        for (unsigned int k = 0; k < 2*qsite_to_index->size(); k++) {
            if (binary_rep[k] && (k < qsite_to_index->size())) {
                qsite_rep.emplace_back(index_to_qsite.value()[k], 'Z');
            }
            else if (binary_rep[k]) {
                qsite_rep.emplace_back(index_to_qsite.value()[k], 'X');
            }
        }
        return qsite_rep;
    }

    // Add new stabilizer plaquette (to be used in corner movement)
    void LogicalQubit::add_stabilizer(unsigned int row, unsigned int col, char shape, char type, GridManager& grid, bool debug) {

        // We require newly measured stabilizers to be along the boundaries and of the appropriate shape, since the purpose of this function is to be used in corner movement
        if (!(((col == col_) && (row > row_) && (row < row_ + dz_) && (shape == 'w')) ||  // Left edge
            ((col == col_ + dx_) && (row > row_) && (row < row_ + dz_) && (shape == 'e')) || // Right edge
            ((row == row_) && (col > col_) && (col < col_ + dx_) && (shape == 'n')) ||  // Top edge
            ((row == row_ + dz_) && (col > col_) && (col < col_ + dx_) && (shape == 's')))) { // Bottom edge
            std::cout << "LogicalQubit::add_stabilizer: Newly added plaquettes must lie along boundaries and be of the appropriate shape, as this function is only meant to be used in corner movement." << std::endl;
            abort();
        }

        // Construct parity check matrix if haven't already
        if (!parity_check_matrix) {
            construct_parity_check_matrix(grid);
        }

        // Construct new Plaquette from the grid 
        // (note validity checks take place within this function, and it also modifies the set of occupied_sites in grid)
        Plaquette new_stabilizer(grid.get_plaquette(row, col, shape, type)); 

        // Construct a new vector in binary symplectic format corresponding to this stabilizer
        std::vector<bool> new_row(2*qsite_to_index->size(), 0);
        for (char qubit : {'a', 'b', 'c', 'd'}) {
            if (qsite_to_index->count(new_stabilizer.get_qsite(qubit))) {
                if (type == 'Z') {
                    new_row[qsite_to_index.value()[new_stabilizer.get_qsite(qubit)]] = 1;
                }
                else if (type == 'X') {
                    new_row[qsite_to_index->size() + qsite_to_index.value()[new_stabilizer.get_qsite(qubit)]] = 1;
                }
            }
        }

        /* Check commutation of new stabilizer with every row of parity_check_matrix (with appended logical operators) */

        // First, transform the new row s.t. the roles of X and Z have been swapped
        std::vector<bool> new_row_transformed(2*qsite_to_index->size());
        for (unsigned int k=0; k<qsite_to_index->size(); k++) {
            new_row_transformed[k] = new_row[k+qsite_to_index->size()];
            new_row_transformed[k+qsite_to_index->size()] = new_row[k];
        }

        // Print row and transformed row (for debugging purposes)
        if (debug) {
            std::cout << std::endl << "Stabilizer being added: ";
            std::copy(new_row.begin(), new_row.end(), std::ostream_iterator<bool>(std::cout));
            std::cout << std::endl;
            std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(new_row), true);
            grid.print_grid(ascii_grid);
        }


        // Then, calculate the binary symplectic product with every row of the parity_check_matrix and track row indices for which it is 1
        std::vector<unsigned int> anticommuting_stabilizers;
        std::vector<unsigned int> anticommuting_logical_ops;
        bool bin_sym_prod;
        for (unsigned int i=0; i<parity_check_matrix->size(); i++) {

            // Calculate binary symplectic prodcut
            bin_sym_prod = 0;
            for (unsigned int k = 0; k < new_row_transformed.size(); k++) {
                bin_sym_prod ^= (parity_check_matrix.value()[i][k] && new_row_transformed[k]);
            }

            // If the product = 1, then the two operators anti-commute
            if (bin_sym_prod) {

                // If this row of parity_check_matrix corresponds with a stabilizer plaquette,
                if (i < parity_check_matrix->size() - 2) {

                    // We don't allow the new stabilizer to anti-commute with any 'f' plaquettes
                    if (i < z_plaquettes.size()) {
                        if (z_plaquettes[i].get_shape() == 'f') {
                            std::cerr << "LogicalQubit::add_stabilizer: Stabilizers that anti-commute with any non-boundary plaquettes are not allowed to be added." << std::endl;
                            abort();
                        }
                    }
                    else if (x_plaquettes[i - z_plaquettes.size()].get_shape() == 'f') {
                        std::cerr << "LogicalQubit::add_stabilizer: Stabilizers that anti-commute with any non-boundary plaquettes are not allowed to be added." << std::endl;
                        abort();  
                    }

                    anticommuting_stabilizers.push_back(i);

                    // Print (for debugging purposes)
                    if (debug) {
                        std::cout << std::endl << "Anti-commuting stabilizer: ";
                        std::copy(parity_check_matrix.value()[i].begin(), parity_check_matrix.value()[i].end(), std::ostream_iterator<bool>(std::cout));
                        std::cout << std::endl;
                        std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(parity_check_matrix.value()[i]), true);
                        grid.print_grid(ascii_grid);
                    }
                }

                // If this row of parity_check_matrix corresponds with a logical operator
                else {
                    anticommuting_logical_ops.push_back(i);

                    // Print (for debugging purposes)
                    if (debug) {
                        std::cout << std::endl << "Anti-commuting logical operator: ";
                        std::copy(parity_check_matrix.value()[i].begin(), parity_check_matrix.value()[i].end(), std::ostream_iterator<bool>(std::cout));
                        std::cout << std::endl;
                        std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(parity_check_matrix.value()[i]), true);
                        grid.print_grid(ascii_grid);
                    }
                }
            }
        }

        // Before modifying the stabilizer vectors, ensure they pass appropriate checks
        test_stabilizers(); 

        // Create an optional alternative logical operator, to be used in cases where neither stored logical operator anti-commuted with the added stabilizer
        std::optional<std::vector<bool>> alternative_logical_operator;

        // Checks on results
        if (anticommuting_stabilizers.size() > 2) {
            std::cerr << "LogicalQubit::add_stabilizer: At most two stabilizers should anti-commute with a measured boundary stabilizer." << std::endl;
            abort();
        }

        else if (anticommuting_logical_ops.size() > 1) {
            std::cerr << "LogicalQubit::add_stabilizer: At most one stored (X or Z) logical operator should anti-commute with a measured boundary stabilizer." << std::endl;
            abort();
        }

        else if ((anticommuting_stabilizers.size() == 0) && (anticommuting_logical_ops.size() == 1)) {
            std::cerr << "LogicalQubit::add_stabilizer: Added stabilizer anti-commutes with a logical operator but with no stabilizers." << std::endl;
            abort();
        }

        else if ((anticommuting_stabilizers.size() == 0) && (anticommuting_logical_ops.size() == 0)) {
            std::cerr << "LogicalQubit::add_stabilizer: Added stabilizer already exists." << std::endl;
            abort();
        }

        else if ((anticommuting_stabilizers.size() == 2) && (anticommuting_logical_ops.size() == 1)) {
            std::cerr << "LogicalQubit::add_stabilizer: Strange situation encountered; two anti-commuting stabilizers and one anti-commuting logical operator." << std::endl;
            abort();            
        }

        // If neither of the stored logical operators anti-commute with the added stabilizer,
        // (depending on the case) we construct an equivalent operator along the opposite edge of the logical qubit.
        else if (anticommuting_logical_ops.size() == 0) {
            std::vector<bool> tmp_logical_operator(2*qsite_to_index->size());
            if (type == 'X') {
                tmp_logical_operator = parity_check_matrix.value()[parity_check_matrix->size()-2];
                for (unsigned int i=0; i<z_plaquettes.size(); i++) {
                    for (unsigned int k=0; k<tmp_logical_operator.size(); k++) {
                        tmp_logical_operator[k] = tmp_logical_operator[k] ^ parity_check_matrix.value()[i][k];
                    }
                }
            }
            else if (type == 'Z') {
                tmp_logical_operator = parity_check_matrix.value()[parity_check_matrix->size()-1];
                for (unsigned int i=0; i<x_plaquettes.size(); i++) {
                    for (unsigned int k=0; k<tmp_logical_operator.size(); k++) {
                        tmp_logical_operator[k] = tmp_logical_operator[k] ^ parity_check_matrix.value()[i + z_plaquettes.size()][k];
                    }
                }
            }
            alternative_logical_operator = std::move(tmp_logical_operator);
            
            // Visualize operator on grid (for debugging purposes)
            if (debug) {
                std::cout << std::endl << "Anti-commuting logical operator (topol. equiv. to the one stored): ";
                std::copy(alternative_logical_operator.value().begin(), alternative_logical_operator.value().end(), std::ostream_iterator<bool>(std::cout));
                std::cout << std::endl;
                std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(alternative_logical_operator.value()), true);
                grid.print_grid(ascii_grid);
            }
        }

        // In case an anti-commuting logical operator needs to be replaced using products with the anti-commuting stabilizers
        std::optional<std::vector<bool>> new_logical_operator;

        // In some cases, a qubit will need to be removed in order to maintain a single logical qubit and retain code distance
        std::optional<unsigned int> overlapping_index;

        // If there is one anti-commuting stabilizer, 
        if (anticommuting_stabilizers.size() == 1) {

            // If the added stabilizer hadn't anti-commuted with one of the logical operators, an alternative logical operator should have been constructed
            if (anticommuting_logical_ops.size() == 0) {
                assert(alternative_logical_operator.has_value());

                // Check if the newly constructed logical operator anti-commutes with the added stabilizer
                bin_sym_prod = 0;
                for (unsigned int k=0; k<alternative_logical_operator.value().size(); k++) {
                    bin_sym_prod ^= (alternative_logical_operator.value()[k] && new_row_transformed[k]);
                }

                if (!bin_sym_prod) {
                    std::cerr << "LogicalQubit::add_stabilizer: Strange situation encountered; one anti-commuting stabilizer and zero anti-commuting logical operators." << std::endl;
                    abort();
                }

                // If the check passes, there is nothing to do in updating the logical operator

            }

            // If there is an anti-commuting logical operator, take its product with the single anti-commuting stabilizer
            else {
                std::vector<bool> tmp_logical_operator(2*qsite_to_index->size());
                for (unsigned int k=0; k<tmp_logical_operator.size(); k++) {
                    tmp_logical_operator[k] = parity_check_matrix.value()[anticommuting_stabilizers[0]][k] ^ parity_check_matrix.value()[anticommuting_logical_ops[0]][k];
                }
                new_logical_operator = std::move(tmp_logical_operator);
            }
        }

        // If there are two anti-commuting stabilizers, 
        else if (anticommuting_stabilizers.size() == 2) {

            // Take their product
            std::vector<bool> tmp_logical_operator(2*qsite_to_index->size());
            for (unsigned int i = 0; i<2*qsite_to_index->size(); i++) {
                tmp_logical_operator[i] = parity_check_matrix.value()[anticommuting_stabilizers[0]][i] ^ parity_check_matrix.value()[anticommuting_stabilizers[1]][i];
            }

            // Find out whether there is an overlapping qubit between this operator and the stored logical operator
            for (unsigned int i = 0; i<2*qsite_to_index->size(); i++) {
                if (tmp_logical_operator[i] == parity_check_matrix.value()[parity_check_matrix->size() - 1 - (type=='X')][i]) {
                    if (!overlapping_index) 
                        overlapping_index = std::make_optional(i);
                    else {
                        std::cerr << "LogicalQubit::add_stabilizer: Multiple overlapping qsites found (path 1)." << std::endl;
                        abort();
                    }
                }
            }

            // Update the tmp_logical_operator if overlap was found
            if (overlapping_index.has_value()) {
                for (unsigned int i = 0; i<2*qsite_to_index->size(); i++) {
                    tmp_logical_operator[i] = tmp_logical_operator[i] ^ parity_check_matrix.value()[parity_check_matrix->size() - 1 - (type=='X')][i];
                }
                new_logical_operator = std::move(tmp_logical_operator);
            }  

            // If none was found, try the alternative logical operator
            else {
                for (unsigned int i = 0; i<2*qsite_to_index->size(); i++) {
                    if (tmp_logical_operator[i] == alternative_logical_operator.value()[i]) {
                        if (!overlapping_index) 
                            overlapping_index = std::make_optional(i);
                        else {
                            std::cerr << "LogicalQubit::add_stabilizer: Multiple overlapping qsites found (path 2)." << std::endl;
                            abort();
                        }
                    }
                }
            }  

            // Now, if there is still no overlapping qubit, then this is an invalid path
            if (!overlapping_index) {
                std::cerr << "LogicalQubit::add_stabilizer: Measured stabilizer does not correspond with corner movement." << std::endl;
                abort();
            }

        }

        // Visualize new logical operator on grid (for debugging purposes)
        if (new_logical_operator.has_value()) {
            if (debug) {
                std::cout << std::endl << "Logical operator replaced with: ";
                std::copy(new_logical_operator.value().begin(), new_logical_operator.value().end(), std::ostream_iterator<bool>(std::cout));
                std::cout << std::endl;
                std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(new_logical_operator.value()), true);
                grid.print_grid(ascii_grid);
            }
        }

        // Visualize qubit to remove (for debugging purposes)
        if (overlapping_index.has_value()) {
            if (debug) {
                std::cout << std::endl << "Qubit flagged for removal at qsite: " << overlapping_index.value() << std::endl;
                std::pair<unsigned int, char> single_qubit = std::make_pair(overlapping_index.value(), type);
                std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator({single_qubit}, true);
                grid.print_grid(ascii_grid);
            }
        }

        // Add the new stabilizer
        // if (type == 'Z') {
        //     z_plaquettes.push_back(std::move(new_stabilizer));
        //     parity_check_matrix->insert(parity_check_matrix->begin() + z_plaquettes.size() - 1, std::move(new_row));
        // }
        // else {
        //     x_plaquettes.push_back(std::move(new_stabilizer));
        //     parity_check_matrix->insert(parity_check_matrix->begin() + z_plaquettes.size() + x_plaquettes.size() - 1, std::move(new_row));
        // }

        // // Loop over anti-commuting operators
        // for (unsigned int i=0; i<indices.size(); i++) {

        //     // Remove the stabilizers that anticommute with the new one
        //     if (i < parity_check_matrix->size() - 2) {
        //         if (i < z_plaquettes.size()) {
        //             z_plaquettes.erase(z_plaquettes.begin() + i);
        //         }
        //         else {
        //             x_plaquettes.erase(x_plaquettes.begin() + i - z_plaquettes.size());
        //         }
        //         parity_check_matrix->erase(parity_check_matrix->begin() + i);
        //     }

        //     // Update the logical operators
        //     else {

        //     }
        // }

        test_stabilizers(); 
        // Next we need to test

        /* We need to make the following updates: 
            - Replace the row(s) of the parity check matrix corresponding with anti-commuting operators
            - Remove any anti-commuting stabilizers from the plaquette vectors and add the new one
            - Make sure that we still have the same number of encoded qubits i.e. that the number of stabilizers still equals the number of qubits minus 1
            - Ensure that the correspondence in ordering of the parity check matrix rows vs. plaquette vector indices is preserved
            - Any removed boundary stabilizers may also need qsites corresp. with 'm' qubits to be de-occupied within the grid
            - A function that notes any decrease in code distance would be very useful
            - Do the logical operators need to be re-calculated?
            - Check that three-qubit stabilizers still work 

        */
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
                 << " " << p.grid()->index_from_coords(p.get_row(), p.get_col(), 1) << std::endl;
        }

        // Print all x stabilizers
        for (const Plaquette& p : x_plaquettes) {
            std::cout << p.get_type() << " " << p.get_row() << " " << p.get_col() << " " << p.get_shape()
                 << " " << p.grid()->index_from_coords(p.get_row(), p.get_col(), 1) << std::endl;
        }
    }

    /* Construct parity check matrix:
        - Rows refer to stabilizers and columns refer to qsites (repeated twice)
        - The final two rows are logical operators (Z and X, respectively)
    */
    void LogicalQubit::construct_parity_check_matrix(const GridManager& grid) {

        // Construct qsites to indices map
        std::map<unsigned int, unsigned int> qsite_to_index_;
        std::set<unsigned int> data = data_qsites();
        std::vector<unsigned int> index_to_qsite_(2*data.size());
        unsigned int counter = 0;
        for (unsigned int site: data) {
            qsite_to_index_[site] = counter;
            index_to_qsite_[counter] = site;
            index_to_qsite_[counter + data.size()] = site;
            counter++;
        }
        qsite_to_index = std::move(qsite_to_index_);
        index_to_qsite = std::move(index_to_qsite_);

        // Construct parity check matrix to represent both stabilizers and observables
        std::vector<std::vector<bool>> parity_check_(x_plaquettes.size() + z_plaquettes.size() + 2, std::vector<bool>(2*data.size(), 0));

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

    // Check validity of parity check matrix 
    bool LogicalQubit::validity_parity_check_matrix() {

        // Make sure that a parity check matrix has been constructed
        if (!parity_check_matrix) {
            return false;
        }

        /* 
            - To check validity, we calculate the upper triangle of a mtx where entry (i, j) is the binary symplectic product between row i and row j of the parity check mtx.
            - This product being zero or one indicates an even or odd number of anti-commuting pairs of single-qubit Paulis in the rows' corresponding Pauli strings.
            - Every entry of the resulting matrix should be zero
            - We also check validity for the logical operators, which are appended to the parity check matrix, and should anticommute with one another but commute with everything else
        */

        bool bin_sym_prod;
        std::vector<bool> tmp_row;
        tmp_row.reserve(2*qsite_to_index->size());
        for (unsigned int i=0; i<parity_check_matrix->size(); i++) {
            for (unsigned int j=i+1; j<parity_check_matrix->size(); j++) {
                bin_sym_prod = 0;
                for (unsigned int k=0; k<qsite_to_index->size(); k++) {
                    tmp_row[k] = parity_check_matrix.value()[j][k+qsite_to_index->size()];
                    tmp_row[k+qsite_to_index->size()] = parity_check_matrix.value()[j][k];
                }
                for (unsigned int k = 0; k < parity_check_matrix.value()[i].size(); k++) {
                    bin_sym_prod ^= (parity_check_matrix.value()[i][k] && tmp_row[k]);
                }

                // The product between the last two rows should be 1 since the Z and X logical operators should anticommute
                if ((i == parity_check_matrix->size() - 2) && (j == parity_check_matrix->size() - 1)) {
                    if (!bin_sym_prod) return bin_sym_prod;
                }

                // Every other pair of rows should have a product of 0
                else {
                    if (bin_sym_prod) return !bin_sym_prod;
                }
            }
        }

        return bin_sym_prod;
    }


    // Print parity check matrix and logical operators
    void LogicalQubit::print_parity_check_matrix(const GridManager& grid) {

        // If the matrix hasn't been constructed yet, construct it
        if (!parity_check_matrix) {
            construct_parity_check_matrix(grid);
        }

        // Print map from column indices to qsites
        std::cout << "Map from column indices to qsites:" << std::endl;
        for (std::pair<unsigned int, unsigned int> pair : qsite_to_index.value()) {
            std::cout << pair.second << " " << pair.first << std::endl;
        }
        for (std::pair<unsigned int, unsigned int> pair : qsite_to_index.value()) {
            std::cout << qsite_to_index->size() + pair.second << " " << pair.first << std::endl;
        }

        // Print parity check matrix and logical operators
        std::cout << std::endl;
        std::cout << "Parity check matrix and logical operators (rows in same order as given by print_stabilizers function):" << std::endl;
        for (unsigned int i = 0; i<parity_check_matrix->size(); i++) {
            if (i == parity_check_matrix->size() - 2) {
                for (unsigned int j=0; j<2*qsite_to_index->size()+1; j++) {
                    std::cout << "-";
                }
                std::cout << std::endl;
            }
            for (unsigned int j = 0; j<qsite_to_index->size(); j++) {
                std::cout << parity_check_matrix.value()[i][j];
            }
            std::cout << "|";
            for (unsigned int j = 0; j<qsite_to_index->size(); j++) {
                std::cout << parity_check_matrix.value()[i][j + qsite_to_index->size()];
            }
            std::cout << std::endl;
        }

        // Check validity of parity check matrix
        bool validity = validity_parity_check_matrix();
        std::cout << std::endl << "Parity check matrix is valid: " << validity << std::endl;
    }

    // Test stabilizers (not fully implemented)
    void LogicalQubit::test_stabilizers() {
        unsigned int num_data_qubits = dx_*dz_;
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
        test_stabilizers();
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