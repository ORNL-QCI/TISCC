#include <TISCC/logicalqubit.hpp>

#include <iomanip>
#include <cassert>
#include <algorithm>
#include <limits>
#include <iterator>
#include <utility>

namespace TISCC 
{
    /* Helper functions */
    // Helper function to calculate the dot product mod 2 of two binary vectors
    bool bin_dot_prod_mod_2(const std::vector<bool>& v1, const std::vector<bool> v2) {
        if (v1.size() != v2.size()) {
            std::cerr << "bin_sym_prod: vectors of unequal length given to bin_sym_prod." << std::endl;
            abort();
        }
        bool bin_prod = 0;
        for (unsigned int k = 0; k < v1.size(); k++) {
            bin_prod ^= (v1[k] && v2[k]);
        }
        return bin_prod;
    }

    // Helper function to calculate the product of two logical operators expressed in binary symplectic format
    std::pair<std::vector<bool>, int> operator_product_binary_format(const std::vector<bool>& v1, const std::vector<bool> v2) {
        if (v1.size() != v2.size()) {
            std::cerr << "operator_product_binary_format: vectors of unequal length given to bin_sym_prod." << std::endl;
            std::cerr << v1.size() << " " << v2.size() << std::endl;
            abort();
        }
        std::vector<bool> operator_product_binary_format(v1.size());
        for (unsigned int k = 0; k < v1.size(); k++) {
            operator_product_binary_format[k] = (v1[k] ^ v2[k]);
        }
        // We employ a convention where Z operators act to the left of X operators and v1 acts to the left of v2.
        // Thus, a minus sign results from there being, at a given site, an X on v1 and a Z on v2 
        int sign = 1;
        for (unsigned int k = 0; k < v1.size()/2; k++) {
            if (v2[k] && v1[k + v1.size()/2]) {
                sign *= -1;
            }
        }

        return std::make_pair(operator_product_binary_format, sign);
    }

    // Helper function to symplectic transform a vector in binary symplectic format
    std::vector<bool> symplectic_transform(const std::vector<bool>& v) {
        if (v.size()%2 != 0) {
            std::cerr << "symplectic_transform: input vector has odd size." << std::endl;
            abort();
        }
        std::vector<bool> v_prime(v.size());
        for (unsigned int k=0; k<v.size()/2; k++) {
            v_prime[k] = v[k+v.size()/2];
            v_prime[k+v.size()/2] = v[k];
        }
        return v_prime;
    }

    // Helper function to return the Pauli weight of an operator expressed in binary symplectic format
    // Note: this currently adds a weight of two if there is any ZX on the same site
    unsigned int pauli_weight(const std::vector<bool>& v) {
        unsigned int weight = 0;
        for (bool elem : v) weight += elem;
        return weight;
    }

    // Transform operators from binary representation to string e.g. 11000101 = Z(ZX)IX = i*ZYIX 
    std::pair<std::string, std::complex<double>> binary_operator_to_pauli_string(const std::vector<bool>& binary_rep) {
        std::complex<double> phase = std::complex<double>(1.0, 0.0);
        std::string result;
        for (unsigned int k = 0; k < binary_rep.size()/2; k++) {
            if (binary_rep[k]) {
                result.push_back('Z');
            }
            else {
                result.push_back('I');
            }
        }
        for (unsigned int k = 0; k < binary_rep.size()/2; k++) {
            if ((binary_rep[binary_rep.size()/2 + k]) && (result[k] == 'I')) {
                result[k] = 'X';
            }
            else if ((binary_rep[binary_rep.size()/2 + k]) && (result[k] == 'Z')) {
                result[k] = 'Y';
                phase *= std::complex<double>(0.0, 1.0);
            }
        }
        return std::make_pair(result, phase);
    }

    // Transform operators from string to binary representation e.g. ZYIX -> -i*(11000101)
    std::pair<std::vector<bool>, std::complex<double>> pauli_string_to_binary_operator(const std::string& pauli_string) {
        std::complex<double> phase = std::complex<double>(1.0, 0.0);
        std::vector<bool> result(2*pauli_string.size(), 0);
        for (unsigned int k = 0; k < pauli_string.size(); k++) {
            if (pauli_string[k] == 'Z') {
                result[k] = 1;
            }
            else if (pauli_string[k] == 'X') {
                result[k+pauli_string.size()] = 1;
            }
            else if (pauli_string[k] == 'Y') {
                result[k] = 1;
                result[k+pauli_string.size()] = 1;
                phase *= std::complex<double>(0.0, -1.0);
            }
            else if (pauli_string[k] != 'I') {
                std::cerr << "pauli_string_to_binary_vector: invalid character found." << std::endl;
                abort();
            }
        }
        return std::make_pair(result, phase);
    }

    /* LogicalQubit member functions */

    /* Stabilizer-related member functions */
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

        // Set default circuit patterns
        for (Plaquette& p : z_plaquettes) {
            p.set_circuit_pattern('Z');
        }
        for (Plaquette& p : x_plaquettes) {
            p.set_circuit_pattern('N');
        }
    }

    // Test stabilizers
    void LogicalQubit::test_stabilizers() {
        double uint_max = std::numeric_limits<unsigned int>::max();

        // Get number of data qubits and number of stabilizers
        unsigned int num_data_qubits = occupied_sites(true).size();
        unsigned int num_stabilizers = z_plaquettes.size() + x_plaquettes.size();

        if (num_data_qubits - num_stabilizers != 1) {
            std::cerr << "LogicalQubit::test_stabilizers: To encode one logical qubit, there should be one fewer stabilizer than data qubit. Currently there are :" << std::endl;
            std::cerr << num_data_qubits << " data qubits and " << num_stabilizers << " stabilizers." << std::endl;
            abort();
        }

        // Check that all plaquette types are correct
        for (const Plaquette& p : z_plaquettes) {
            if (p.get_operator_type() != 'Z') {
                std::cerr << "LogicalQubit::test_stabilizers: Wrong plaquette type in z_plaquettes." << std::endl;
                abort();
            }
        }
        for (const Plaquette& p : x_plaquettes) {
            if (p.get_operator_type() != 'X') {
                std::cerr << "LogicalQubit::test_stabilizers: Wrong plaquette type in x_plaquettes." << std::endl;
                abort();
            }
        }

        // Create vector of all plaquettes
        std::vector<Plaquette> all_plaquettes;
        all_plaquettes.insert(all_plaquettes.end(), z_plaquettes.begin(), z_plaquettes.end());
        all_plaquettes.insert(all_plaquettes.end(), x_plaquettes.begin(), x_plaquettes.end());

        // Check that measure qubits are all unique
        std::vector<unsigned int> measure_qubits;
        for (Plaquette p : all_plaquettes) {
            for (unsigned int q : measure_qubits) {
                if (p.get_qsite('m') == q) {std::cerr << "LogicalQubit::test_stabilizers: Found duplicate measure qubit among plaquettes."; abort();}
            }
            measure_qubits.push_back(p.get_qsite('m'));
        }

        // Check consistency of stabilizer type with its qsites
        for (Plaquette p: all_plaquettes) { 
            bool n_invalid = (p.get_shape() == 'n' && ((p.get_qsite('a') != uint_max) || (p.get_qsite('b') != uint_max) || (p.get_qsite('c') == uint_max) || (p.get_qsite('d') == uint_max)));
            bool s_invalid = (p.get_shape() == 's' && ((p.get_qsite('c') != uint_max) || (p.get_qsite('d') != uint_max) || (p.get_qsite('a') == uint_max) || (p.get_qsite('b') == uint_max)));
            bool e_invalid = (p.get_shape() == 'e' && ((p.get_qsite('b') != uint_max) || (p.get_qsite('d') != uint_max) || (p.get_qsite('a') == uint_max) || (p.get_qsite('c') == uint_max)));
            bool w_invalid = (p.get_shape() == 'w' && ((p.get_qsite('a') != uint_max) || (p.get_qsite('c') != uint_max) || (p.get_qsite('b') == uint_max) || (p.get_qsite('d') == uint_max)));

            if (n_invalid || s_invalid || e_invalid || w_invalid) {
                std::cerr << "LogicalQubit::test_stabilizers: Inconsistent stabilizer found." << std::endl; abort();
            }
        }

        // Ensure consistency with rows of parity check mtx
        for (unsigned int i=0; i<parity_check_matrix.size() - 2; i++) {
            bool consistent = 1;
            unsigned int weight = 0;
            if (i < z_plaquettes.size()) {
                for (char qubit : {'a', 'b', 'c', 'd'}) {
                    if (z_plaquettes[i].get_qsite(qubit) != uint_max) {
                        weight++;
                        consistent = consistent && parity_check_matrix[i][qsite_to_index[z_plaquettes[i].get_qsite(qubit)]];
                    }
                }
            }
            else {
                for (char qubit : {'a', 'b', 'c', 'd'}) {
                    if (x_plaquettes[i - z_plaquettes.size()].get_qsite(qubit) != uint_max) {
                        weight++;
                        consistent = consistent && parity_check_matrix[i][qsite_to_index[x_plaquettes[i - z_plaquettes.size()].get_qsite(qubit)] + qsite_to_index.size()];
                    }
                } 
            }
            if (weight != pauli_weight(parity_check_matrix[i])) {consistent = 0;}
            if (!consistent) {
                std::cerr << "LogicalQubit::test_stabilizers: stabilizer row " << i << " out of " << parity_check_matrix.size() - 2 << " in parity check matrix inconsistent with corresponding stabilizer." << std::endl;
                abort();
            }
        }
    }

    void LogicalQubit::print_stabilizers() const {

        // Print all z stabilizers
        for (const Plaquette& p : z_plaquettes) {
            std::cout << p.get_operator_type() << " " << p.get_row() << " " << p.get_col() << " " << p.get_shape()
                 << " " << p.grid()->index_from_coords(p.get_row(), p.get_col(), 1) << std::endl;
        }

        // Print all x stabilizers
        for (const Plaquette& p : x_plaquettes) {
            std::cout << p.get_operator_type() << " " << p.get_row() << " " << p.get_col() << " " << p.get_shape()
                 << " " << p.grid()->index_from_coords(p.get_row(), p.get_col(), 1) << std::endl;
        }
    }

    /* Parity-check matrix related member functions */
    // Transform operators from binary representation to pair<qsite unsigned int, Pauli char>
    // ** Note: this function will convert ZX directly to Y without tracking any phase.
    std::vector<std::pair<unsigned int, char>> LogicalQubit::binary_operator_to_qsites(const std::vector<bool>& binary_rep) {
        std::vector<std::pair<unsigned int, char>> qsite_rep;
        for (unsigned int k = 0; k < qsite_to_index.size(); k++) {
            if (binary_rep[k] && binary_rep[k+qsite_to_index.size()]) {
                qsite_rep.emplace_back(index_to_qsite[k], 'Y');
            }
            else if (binary_rep[k]) {
                qsite_rep.emplace_back(index_to_qsite[k], 'Z');
            }
            else if (binary_rep[k+qsite_to_index.size()]) {
                qsite_rep.emplace_back(index_to_qsite[k], 'X');
            }
        }
        return qsite_rep;
    }

    /* Construct parity check matrix:
        - Rows refer to stabilizers and columns refer to qsites (repeated twice)
        - The final two rows are logical operators (Z and X, respectively)
    */
    void LogicalQubit::construct_parity_check_matrix(const GridManager& grid) {

        // Construct qsites to indices map
        std::map<unsigned int, unsigned int> qsite_to_index_;
        std::set<unsigned int> data = occupied_sites(true);
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
                if (qsite_to_index.count(p.get_qsite(qubit))) 
                    parity_check_[counter][qsite_to_index[p.get_qsite(qubit)]] = 1;
            }
            counter++;
        }

       // Iterate over x stabilizers
        for (const Plaquette& p : x_plaquettes) {
            for (char qubit : {'a', 'b', 'c', 'd'}) {
                if (qsite_to_index.count(p.get_qsite(qubit))) 
                    parity_check_[counter][data.size() + qsite_to_index[p.get_qsite(qubit)]] = 1;
            }
            counter++;
        }

        // Obtain observables
        for (unsigned int qsite : data) {
            // The col_'th column of the grid contains the left-most column of data qubits contained in the surface code patch
            if (grid.get_col(qsite) == col_) {
                parity_check_[counter][qsite_to_index[qsite]] = 1;
            }
            // The (row_+1)'th row of the grid contains the top-most row of data qubits contained in the surface code patch
            if (grid.get_row(qsite) == row_+1) {
                parity_check_[counter+1][data.size() + qsite_to_index[qsite]] = 1;
            }
        }

        parity_check_matrix = std::move(parity_check_);
    }

    // Check validity of parity check matrix 
    bool LogicalQubit::validity_parity_check_matrix() const {
        /* 
            - To check validity, we calculate the upper triangle of a mtx where entry (i, j) is the binary symplectic product between row i and row j of the parity check mtx.
            - This product being zero or one indicates an even or odd number of anti-commuting pairs of single-qubit Paulis in the rows' corresponding Pauli strings.
            - Every entry of the resulting matrix should be zero
            - We also check validity for the logical operators, which are appended to the parity check matrix, and should anticommute with one another but commute with everything else
        */

        bool bin_sym_prod;
        std::vector<bool> tmp_row(2*qsite_to_index.size());
        for (unsigned int i=0; i<parity_check_matrix.size(); i++) {
            for (unsigned int j=i+1; j<parity_check_matrix.size(); j++) {
                for (unsigned int k=0; k<qsite_to_index.size(); k++) {
                    tmp_row[k] = parity_check_matrix[j][k+qsite_to_index.size()];
                    tmp_row[k+qsite_to_index.size()] = parity_check_matrix[j][k];
                }
                bin_sym_prod = bin_dot_prod_mod_2(parity_check_matrix[i], tmp_row);

                // The product between the last two rows should be 1 since the Z and X logical operators should anticommute
                if ((i == parity_check_matrix.size() - 2) && (j == parity_check_matrix.size() - 1)) {
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
    void LogicalQubit::print_parity_check_matrix() const {

        // Print map from column indices to qsites
        std::cout << "Map from column indices to qsites:" << std::endl;
        for (std::pair<unsigned int, unsigned int> pair : qsite_to_index) {
            std::cout << pair.second << " " << pair.first << std::endl;
        }
        for (std::pair<unsigned int, unsigned int> pair : qsite_to_index) {
            std::cout << qsite_to_index.size() + pair.second << " " << pair.first << std::endl;
        }

        // Print parity check matrix and logical operators
        std::cout << std::endl;
        std::cout << "Parity check matrix and logical operators (rows in same order as given by print_stabilizers function):" << std::endl;
        for (unsigned int i = 0; i<parity_check_matrix.size(); i++) {
            if (i == parity_check_matrix.size() - 2) {
                for (unsigned int j=0; j<2*qsite_to_index.size()+1; j++) {
                    std::cout << "-";
                }
                std::cout << std::endl;
            }
            for (unsigned int j = 0; j<qsite_to_index.size(); j++) {
                std::cout << parity_check_matrix[i][j];
            }
            std::cout << "|";
            for (unsigned int j = 0; j<qsite_to_index.size(); j++) {
                std::cout << parity_check_matrix[i][j + qsite_to_index.size()];
            }
            std::cout << std::endl;
        }

        // Check validity of parity check matrix
        bool validity = validity_parity_check_matrix();
        std::cout << std::endl << "Parity check matrix is valid: " << validity << std::endl;
    }

    /* Other initializations and constructor */
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

    LogicalQubit::LogicalQubit(unsigned int dx, unsigned int dz, unsigned int row, unsigned int col, GridManager& grid) : 
        dx_(dx), dx_init_(dx), dz_(dz), dz_init_(dz), row_(row), col_(col), default_arrangement_(true), TI_model() { 
        init_stabilizers(dx, dz, row, col, grid);
        construct_parity_check_matrix(grid);
        init_circuits();
        test_stabilizers();
    }

    /* Derived quantities and accessors */
    // Default edge logical operator
    const std::vector<bool>& LogicalQubit::get_logical_operator_default_edge(char type) const {
        return parity_check_matrix[parity_check_matrix.size() - 1 - (type == 'Z')];
    }

    // Opposite edge logical operator
    std::vector<bool> LogicalQubit::get_logical_operator_opposite_edge(char type) const {
        std::vector<bool> tmp_logical_operator = get_logical_operator_default_edge(type);
        if (type == 'Z') {
            for (unsigned int i=0; i<z_plaquettes.size(); i++) {
                tmp_logical_operator = operator_product_binary_format(tmp_logical_operator, parity_check_matrix[i]).first;
            }
        }
        else if (type == 'X') {
            for (unsigned int i=0; i<x_plaquettes.size(); i++) {
                tmp_logical_operator = operator_product_binary_format(tmp_logical_operator, parity_check_matrix[i + z_plaquettes.size()]).first;
            }
        }
        return tmp_logical_operator;
    }

    // A product between stabilizer measurements taken at these qsites with exp. val. of the current logical op yields correct logical exp. val.
    std::vector<unsigned int> LogicalQubit::get_logical_deformation_qsites(char type) const {
        if (type == 'X') {
            return x_deformation_qsites;
        }
        else if (type == 'Z') {
            return z_deformation_qsites;
        }
        else {std::cerr << "LogicalQubit::get_logical_deformation_qsites: invalid input given." << std::endl; abort();}
    }

    // Function to return all qsites occupied by the surface code (or just_data_qubits)
    std::set<unsigned int> LogicalQubit::occupied_sites(bool just_data_qubits) const {
        std::set<unsigned int> sites;
        for (char qubit : {'a', 'b', 'c', 'd', 'm'}) {
            if ((just_data_qubits) && (qubit=='m')) {
                continue;
            }
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

    // Recalculate code distance using logical operators
    void LogicalQubit::recalculate_code_distance() {
        dx_ = std::min(pauli_weight(get_logical_operator_default_edge('X')), pauli_weight(get_logical_operator_opposite_edge('X')));
        dz_ = std::min(pauli_weight(get_logical_operator_default_edge('Z')), pauli_weight(get_logical_operator_opposite_edge('Z')));
    }

    // Update the logical deformation vectors with qsite
    void LogicalQubit::add_logical_deformation_qsites(char type, unsigned int qsite) {
        if (type == 'X') {x_deformation_qsites.push_back(qsite);}
        else if (type == 'Z') {z_deformation_qsites.push_back(qsite);}
        else {std::cerr << "LogicalQubit::add_logical_deformation_qsites: invalid input given." << std::endl; abort();}
    }

    // Apply a given ``qubit-level'' instruction to all plaquettes in a given vector and add corresponding HW_Instructions to hw_master
    double LogicalQubit::apply_instruction(const Instruction& instr, Plaquette& p, double time, unsigned int step, 
        const GridManager& grid, std::vector<HW_Instruction>& hw_master) {

        // Create tmp variables for time
        // ** Initialization to zero causes apply_instruction to return 0 if an Idle was given and an updated time in any other case
        double time_tmp;
        double time_to_return = time;

        // Don't explicitly apply an "Idle" operation
        if (instr.get_name() != "Idle") {

            // Add single-qubit instructions and track updated time
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

            // Track the time update corresp. to the plaquette for which the instruction took the longest
            if (time_tmp > time_to_return) time_to_return = time_tmp;

        }

        return time_to_return;
    }

    double LogicalQubit::idle(unsigned int cycles, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {

        // The four circuit types should have the same number of instructions and can be run in parallel on the applicable stabilizers
        if (!(Z_Circuit_Z_Type.size() == X_Circuit_N_Type.size()) && (Z_Circuit_Z_Type.size() == Z_Circuit_N_Type.size())
            && (Z_Circuit_Z_Type.size() == X_Circuit_Z_Type.size())) {
                std::cerr << Z_Circuit_Z_Type.size() << " " << X_Circuit_N_Type.size() << " " << Z_Circuit_N_Type.size() << " " << X_Circuit_Z_Type.size() << std::endl;
                abort();
        }
        unsigned int num_instructions = Z_Circuit_Z_Type.size();

        // Create a vector of all plaquettes for convenience
        std::vector<Plaquette> all_plaquettes;
        all_plaquettes.insert(all_plaquettes.end(), z_plaquettes.begin(), z_plaquettes.end());
        all_plaquettes.insert(all_plaquettes.end(), x_plaquettes.begin(), x_plaquettes.end());

        // Loop over surface code cycles
        for (unsigned int cycle=0; cycle < cycles; cycle++) {

            // Loop over Instructions 
            for (unsigned int i=0; i<num_instructions; i++) {

                // Loop over all plaquettes and apply instruction
                double latest_time_returned = time;
                double time_tmp;
                for (Plaquette& p : all_plaquettes) {

                    if (p.get_operator_type() == 'Z' && p.get_circuit_pattern() == 'Z') {
                        time_tmp = apply_instruction(Z_Circuit_Z_Type[i], p, time, i, grid, hw_master);
                    }

                    else if (p.get_operator_type() == 'Z' && p.get_circuit_pattern() == 'N') {
                        time_tmp = apply_instruction(Z_Circuit_N_Type[i], p, time, i, grid, hw_master);
                    }

                    else if (p.get_operator_type() == 'X' && p.get_circuit_pattern() == 'Z') {
                        time_tmp = apply_instruction(X_Circuit_Z_Type[i], p, time, i, grid, hw_master);
                    }

                    else if (p.get_operator_type() == 'X' && p.get_circuit_pattern() == 'N') {
                        time_tmp = apply_instruction(X_Circuit_N_Type[i], p, time, i, grid, hw_master);
                    }

                    else {
                        std::cerr << "LogicalQubit::idle: Unrecognized operator type and/or circuit pattern." << std::endl;
                        abort();
                    }

                    if (time_tmp > latest_time_returned) {latest_time_returned = time_tmp;}

                }

                time = latest_time_returned;

            }
            
        }

        // Sort master list of hardware instructions according to overloaded operator<
        std::stable_sort(hw_master.begin(), hw_master.end());

        return time;

    }

    double LogicalQubit::transversal_op(const std::string& op, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {
        
        // Get all occupied sites
        std::set<unsigned int> sites = occupied_sites(true);

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

    double LogicalQubit::apply_pauli(char pauli, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {

        // Obtain appropriate logical operator
        std::vector<bool> logical_operator;
        if (pauli == 'Z') {logical_operator = get_logical_operator_default_edge('Z');}
        else if (pauli == 'X') {logical_operator = get_logical_operator_default_edge('X');}
        else if (pauli == 'Y') {logical_operator = operator_product_binary_format(get_logical_operator_default_edge('Z'), get_logical_operator_default_edge('X')).first;}
        else {std::cerr << "LogicalQubit::apply_pauli: Invalid Pauli operator given." << std::endl; abort();}

        // Loop over sites on the operator and apply Pauli operators. 
        // ** Y is covered by applying Z_{pi/2}*X_{pi/2}. (Z_{pi/2}*X_{pi/2} = -i*Z -i*X = - Z*X = -iY = -Y_{pi/2}) 
        std::vector<double> times(qsite_to_index.size(), time);
        for (unsigned int i=logical_operator.size()-1; i < logical_operator.size(); i--) {
            if ((i < qsite_to_index.size()) && (logical_operator[i])) {
                times[i] = TI_model.add_Z(index_to_qsite[i], times[i], 0, grid, hw_master);
            }

            else if ((i >= qsite_to_index.size()) && (logical_operator[i])) {
                times[i-qsite_to_index.size()] = TI_model.add_X(index_to_qsite[i], times[i-qsite_to_index.size()], 0, grid, hw_master);
            }
        }

        // Sort master list of hardware instructions according to overloaded operator<
        std::stable_sort(hw_master.begin(), hw_master.end());

        return *std::max_element(times.begin(), times.end());
    }

    // This is not a good state injection protocol; it is just an easy way to generate the positive Y eigenstate in the case of zero noise 
    double LogicalQubit::inject_y_state(const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {
        double time_tmp;
        time = transversal_op("prepz", grid, hw_master, time);
        std::vector<bool> logical_operator = get_logical_operator_default_edge('X');
        for (unsigned int i=logical_operator.size()/2; i<logical_operator.size(); i++) {
            if (logical_operator[i]) {
                time_tmp = TI_model.add_H(index_to_qsite[i], time, 0, grid, hw_master);
            }
        }
        time = time_tmp;
        // We locate the site on which Y_L supports a Y. Note Z_{pi/4} |+> = e^{-i*pi/4} S |+> = e^{-i*pi/4} * |Y, +>. We ignore the global phase.
        auto y_string = binary_operator_to_pauli_string(operator_product_binary_format(get_logical_operator_default_edge('Z'), get_logical_operator_default_edge('X')).first);
        bool y_found = 0;
        for (unsigned int i=0; i<y_string.first.size(); i++) {
            if (y_string.first[i] == 'Y') {
                if (y_found == 1) {
                    std::cerr << "LogicalQubit::inject_y_state: found > 1 y operator in y_string." << std::endl;
                    abort();
                }
                time_tmp = TI_model.add_Z_pi4(index_to_qsite[i], time, 0, grid, hw_master);
                y_found = 1;
            }
        }

        return time_tmp;
    }

    // Helper function to return the data qubits from this patch that are NOT occupied by two others
    std::set<unsigned int> LogicalQubit::get_strip(LogicalQubit& lq1, LogicalQubit& lq2) {
        std::set<unsigned int> data = occupied_sites(true);
        std::set<unsigned int> lq1_data = lq1.occupied_sites(true);
        std::set<unsigned int> lq2_data = lq2.occupied_sites(true);
        std::set<unsigned int> strip;
        for (unsigned int lq_site : data) {
            if ((lq1_data.find(lq_site) == lq1_data.end()) && (lq2_data.find(lq_site) == lq2_data.end())) {
                strip.insert(lq_site);
            }
        }
        return strip;
    }

    // Obtain pair<qsite unsigned int, Pauli char> for each stabilizer measure qubit for the purpose of labeling on the grid
    std::vector<std::pair<unsigned int, char>> LogicalQubit::syndrome_measurement_qsites() {
        std::vector<Plaquette> all_plaquettes;
        all_plaquettes.insert(all_plaquettes.end(), z_plaquettes.begin(), z_plaquettes.end());
        all_plaquettes.insert(all_plaquettes.end(), x_plaquettes.begin(), x_plaquettes.end());

        std::vector<std::pair<unsigned int, char>> syndrome_measurement_qsites;
        for (const Plaquette& p : all_plaquettes) {
            syndrome_measurement_qsites.emplace_back(p.get_qsite('m'), p.get_operator_type());
        }
        return syndrome_measurement_qsites;
    }

    // Swap roles of x and z for this patch (used during Hadamard and patch rotation)
    void LogicalQubit::xz_swap(const GridManager& grid) {

        // Construct new parity check matrix
        std::vector<std::vector<bool>> new_parity_check_matrix;
        std::vector<bool> tmp_row;
        for (unsigned int i = 0; i<x_plaquettes.size(); i++) {
            tmp_row = symplectic_transform(parity_check_matrix[i + z_plaquettes.size()]);
            new_parity_check_matrix.push_back(tmp_row);
        }
        for (unsigned int i=0; i<z_plaquettes.size(); i++) {
            tmp_row = symplectic_transform(parity_check_matrix[i]);
            new_parity_check_matrix.push_back(tmp_row);
        }
        new_parity_check_matrix.push_back(symplectic_transform(parity_check_matrix[parity_check_matrix.size() - 1]));
        new_parity_check_matrix.push_back(symplectic_transform(parity_check_matrix[parity_check_matrix.size() - 2]));   
        parity_check_matrix = std::move(new_parity_check_matrix);
        assert(validity_parity_check_matrix()==1);

        // Transform stabilizers while keeping circuit patterns the same
        std::vector<Plaquette> new_z_plaquettes;
        std::vector<Plaquette> new_x_plaquettes;
        for (Plaquette& p : z_plaquettes) {
            p.change_operator_type('X');
            new_x_plaquettes.push_back(p);
        }
        for (Plaquette& p : x_plaquettes) {
            p.change_operator_type('Z');
            new_z_plaquettes.push_back(p);
        }
        z_plaquettes = std::move(new_z_plaquettes);
        x_plaquettes = std::move(new_x_plaquettes);
        test_stabilizers();

        // Update code distances
        recalculate_code_distance();

        // We have changed the stabilizer arrangement
        default_arrangement_ = false;

    }

    // The result of this process is the same as if we flipped the patch upside down and then did xz_swap
    float LogicalQubit::flip_patch(GridManager& grid, std::vector<HW_Instruction> hw_master, float time, bool debug) {

        // If the stabilizers have been altered at all, don't allow merge
        if (!default_arrangement_) {
            std::cerr << "LogicalQubit::flip_patch: flip_patch not allowed for qubits that do not have the default stabilizer arrangement." << std::endl;
            abort();
        }

        // Recall that order matters
        extend_logical_operator_clockwise('X', "opposite", dz_init_ - 1, grid, hw_master, time, debug);
        extend_logical_operator_clockwise('Z', "opposite", dx_init_ - 1, grid, hw_master, time, debug); 
        extend_logical_operator_clockwise('X', "default", dz_init_ - 1, grid, hw_master, time, debug); 
        extend_logical_operator_clockwise('Z', "default", dx_init_ - 1, grid, hw_master, time, debug);       

        // We have changed the stabilizer arrangement
        default_arrangement_ = false;

        return time;

    }

    // Translate patch s rows "South" or e rows "East" on the underlying grid
    double LogicalQubit::translate_patch(int s, int e, const GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) {

        std::cerr << "LogicalQubit::translate_patch: not fully implemented." << std::endl;
        abort();

        // Handle trivial case
        if ((s == 0) && (e == 0 )) return time;

        // Rows and columns needed for the patch itself
        // **Note: In the standard configuration, we need one extra row and column than the code distance to avoid interference between boundary stabilizers of adjacent patches
        int patch_rows = dz_init_ + 1;
        int patch_cols = dx_init_ + 1;

        // Check that minimum number of rows and columns needed to execute translation is satisfied by the grid
        unsigned int min_rows = row_ + patch_rows;
        unsigned int min_cols = col_ + patch_cols;
        if (s > 0) {min_rows += s;}
        if (e > 0) {min_cols += e;}

        if ((min_rows > grid.get_nrows()) || min_cols > grid.get_ncols()) {
            std::cerr << "LogicalQubit::translate_patch: Ordered patch translation requires larger grid." << std::endl; abort();
        }

        // Check that negative translations are possible given grid constraints
        if (((row_ + s) < 0) || ((col_ + e) < 0)) {
            std::cerr << "LogicalQubit::translate_patch: Ordered patch translation requires larger grid." << std::endl; abort();
        }  

        /* For now, we only allow e == 1 and s == 0 */
        if (!((e==1)&&(s==0))) {std::cerr << "LogicalQubit::translate_patch: only implemented for one shift rightwards." << std::endl; abort();}



        // /* Update plaquettes and parity check matrix */
        // std::vector<Plaquette> all_plaquettes;
        // all_plaquettes.insert(all_plaquettes.end(), z_plaquettes.begin(), z_plaquettes.end());
        // all_plaquettes.insert(all_plaquettes.end(), x_plaquettes.begin(), x_plaquettes.end());

        // // Loop over data qsites
        // for (auto& data_qsite_index_pair : qsite_to_index) {

        //     // Get stabilizer indices for this qsite
        //     std::vector<unsigned int> stabilizer_indices;

        //     // Find all the stabilizers that have support on this qsite
        //     unsigned int i=0; 
        //     for (const std::vector<bool>& stabilizer_vec : parity_check_matrix) {
        //         if (stabilizer_vec[data_qsite_index_pair.second] || stabilizer_vec[data_qsite_index_pair.second + z_plaquettes.size()]) {
        //             stabilizer_indices.push_back(i);
        //         }
        //         i++;
        //     }

        //     // Update plaquettes
        //     unsigned int new_qsite = grid.shift_qsite(data_qsite_index_pair.first, s, e);
        //     for (unsigned int index : stabilizer_indices) {
        //         bool found=0;
        //         for (char qubit : {'a', 'b', 'c', 'd'}) {
        //             if (all_plaquettes[index].get_qsite(qubit) == data_qsite_index_pair.first) {
        //                 found = 1;
        //                 all_plaquettes[index].move_to_site(qubit, new_qsite);
        //             }
        //         } 
        //     }

        //     // 
        //     data_qsite_index_pair.first = new_qsite;

        // }

        // // Retrieve set of qsites occupied by this patch
        // std::set<unsigned int> occ_sites = occupied_sites(false);


        // // 

        // // Incorporate necessary Move operations on the Grid


        // // Make necessary changes to stabilizers


        // for (Plaquette& p : all_plaquettes) {
        //     for (char qubit : {'a', 'b', 'c', 'd'}) {
        //         p.move_to_site(qubit, grid.shift_qsite(p.get_qsite(qubit), s, e));
        //     } 
        // }

        // // Make necessary changes to parity check matrix and related variables
        return time;
    }

   // Add new stabilizer plaquette (to be used in corner movement)
    double LogicalQubit::add_stabilizer(unsigned int row, unsigned int col, char shape, char type, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time, bool debug) {

        // We require stabilizers of a known type (X or Z)
        char opp_type;
        if (type == 'Z') {
            opp_type = 'X';
        }
        else if (type == 'X') {
            opp_type = 'Z';
        }
        else {
            std::cerr << "LogicalQubit::add_stabilizer: Invalid stabilizer type (X or Z) given." << std::endl;
            abort(); 
        }

        // We require newly measured stabilizers to be along the boundaries and of the appropriate shape, since the purpose of this function is to be used in corner movement
        if (!(((col == col_) && (row > row_) && (row < row_ + dz_init_) && (shape == 'w')) ||  // Left edge
            ((col == col_ + dx_init_) && (row > row_) && (row < row_ + dz_init_) && (shape == 'e')) || // Right edge
            ((row == row_) && (col > col_) && (col < col_ + dx_init_) && (shape == 'n')) ||  // Top edge
            ((row == row_ + dz_init_) && (col > col_) && (col < col_ + dx_init_) && (shape == 's')))) { // Bottom edge
            std::cout << "LogicalQubit::add_stabilizer: Newly added plaquettes must lie along boundaries and be of the appropriate shape, as this function is only meant to be used in corner movement." << std::endl;
            abort();
        }

        // In some cases, a hardware instruction is added, so time needs to be tracked
        double time_tmp = time;

        // Create vector of all plaquettes for convenience
        std::vector<Plaquette> all_plaquettes;
        all_plaquettes.insert(all_plaquettes.end(), z_plaquettes.begin(), z_plaquettes.end());
        all_plaquettes.insert(all_plaquettes.end(), x_plaquettes.begin(), x_plaquettes.end());

        /* Add the stabilizer to the grid:
            - Depending on the case, a data qubit need to be added. 
            - We tease out that case using differences in the set of occupied sites.
        */

        // Create optional added_qubit and set of (just data) qsites previously occupied by this lq 
        std::optional<unsigned int> added_qubit;
        std::set<unsigned int> previously_occupied_sites = occupied_sites(true);

        // Construct new Plaquette from the grid without explicitly adding it yet
        // (note that more validity checks take place within this function)
        // ** The set of occupied_sites on grid is updated to include any new qubits
        Plaquette new_stabilizer(grid.get_plaquette(row, col, shape, type)); 

        // Update added_qubit if needed
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();  
        for (char qubit : {'a', 'b', 'c', 'd'}) {
            if ((new_stabilizer.get_qsite(qubit) != uint_max) && !previously_occupied_sites.count(new_stabilizer.get_qsite(qubit))) {
                if (!added_qubit.has_value()) {
                    added_qubit = new_stabilizer.get_qsite(qubit);
                }
                else {
                    std::cerr << "LogicalQubit::add_stabilizer: more than 1 data qubit from the added stabilizer was not previously loaded on the grid." << std::endl;
                    abort();
                }
            }
        }

        // Make sure that it is represented in the qsites_to_indices map
        if (added_qubit.has_value() && !qsite_to_index.count(added_qubit.value())) {
            std::cerr << "LogicalQubit::add_stabilizer: Qubit to be added not represented in qsites_to_indices map." << std::endl;
            abort();
        }

        // Construct a new vector in binary symplectic format corresponding to this stabilizer
        std::vector<bool> new_row(2*qsite_to_index.size(), 0);
        for (char qubit : {'a', 'b', 'c', 'd'}) {
            if (qsite_to_index.count(new_stabilizer.get_qsite(qubit))) {
                if (type == 'Z') {
                    new_row[qsite_to_index[new_stabilizer.get_qsite(qubit)]] = 1;
                }
                else if (type == 'X') {
                    new_row[qsite_to_index.size() + qsite_to_index[new_stabilizer.get_qsite(qubit)]] = 1;
                }
            }
        }

        /* Check commutation of new stabilizer with every row of parity_check_matrix (and appended logical operators) */

        // First, transform the new row s.t. the roles of X and Z have been swapped
        std::vector<bool> new_row_transformed = symplectic_transform(new_row);

        // Then, calculate the binary product (mod 2) with every row of the parity_check_matrix and track row indices for which it is 1
        std::vector<unsigned int> anticommuting_stabilizers;
        std::vector<unsigned int> anticommuting_logical_ops;
        for (unsigned int i=0; i<parity_check_matrix.size(); i++) {

            // If the binary symplectic product = 1, then the two operators anti-commute
            if (bin_dot_prod_mod_2(parity_check_matrix[i], new_row_transformed)) {

                // If this row of parity_check_matrix corresponds with a stabilizer plaquette,
                if (i < parity_check_matrix.size() - 2) {

                    // We don't allow the new stabilizer to anti-commute with any 'f' plaquettes
                    if (all_plaquettes[i].get_shape() == 'f') {

                        // If we are adding a qubit, this would not yet be reflected in the parity check mtx.
                        if (added_qubit.has_value()) {

                            // Calculate binary symplectic product with updated stabilizer; update parity check mtx if necessary
                            std::vector<bool> plaquette_with_added_qubit = parity_check_matrix[i];
                            plaquette_with_added_qubit[qsite_to_index[added_qubit.value()] + (type == 'Z')*qsite_to_index.size()] = 1;
                            if (bin_dot_prod_mod_2(plaquette_with_added_qubit, new_row_transformed)) {
                                std::cerr << "LogicalQubit::add_stabilizer: Stabilizers that anti-commute with any non-boundary plaquettes are not allowed to be added." << std::endl;
                                abort();   
                                // throw std::runtime_error("LogicalQubit::add_stabilizer: Stabilizers that anti-commute with any non-boundary plaquettes are not allowed to be added.");
                            }
                            else {

                                // Update parity check matrix and stabilizer vector
                                parity_check_matrix[i] = plaquette_with_added_qubit;
                                if (i < z_plaquettes.size()) {
                                    z_plaquettes[i] = grid.get_plaquette(z_plaquettes[i].get_row(), z_plaquettes[i].get_col(), 'f', 'Z');
                                }
                                else {
                                    x_plaquettes[i - z_plaquettes.size()] = grid.get_plaquette(x_plaquettes[i - z_plaquettes.size()].get_row(), x_plaquettes[i - z_plaquettes.size()].get_col(), 'f', 'X');
                                }
                            }
                        }

                        else {
                            std::cerr << "LogicalQubit::add_stabilizer: Stabilizers that anti-commute with any non-boundary plaquettes are not allowed to be added." << std::endl;
                            abort(); 
                            // throw std::runtime_error("LogicalQubit::add_stabilizer: Stabilizers that anti-commute with any non-boundary plaquettes are not allowed to be added."); 
                        }
  
                    }

                    else {
                        anticommuting_stabilizers.push_back(i);
                    }
                }

                // If this row of parity_check_matrix corresponds with a logical operator
                else {

                    // If we are adding a qubit, this would not yet be reflected in the parity check mtx.
                    if (added_qubit.has_value()) {

                        // Calculate binary symplectic product with updated logical; update parity check mtx if logical now commutes
                        std::vector<bool> logical_with_added_qubit = parity_check_matrix[i];
                        logical_with_added_qubit[qsite_to_index[added_qubit.value()] + (type == 'Z')*qsite_to_index.size()] = 1;
                        if (!bin_dot_prod_mod_2(logical_with_added_qubit, new_row_transformed)) {
                            parity_check_matrix[i] = logical_with_added_qubit;
                        }
                        else {
                            anticommuting_logical_ops.push_back(i);
                        }
                    }

                    else {
                        anticommuting_logical_ops.push_back(i);
                    }
                }
            }
        }

        // Print stabilizers/qubits to be added/removed and logical ops to be updated
        if (debug) {
            std::cout << std::endl << "Stabilizer being added: ";
            std::copy(new_row.begin(), new_row.end(), std::ostream_iterator<bool>(std::cout));
            std::cout << std::endl;
            std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(new_row), true);
            grid.print_grid(ascii_grid);

            if (added_qubit.has_value()) {
                std::cout << std::endl << "Qubit to be added at qsite: " << added_qubit.value() << std::endl;
                std::pair<unsigned int, char> single_qubit = std::make_pair(added_qubit.value(), opp_type);
                ascii_grid = grid.ascii_grid_with_operator({single_qubit}, true);
                grid.print_grid(ascii_grid);
            }

            for (unsigned int index: anticommuting_stabilizers) {
                std::cout << std::endl << "Stabilizer being removed: ";
                std::copy(parity_check_matrix[index].begin(), parity_check_matrix[index].end(), std::ostream_iterator<bool>(std::cout));
                std::cout << std::endl;
                ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(parity_check_matrix[index]), true);
                grid.print_grid(ascii_grid);
            }

            for (unsigned int index: anticommuting_logical_ops) {
                std::cout << std::endl << "Anti-commuting logical operator to be updated: ";
                std::copy(parity_check_matrix[index].begin(), parity_check_matrix[index].end(), std::ostream_iterator<bool>(std::cout));
                std::cout << std::endl;
                ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(parity_check_matrix[index]), true);
                grid.print_grid(ascii_grid);
            }
        }

        // Create an optional alternative logical operator, to be used in cases where neither stored logical operator anti-commuted with the added stabilizer
        std::optional<std::vector<bool>> logical_operator_opposite_edge;

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
            if (!added_qubit.has_value()) {
                std::cerr << "LogicalQubit::add_stabilizer: Added stabilizer already exists." << std::endl;
                abort();
            }
        }

        else if ((anticommuting_stabilizers.size() == 2) && (anticommuting_logical_ops.size() == 1)) {
            std::cerr << "LogicalQubit::add_stabilizer: Strange situation encountered; two anti-commuting stabilizers and one anti-commuting logical operator." << std::endl;
            abort();            
        }

        // If neither of the stored logical operators anti-commute with the added stabilizer,
        // (depending on the case) we construct an equivalent operator along the opposite edge of the logical qubit.
        else if (anticommuting_logical_ops.size() == 0) {
            std::vector<bool> tmp_logical_operator = get_logical_operator_opposite_edge(opp_type);
            logical_operator_opposite_edge = std::move(tmp_logical_operator);
        }

        // In case an anti-commuting logical operator needs to be replaced using products with the anti-commuting stabilizers
        std::optional<std::vector<bool>> new_logical_operator_default_edge;
        std::optional<std::vector<bool>> new_logical_operator_opposite_edge;

        // In some cases, a qubit will need to be removed in order to maintain a single logical qubit and retain code distance
        std::optional<unsigned int> overlapping_index;

        // If there is one anti-commuting stabilizer, 
        if (anticommuting_stabilizers.size() == 1) {

            // If the added stabilizer hadn't anti-commuted with one of the logical operators, an alternative logical operator should have been constructed
            if (anticommuting_logical_ops.size() == 0) {
                assert(logical_operator_opposite_edge.has_value());

                // Check if the newly constructed logical operator anti-commutes with the added stabilizer
                if (!bin_dot_prod_mod_2(logical_operator_opposite_edge.value(), new_row_transformed)) {
                    std::cerr << "LogicalQubit::add_stabilizer: Strange situation encountered; one anti-commuting stabilizer and zero anti-commuting logical operators." << std::endl;
                    abort();
                }

                // If the check passes, there is nothing to do in updating the logical operator

            }

            // If there is an anti-commuting logical operator along the default edge, take its product with the single anti-commuting stabilizer (both assumed to ONLY support Z or ONLY support X)
            else {
                std::vector<bool> tmp_logical_operator = operator_product_binary_format(parity_check_matrix[anticommuting_stabilizers[0]], parity_check_matrix[anticommuting_logical_ops[0]]).first;
                new_logical_operator_default_edge = std::move(tmp_logical_operator);
                add_logical_deformation_qsites(opp_type, all_plaquettes[anticommuting_stabilizers[0]].get_qsite('m'));
            }
        }

        // If there are two anti-commuting stabilizers, 
        else if (anticommuting_stabilizers.size() == 2) {

            // Take their product (both assumed to ONLY support Z or ONLY support X)
            std::vector<bool> tmp_logical_operator = operator_product_binary_format(parity_check_matrix[anticommuting_stabilizers[0]], parity_check_matrix[anticommuting_stabilizers[1]]).first;

            // if (debug) {
            //     std::cout << std::endl << "New logical operator comprised of two anti-commuting stabilizers: ";
            //     std::copy(tmp_logical_operator.begin(), tmp_logical_operator.end(), std::ostream_iterator<bool>(std::cout));
            //     std::cout << std::endl;
            //     std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(tmp_logical_operator), true);
            //     grid.print_grid(ascii_grid);
            // }

            // Find out whether there is an overlapping qubit between this operator and the stored logical operator
            for (unsigned int i = 0; i<2*qsite_to_index.size(); i++) {
                if (tmp_logical_operator[i] && (parity_check_matrix[parity_check_matrix.size() - 1 - (type=='X')][i])) {
                    if (!overlapping_index) 
                        overlapping_index = std::make_optional(index_to_qsite[i]);
                    else {
                        std::cerr << "LogicalQubit::add_stabilizer: Multiple overlapping qsites found (path 1)." << std::endl;
                        abort();
                    }
                }
            }

            // Update the tmp_logical_operator if overlap was found (both assumed to ONLY support Z or ONLY support X)
            if (overlapping_index.has_value()) {
                tmp_logical_operator = operator_product_binary_format(tmp_logical_operator, parity_check_matrix[parity_check_matrix.size() - 1 - (type=='X')]).first;
                new_logical_operator_default_edge = std::move(tmp_logical_operator);
                add_logical_deformation_qsites(opp_type, all_plaquettes[anticommuting_stabilizers[0]].get_qsite('m'));
                add_logical_deformation_qsites(opp_type, all_plaquettes[anticommuting_stabilizers[1]].get_qsite('m'));

                // Visualize new logical operator on grid (for debugging purposes)
                if (debug) {
                    std::cout << std::endl << "Default edge logical operator replaced with: ";
                    std::copy(new_logical_operator_default_edge.value().begin(), new_logical_operator_default_edge.value().end(), std::ostream_iterator<bool>(std::cout));
                    std::cout << std::endl;
                    std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(new_logical_operator_default_edge.value()), true);
                    grid.print_grid(ascii_grid);
                }
            } 

            // If none was found, try the opposite edge logical operator
            else {
                for (unsigned int i = 0; i<2*qsite_to_index.size(); i++) {
                    if (tmp_logical_operator[i] && (logical_operator_opposite_edge.value()[i])) {
                        if (!overlapping_index) 
                            overlapping_index = std::make_optional(index_to_qsite[i]);
                        else {
                            std::cerr << "LogicalQubit::add_stabilizer: Multiple overlapping qsites found (path 2)." << std::endl;
                            abort();
                        }
                    }
                }

                // Update the tmp_logical_operator if overlap was found (both assumed to ONLY support Z or ONLY support X)
                if (overlapping_index.has_value()) {
                    tmp_logical_operator = operator_product_binary_format(tmp_logical_operator, logical_operator_opposite_edge.value()).first;
                    new_logical_operator_opposite_edge = std::move(tmp_logical_operator);

                    // Visualize new logical operator on grid (for debugging purposes)
                    if (debug) {
                        std::cout << std::endl << "Opposite edge logical operator replaced with: ";
                        std::copy(new_logical_operator_opposite_edge.value().begin(), new_logical_operator_opposite_edge.value().end(), std::ostream_iterator<bool>(std::cout));
                        std::cout << std::endl;
                        std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(new_logical_operator_opposite_edge.value()), true);
                        grid.print_grid(ascii_grid);
                    }
                } 
            }

            // Now, if there is still no overlapping qubit, then this is an invalid path
            if (!overlapping_index) {
                std::cerr << "LogicalQubit::add_stabilizer: Measured stabilizer does not correspond with corner movement." << std::endl;
                abort();
                // throw std::runtime_error("LogicalQubit::add_stabilizer: Measured stabilizer does not correspond with corner movement.");
            }

        }

        // If the default-edge logical operator of the same type as the added stabilizer has overlap with it, it should be replaced (both assumed to ONLY support Z or ONLY support X)
        std::optional<std::vector<bool>> new_same_type_logical_operator;
        std::vector<bool> tmp_row = operator_product_binary_format(parity_check_matrix[parity_check_matrix.size() - 1 - (type=='Z')], new_row).first;
        if (pauli_weight(tmp_row) <= pauli_weight(parity_check_matrix[parity_check_matrix.size() - 1 - (type=='Z')])) {
            new_same_type_logical_operator = std::move(tmp_row);
            add_logical_deformation_qsites(type, new_stabilizer.get_qsite('m'));
        }

        // Replace the old logical operators with new ones in the parity check mtx
        if (new_logical_operator_default_edge.has_value()) {
            parity_check_matrix[parity_check_matrix.size() - 1 - (type=='X')] = new_logical_operator_default_edge.value();
        }
        if (new_same_type_logical_operator.has_value()) {
            parity_check_matrix[parity_check_matrix.size() - 1 - (type=='Z')] = new_same_type_logical_operator.value();
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

        // Remove the overlapping qubit if necessary
        if (overlapping_index.has_value()) {
            if (type=='X') {
                for (unsigned int i=0; i<x_plaquettes.size(); i++) {
                    if (x_plaquettes[i].remove_supported_qsite(overlapping_index.value())) {
                        parity_check_matrix[i + z_plaquettes.size()][qsite_to_index[overlapping_index.value()] + qsite_to_index.size()] = 0;
                    }
                }
                time_tmp = TI_model.add_H(overlapping_index.value(), time, 0, grid, hw_master);
                time_tmp = TI_model.add_meas(overlapping_index.value(), time_tmp, 1, grid, hw_master);
            }
            else {
                for (unsigned int i=0; i<z_plaquettes.size(); i++) {
                    if (z_plaquettes[i].remove_supported_qsite(overlapping_index.value())) {
                        parity_check_matrix[i][qsite_to_index[overlapping_index.value()]] = 0;
                    }
                }
                time_tmp = TI_model.add_meas(overlapping_index.value(), time, 0, grid, hw_master);
            }

            if (parity_check_matrix[parity_check_matrix.size() - 1 - (type=='Z')][qsite_to_index[overlapping_index.value()] + (type=='X')*qsite_to_index.size()]) {
                parity_check_matrix[parity_check_matrix.size() - 1 - (type=='Z')][qsite_to_index[overlapping_index.value()] + (type=='X')*qsite_to_index.size()] = 0;
                add_logical_deformation_qsites(type, overlapping_index.value());
            }

        }

        // If we had to add a qubit, add appropriate hardware instructions to initialize it
        if (added_qubit.has_value()) {
            time_tmp = TI_model.add_init(added_qubit.value(), time, 0, grid, hw_master);
            if (type=='Z') {
                time_tmp = TI_model.add_H(added_qubit.value(), time, 1, grid, hw_master);
            }
        }

        // Remove anti-commuting stabilizers from plaquette vectors
        for (unsigned int i=0; i<anticommuting_stabilizers.size(); i++) {
            if (anticommuting_stabilizers[i] < z_plaquettes.size()) {
                z_plaquettes.erase(z_plaquettes.begin() + anticommuting_stabilizers[i]);   
            }
            else {
                x_plaquettes.erase(x_plaquettes.begin() + anticommuting_stabilizers[i] - z_plaquettes.size());
            }
            parity_check_matrix.erase(parity_check_matrix.begin() + anticommuting_stabilizers[i]);

            for (unsigned int j=i+1; j<anticommuting_stabilizers.size(); j++) {
                if (anticommuting_stabilizers[j] > anticommuting_stabilizers[i]) 
                    anticommuting_stabilizers[j]--;
            }
        }

        // Finally, add the new stabilizer
        if (type == 'Z') {
            z_plaquettes.push_back(std::move(new_stabilizer));
            parity_check_matrix.insert(parity_check_matrix.begin() + z_plaquettes.size() - 1, std::move(new_row));
        }
        else {
            x_plaquettes.push_back(std::move(new_stabilizer));
            parity_check_matrix.insert(parity_check_matrix.begin() + z_plaquettes.size() + x_plaquettes.size() - 1, std::move(new_row));
        }

        // Test validity of stabilizers and parity check mtx
        test_stabilizers(); 
        if (!validity_parity_check_matrix()) {
            std::cerr << "LogicalQubit::add_stabilizer: Final parity check matrix invalid." << std::endl;
            abort();
        }

        // Print final grid and final default-edge logical operators
        if (debug) {
            std::cout << std::endl << "Final grid: " << std::endl;
            std::vector<std::string> ascii_grid = grid.ascii_grid(true);
            grid.print_grid(ascii_grid);

            std::cout << std::endl << "Final default-edge logical Z operator: ";
            std::copy(get_logical_operator_default_edge('Z').begin(), get_logical_operator_default_edge('Z').end(), std::ostream_iterator<bool>(std::cout));
            std::cout << std::endl;
            ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(get_logical_operator_default_edge('Z')), true);
            grid.print_grid(ascii_grid);

            std::cout << std::endl << "Final default-edge logical X operator: ";
            std::copy(get_logical_operator_default_edge('X').begin(), get_logical_operator_default_edge('X').end(), std::ostream_iterator<bool>(std::cout));
            std::cout << std::endl;
            ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(get_logical_operator_default_edge('X')), true);
            grid.print_grid(ascii_grid);

            std::cout << std::endl << "Final opposite-edge logical Z operator: ";
            std::vector<bool> opp_edge_z = get_logical_operator_opposite_edge('Z');
            std::copy(opp_edge_z.begin(), opp_edge_z.end(), std::ostream_iterator<bool>(std::cout));
            std::cout << std::endl;
            ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(opp_edge_z), true);
            grid.print_grid(ascii_grid);

            std::cout << std::endl << "Final opposite-edge logical X operator: ";
            std::vector<bool> opp_edge_x = get_logical_operator_opposite_edge('X');
            std::copy(opp_edge_x.begin(), opp_edge_x.end(), std::ostream_iterator<bool>(std::cout));
            std::cout << std::endl;
            ascii_grid = grid.ascii_grid_with_operator(binary_operator_to_qsites(opp_edge_x), true);
            grid.print_grid(ascii_grid);
        }

        // Set this variable to indicate that the stabilizer arrangement has been altered
        default_arrangement_ = false;

        // Recalculate code distance
        recalculate_code_distance();
        if (debug) {
            std::cout << std::endl << "New code distances: dx = " << dx_ << " and dz = " << dz_ << "." << std::endl;
        }

        return time_tmp;
    }
    
    // Implements corner movements by figuring out which stabilizers to measure when extending a logical operator `clockwise'
    /* **Note: add_stabilizer will update the default (rather than opposite) logical operator that has support on the added stabilizer in ambiguous cases,   
        which may become relevant in the case where the operator to be extended has support on a single qubit */
    double LogicalQubit::extend_logical_operator_clockwise(char type, std::string_view edge_type, unsigned int weight_to_add, 
        GridManager& grid, std::vector<HW_Instruction>& hw_master, double time, bool debug) {

        // We require stabilizers of a known type (X or Z)
        char opp_type;
        if (type == 'Z') {
            opp_type = 'X';
        }
        else if (type == 'X') {
            opp_type = 'Z';
        }
        else {
            std::cerr << "LogicalQubit::add_stabilizer: Invalid stabilizer type (X or Z) given." << std::endl;
            abort(); 
        }

        // Hardware instructions may be added in add_stabilizer, so we need to track time
        double time_tmp;
        double time_to_return = time;

        // Grab logical operators of 'type' on default and opposite edge
        std::vector<bool> logical_operator;
        if (edge_type == "default") {logical_operator = get_logical_operator_default_edge(type);}
        else if (edge_type == "opposite") {logical_operator = get_logical_operator_opposite_edge(type);}
        else {std::cerr << "LogicalQubit::extend_logical_operator_clockwise: invalid input for edge_type." << std::endl; abort();}
        unsigned int max_weight_to_add = 2*(dx_init_ - 1) + 2*(dz_init_ - 1) - pauli_weight(logical_operator);

        // Find any qsite supported on the logical operator
        std::optional<unsigned int> qsite;
        for (unsigned int i=0; i<logical_operator.size(); i++) {
            if (logical_operator[i]) {
                if (grid.is_occupied(index_to_qsite[i])) { 
                    qsite = index_to_qsite[i];
                    break;
                }
                else {
                    std::cerr << "LogicalQubit::extend_logical_operator_clockwise: Logical operator supports qsite that isn't occupied on the grid." << std::endl;
                    abort();
                }

            } 
        }

        if (!qsite.has_value()) {
            std::cerr << "LogicalQubit::extend_logical_operator_clockwise: Logical operator in question has no support on any qsites." << std::endl;
            abort();
        }

        // Find the grid spec of this qsite
        unsigned int row = grid.get_row(qsite.value());
        unsigned int col = grid.get_col(qsite.value());
        unsigned int index = grid.get_idx(qsite.value());

        // Loop over single-qubit operators to add in extending the logical operator
        for (unsigned int i=0; i<weight_to_add; i++) {

            // Don't let it try to add above the maximum Pauli weight before the operator eats its tail
            if (i >= max_weight_to_add) {
                break;
            }

            // We make sure we aren't adding a qsite already supported on the logical
            unsigned int pc_column_index = qsite_to_index[grid.index_from_coords(row, col, index)] + (type == 'X')*qsite_to_index.size();
            while (logical_operator[pc_column_index]) {

                // Figure out which boundary the present qsite lies on. Assign corners assuming movement will be clockwise. Set coords to the next (occupied) clockwise boundary qsite.
                if ((col == col_) && (row != row_ + 1)) {
                    // left edge 
                    row--;
                    if (!grid.is_occupied(grid.index_from_coords(row, col, index))) {
                        col++;
                    }
                }

                else if ((row == row_ + 1) && (col != col_ + dx_init_ - 1)) {
                    // top edge
                    col++;
                    if (!grid.is_occupied(grid.index_from_coords(row, col, index))) {
                        row++;
                    }
                }

                else if ((col == col_ + dx_init_ - 1) && (row != row_ + dz_init_)) {
                    // right edge
                    row++;
                    if (!grid.is_occupied(grid.index_from_coords(row, col, index))) {
                        col--;
                    }
                }

                else if ((row == row_ + dz_init_) && (col != col_)) {
                    // bottom edge
                    col--;
                    if (!grid.is_occupied(grid.index_from_coords(row, col, index))) {
                        row--;
                    }
                }

                else {
                    std::cerr << "LogicalQubit::extend_logical_operator_clockwise: Logical operator has support on non-boundary qubit." << std::endl;
                    abort();
                }

                qsite = grid.index_from_coords(row, col, index);
                pc_column_index = qsite_to_index[qsite.value()] + (type == 'X')*qsite_to_index.size();

            }

            // Create vector of all plaquettes for convenience
            std::vector<Plaquette> all_plaquettes;
            all_plaquettes.insert(all_plaquettes.end(), z_plaquettes.begin(), z_plaquettes.end());
            all_plaquettes.insert(all_plaquettes.end(), x_plaquettes.begin(), x_plaquettes.end());

            // Find any boundary stabilizer that has on this qsite corresponding to 'type'
            std::optional<unsigned int> supporting_stabilizer_index;
            for (unsigned int i=0; i<all_plaquettes.size(); i++) {
                if ((parity_check_matrix[i][pc_column_index]) &&
                    (all_plaquettes[i].get_shape() != 'f')) {
                    supporting_stabilizer_index = i;
                    break;
                }
            }

            if (!supporting_stabilizer_index.has_value()) {
                std::cerr << "LogicalQubit::extend_logical_operator_clockwise: No supporting boundary stabilizer found at (qsite, pc_column): (" << qsite.value() << ", " << pc_column_index << ")" << std::endl;
                abort();              
            }

            // Set up basic parameters for input to add_stabilizer
            char new_stab_type = opp_type;
            char new_stab_shape = all_plaquettes[supporting_stabilizer_index.value()].get_shape();
            unsigned int new_stab_row = all_plaquettes[supporting_stabilizer_index.value()].get_row();
            unsigned int new_stab_col = all_plaquettes[supporting_stabilizer_index.value()].get_col();
            int shift_horizontal = 0;
            int shift_vertical = 0;
            if (new_stab_shape == 'n') {shift_horizontal = 1;}
            else if (new_stab_shape == 'e') {shift_vertical = 1;}
            else if (new_stab_shape == 's') {shift_horizontal = -1;}
            else if (new_stab_shape == 'w') {shift_vertical = -1;}
            else {std::cerr << "LogicalQubit::extend_logical_operator_clockwise: invalid stabilizer shape." << std::endl;}

            // We will handle two cases:
            //  (1) The found boundary stabilizer has no support on the current logical
            //  (2) The found boundary stabilizer has one qubit of support on the current logical

            // Grab the other qsite associated with this boundary stabilizer
            std::optional<unsigned int> other_qsite;
            unsigned int pc_column_index_other;
            for (char qubit : {'a', 'b', 'c', 'd'}) {
                if ((all_plaquettes[supporting_stabilizer_index.value()].get_qsite(qubit) != std::numeric_limits<unsigned int>::max()) &&
                    all_plaquettes[supporting_stabilizer_index.value()].get_qsite(qubit) != qsite.value()) {
                    if (!other_qsite.has_value()) {
                        other_qsite = all_plaquettes[supporting_stabilizer_index.value()].get_qsite(qubit);
                    }
                    else {
                        for (char qubit : {'a', 'b', 'c', 'd', 'm'}) {
                            std::cout << qubit << " " << all_plaquettes[supporting_stabilizer_index.value()].get_qsite(qubit) << std::endl;
                        }
                        std::cerr << "LogicalQubit::extend_logical_operator_clockwise: Supporting stabilizer invalid (path 1)." << std::endl;
                        abort();
                    }
                }
            }
            
            if (!other_qsite.has_value()) {
                std::cerr << "LogicalQubit::extend_logical_operator_clockwise: Supporting stabilizer invalid (path 2)." << std::endl;
            }

            else {
                pc_column_index_other = qsite_to_index[other_qsite.value()] + (type=='X')*qsite_to_index.size();
            }

            // Check whether logical operator has support here and add_stabilizer depending on the case
            if (logical_operator[pc_column_index_other]) {
                time_tmp = add_stabilizer(new_stab_row+shift_vertical, new_stab_col+shift_horizontal, new_stab_shape, new_stab_type, grid, hw_master, time, debug);           
            }
            else {
                time_tmp = add_stabilizer(new_stab_row-shift_vertical, new_stab_col-shift_horizontal, new_stab_shape, new_stab_type, grid, hw_master, time, debug); 
            }

            /* Something to note; it is unclear whether there will be cases where qubits will need to be added and removed in the same corner movement. 
            In that case, will have measure and init of the same qubit at the same time in the hardware circuit. */

            // Update logical operator
            unsigned int weight_diff = pauli_weight(logical_operator);
            if (edge_type == "default") {logical_operator = get_logical_operator_default_edge(type);}
            else if (edge_type == "opposite") {logical_operator = get_logical_operator_opposite_edge(type);}
            weight_diff = pauli_weight(logical_operator) - weight_diff;
            i += weight_diff - 1;
            if (weight_diff == 0) {
                std::cerr << "LogicalQubit::extend_logical_operator_clockwise: logical operator of desired edge type did not change weight. add_stabilizer may have chosen the default edge operator to modify." << std::endl;
            }

            // Print for debugging purposes
            if (debug) {
                std::vector<std::string> ascii_grid = grid.ascii_grid_with_operator(syndrome_measurement_qsites(), true);
                grid.print_grid(ascii_grid);
            }

          if (time_tmp > time_to_return) time_to_return = time_tmp;

        }

        return time_to_return;

    }

    // Reset stabilizer circuits to default values
    void LogicalQubit::reset_stabilizer_circuit_patterns() {
        // Set default circuit patterns
        for (Plaquette& p : z_plaquettes) {
            p.set_circuit_pattern('Z');
        }
        for (Plaquette& p : x_plaquettes) {
            p.set_circuit_pattern('N');
        }
    }

    // Construct and return a logical qubit that represents the merged product of two input logical qubits
    LogicalQubit* merge(LogicalQubit& lq1, LogicalQubit& lq2, GridManager& grid) {

        // If the stabilizers have been altered at all, don't allow merge
        if ((!lq1.default_arrangement()) || (!lq2.default_arrangement())) {
            std::cerr << "merge: Merge not allowed for qubits that do not have the default stabilizer arrangement." << std::endl;
            abort();
        }

        // Determine whether to merge horizontally or vertically and set parameters

        // If they are horizontally displaced,
        if (lq1.get_row() == lq2.get_row()) {

            // We require them to have the same code distance
            assert(lq1.get_dz() == lq2.get_dz());

            // Merge horizontally
            unsigned int extra_strip = 0;
            if (lq1.get_dx()%2 == 0) {extra_strip = 1;}
            assert(lq2.get_col() == lq1.get_col() + lq1.get_dx() + 1 + extra_strip);
            return new LogicalQubit(lq1.get_dx() + lq2.get_dx() + 1 + extra_strip, lq1.get_dz(), lq1.get_row(), lq1.get_col(), grid);

        }

        // Otherwise, if they are vertically displaced,
        else if (lq1.get_col() == lq2.get_col()) {

            // We require them to have the same x code distance
            assert(lq1.get_dx() == lq2.get_dx());

            // Merge vertically
            unsigned int extra_strip = 0;
            if (lq1.get_dz()%2 == 0) {extra_strip = 1;} 
            assert(lq2.get_row() == lq1.get_row() + lq1.get_dz() + 1 + extra_strip);
            return new LogicalQubit(lq1.get_dx(), lq1.get_dz() + lq2.get_dz() + 1 + extra_strip, lq1.get_row(), lq1.get_col(), grid); 

        }

        else {
            std::cerr << "merge: this operation must take place between logical qubits either vertically or horizontally separated, but not both." << std::endl;
            abort();
        }



    }
}