#include <TISCC/hardwaremodel.hpp>

#include <limits>
#include <cassert>
#include <algorithm>


namespace TISCC 
{
    // Initialize hash table to map native TI operations to times (in microseconds)
    void HardwareModel::init_TI_ops() {

        // Single-site ops
        TI_ops["Prepare_Z"] = 10;
        TI_ops["Measure_Z"] = 120;
        TI_ops["X_pi/2"] = 10;
        TI_ops["Y_pi/2"] = 10;
        TI_ops["Z_pi/2"] = 3;
        TI_ops["X_pi/4"] = 10;
        TI_ops["Y_pi/4"] = 10;
        TI_ops["Z_pi/4"] = 3;
        TI_ops["Z_pi/8"] = 3;
        TI_ops["X_-pi/4"] = 10;
        TI_ops["Y_-pi/4"] = 10;
        TI_ops["Z_-pi/4"] = 3;
        TI_ops["Z_-pi/8"] = 3;

        // Move is currently per-site for two types of sites: trapping zones and junctions
        double move_velocity = 80; // m/s, see https://arxiv.org/pdf/2301.05279.pdf
        double junction_move_velocity = 4; // m/s, see https://arxiv.org/pdf/2301.05279.pdf (it seems 4 m/s is safe but up to 6 m/s is possible)
        TI_ops["Move"] = trap_width / move_velocity; // us
        TI_ops["Junction"] = trap_width / junction_move_velocity;

        // Includes Merge, Cool, Interact, and Split 
        TI_ops["ZZ"] = 2000;

        // Initialize map to Stim instructions while we're at it 
        TI_ops_to_stim["Prepare_Z"] = "RZ";
        TI_ops_to_stim["Measure_Z"] = "MZ";
        TI_ops_to_stim["X_pi/2"] = "X";
        TI_ops_to_stim["Y_pi/2"] = "Y";
        TI_ops_to_stim["Z_pi/2"] = "Z";
        TI_ops_to_stim["X_pi/4"] = "SQRT_X";
        TI_ops_to_stim["Y_pi/4"] = "SQRT_Y";
        TI_ops_to_stim["Z_pi/4"] = "SQRT_Z";
        TI_ops_to_stim["X_-pi/4"] = "SQRT_X_DAG";
        TI_ops_to_stim["Y_-pi/4"] = "SQRT_Y_DAG";
        TI_ops_to_stim["Z_-pi/4"] = "SQRT_Z_DAG";
        TI_ops_to_stim["ZZ"] = "SQRT_ZZ";

    }

    // Helper function to add the Prepare_Z HW_Instruction to a circuit given a plaquette and a qubit label
    double HardwareModel::add_init(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const {

        // If the plaquette in question does not have the present qubit label, do nothing
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        if (p.get_qsite(qubit) == uint_max) {
            return time;
        }

        // Perform validity check
        if ((*p.grid())[p.get_qsite(qubit)] != 'O') {
            std::cerr << "HardwareModel::add_init: Can only Prepare Z at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instruction onto the circuit
        circuit.push_back(HW_Instruction("Prepare_Z", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_operator_type()));

        // Return updated time
        return time + TI_ops.at("Prepare_Z");

    }

    // Helper function to add the Measure_Z HW_Instruction to a circuit given a plaquette and a qubit label
    double HardwareModel::add_meas(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const {

        // If the plaquette in question does not have the present qubit label, do nothing
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        if (p.get_qsite(qubit) == uint_max) {
            return time;
        }

        // Perform validity check
        if ((*p.grid())[p.get_qsite(qubit)] != 'O') {
            std::cerr << "HardwareModel::add_meas: Can only apply Measure Z at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        circuit.push_back(HW_Instruction("Measure_Z", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_operator_type()));

        // Return updated time
        return time + TI_ops.at("Measure_Z");

    }

    // Helper function to add H gate in terms of native TI gates to a circuit given a plaquette and a qubit label
    double HardwareModel::add_H(const Plaquette& p, char qubit, double time, unsigned int step, std::vector<HW_Instruction>& circuit) const {
        
        // If the plaquette in question does not have the present qubit label, do nothing
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        if (p.get_qsite(qubit) == uint_max) {
            return time;
        }

        // Perform validity check
        if ((*p.grid())[p.get_qsite(qubit)] != 'O') {
            std::cerr << "HardwareModel::add_H: Can only apply Hadamard gates at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        circuit.push_back(HW_Instruction("Y_-pi/4", p.get_qsite(qubit), uint_max, time, step, qubit, ' ', p.get_shape(), p.get_operator_type()));
        circuit.push_back(HW_Instruction("Z_pi/2", p.get_qsite(qubit), uint_max, time + TI_ops.at("Y_-pi/4"), step, qubit, ' ', p.get_shape(), p.get_operator_type()));

        // Return updated time
        return time + TI_ops.at("Y_-pi/4") + TI_ops.at("Z_pi/2");

    }

    // Helper function to add the Prepare_Z HW_Instruction to a circuit given a qsite
    double HardwareModel::add_init(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_init: Can only Prepare Z at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instruction onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Prepare_Z", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Prepare_Z");

    }

    double HardwareModel::add_X(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_X: Can only apply X gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("X_pi/2", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("X_pi/2");

    }

    double HardwareModel::add_X_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_X_pi4: Can only apply sqrt(X) gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("X_pi/4", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("X_pi/4");

    }

    double HardwareModel::add_X_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_X_mpi4: Can only apply sqrt(X) gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("X_-pi/4", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("X_-pi/4");

    }

    double HardwareModel::add_Y(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_Y: Can only apply Y gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Y_pi/2", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Y_pi/2");

    }

    double HardwareModel::add_Y_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_Y_pi4: Can only apply sqrt(Y) gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Y_pi/4", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Y_pi/4");

    }

    double HardwareModel::add_Y_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_Y_mpi4: Can only apply sqrt(Y) gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Y_-pi/4", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Y_-pi/4");

    }

    double HardwareModel::add_Z(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_Z: Can only apply Z gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Z_pi/2", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Z_pi/2");

    }

    double HardwareModel::add_Z_pi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_Z_pi4: Can only apply sqrt(Z) gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Z_pi/4", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Z_pi/4");

    }

    double HardwareModel::add_Z_mpi4(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_Z_mpi4: Can only apply sqrt(Z) gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Z_-pi/4", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Z_-pi/4");

    }

    double HardwareModel::add_Z_pi8(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_Z_pi8: Can only apply sqrt(Z) gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Z_pi/8", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Z_pi/8");

    }

    double HardwareModel::add_Z_mpi8(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {

        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_Z_mpi8: Can only apply sqrt(Z) gate at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Z_-pi/8", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Z_-pi/8");

    }

    // Helper function to add H gate in terms of native TI gates to a circuit given a qsite
    double HardwareModel::add_H(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {
        
        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_H: Can only apply Hadamard gates at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instructions onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Y_-pi/4", site, uint_max, time, step, 'X', ' ', 'X', 'X'));
        circuit.push_back(HW_Instruction("Z_pi/2", site, uint_max, time + TI_ops.at("Y_-pi/4"), step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Y_-pi/4") + TI_ops.at("Z_pi/2");

    }

    double HardwareModel::add_meas(unsigned int site, double time, unsigned int step, const GridManager& grid, std::vector<HW_Instruction>& circuit) const {
        // Perform validity check
        if (grid[site] != 'O') {
            std::cerr << "HardwareModel::add_meas: Can only apply Measure Z at 'O' QSites." << std::endl;
            abort();
        }

        // Push corresponding HW_Instruction onto the circuit
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        circuit.push_back(HW_Instruction("Measure_Z", site, uint_max, time, step, 'X', ' ', 'X', 'X'));

        // Return updated time
        return time + TI_ops.at("Measure_Z");
    }

    // Helper function to add CNOT gate in terms of native TI gates to a circuit given a plaquette and qubit labels
    double HardwareModel::add_CNOT(Plaquette& p, char control, char target, double time, unsigned int step, const GridManager& grid, 
        std::vector<HW_Instruction>& circuit) const {

        // If the plaquette in question does not have the present qubit labels, do nothing
        unsigned int uint_max = std::numeric_limits<unsigned int>::max();
        if ((p.get_qsite(control) == uint_max) || (p.get_qsite(target) == uint_max)) {
            return time;
        }
        
        // This function currently requires one of the qubits to be 'm'
        if (control != 'm' && target != 'm') {
            std::cerr << "HardwareModel::add_CNOT: Currently requires either the control or the target to be 'm'." << std::endl;
            abort();
        }

        // This function currently depends on both control and target beginning at their 'home' sites on the plaquette
        if ((!p.is_home(control)) || (!p.is_home(target))) {
            std::cerr << "HardwareModel::add_CNOT: Currently requires control and target to be at their home site on the plaquette." << std::endl;
            abort();
        }

        // A qubit cannot perform a CNOT gate with itself
        if (control == target) {
            std::cerr << "HardwareModel::add_CNOT: A qubit cannot perform a CNOT gate with itself." << std::endl;
            abort();
        }

        // Obtain data qubit label
        char data_qubit = control;
        if (control == 'm') {
            data_qubit = target;
        }

        // Perform initial rotation
        circuit.push_back(HW_Instruction("Y_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_operator_type()));
        time += TI_ops.at("Y_-pi/4");

        // Retrieve path to data qubit and apply move_along_path
        std::vector<unsigned int> path = grid.get_path(p.get_qsite('m'), p.get_qsite(data_qubit));
        move_along_path(p, step, circuit, time, path, grid);

        // Apply ZZ operation
        circuit.push_back(HW_Instruction("ZZ", p.get_qsite(data_qubit), p.get_qsite('m'), time, step, data_qubit, 'm', p.get_shape(), p.get_operator_type()));
        time += TI_ops.at("ZZ");

        // Move the measure qubit back to its home base for further operations
        reverse(path.begin(), path.end());
        move_along_path(p, step, circuit, time, path, grid);

        // Before final rotations, ensure control and target are home
        assert(p.is_home(control) && p.is_home(target)); 

        // Perform final rotations
        circuit.push_back(HW_Instruction("Z_-pi/4", p.get_qsite(control), uint_max, time, step, control, ' ', p.get_shape(), p.get_operator_type()));
        circuit.push_back(HW_Instruction("X_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_operator_type()));
        time += TI_ops.at("X_-pi/4");
        circuit.push_back(HW_Instruction("Z_-pi/4", p.get_qsite(target), uint_max, time, step, target, ' ', p.get_shape(), p.get_operator_type()));

        // Return updated time
        return time + TI_ops.at("Z_-pi/4");
    }

    // Translate all qubits left one column on the grid
    double HardwareModel::shift_left(const std::set<unsigned int>& qsites, unsigned int patch_cols, GridManager& grid, std::vector<HW_Instruction>& hw_master, double time) const {

        // For use in time accounting
        double time_tmp;

        // Move qubits out of the way (south-west next to measure qubit) in order to shift the left-wards strip of data qubits over
        for (unsigned int q : qsites) {

            // Check that it is occupied
            if (!grid.is_occupied(q)) {
                std::cerr << "HardwareModel::shift_left: Movement ordered on unoccupied qsite." << std::endl; abort();
            }

            // Get path next to measure qubit
            unsigned int row = grid.get_row(q);
            unsigned int col = grid.get_col(q);
            std::vector<unsigned int> path = grid.get_path(grid.index_from_coords(row, col, 1), q);
            path.push_back(q);
            reverse(path.begin(), path.end());
            path.pop_back();

            // Move along path
            // **Note: This changes the set of occupied sites on the grid
            time_tmp = move_along_path(q, 0, hw_master, time, path, grid);

        }

        // std::vector<std::string> ascii_grid = grid.ascii_grid(true);
        // grid.print_grid(ascii_grid);

        // Update time
        time = time_tmp;

        /* Shuttle the occupied sites (if they are indeed occupied) through to the right */

        // First, get the strip to the left 
        std::set<unsigned int> strip;
        for (unsigned int q : qsites) {
            if (qsites.find(grid.shift_qsite(q, 0, -1)) == qsites.end())
                strip.insert(grid.shift_qsite(q, 0, -1));
        }

        for (unsigned int q : strip) {

            time_tmp = time;

            // Check that it is occupied
            if (!grid.is_occupied(q)) continue;

            for (unsigned int i = 0; i<patch_cols; i++) {

                // Get path next to measure qubit
                std::vector<unsigned int> path = grid.get_path(grid.shift_qsite(q, 0, i), grid.shift_qsite(q, 0, i+1));
                path.push_back(grid.shift_qsite(q, 0, i+1));

                // Move along path
                // **Note: This changes the set of occupied sites on the grid
                time_tmp = move_along_path(grid.shift_qsite(q, 0, i), 0, hw_master, time_tmp, path, grid);

            }

        }

        // ascii_grid = grid.ascii_grid(true);
        // grid.print_grid(ascii_grid);

        // Update time
        time = time_tmp;

        // Then move our data qubits north-west to the position we want
        for (unsigned int q : qsites) {

            // Get path to final site
            unsigned int row = grid.get_row(q);
            unsigned int col = grid.get_col(q);
            std::vector<unsigned int> path = grid.get_path(grid.index_from_coords(row, col, 1), grid.shift_qsite(q, 0, -1));
            path.push_back(grid.shift_qsite(q, 0, -1));
            reverse(path.begin(), path.end());
            path.pop_back();
            reverse(path.begin(), path.end());

            // Move along path
            // **Note: This changes the set of occupied sites on the grid
            time_tmp = move_along_path(path[0], 0, hw_master, time, path, grid);
        }

        // ascii_grid = grid.ascii_grid(true);
        // grid.print_grid(ascii_grid);

        /* Before any subsequent idle operation, will need to make sure measure qubits are there */
        return time_tmp;
    }

    // Move the measure qubit of a given plaquette along a path (to the closest site adjacent to a data qubit)
    void HardwareModel::move_along_path(Plaquette& p, unsigned int step, std::vector<HW_Instruction>& circuit, double& time,
        const std::vector<unsigned int>& path, const GridManager& grid) const {

        // Validity check 
        if (p.get_qsite('m') != path[0]) {
            std::cerr << "HardwareModel::move_along_path: current site inconsistent with path given." << std::endl;
            abort();
        }

        // Construct HW_Instructions to Move 'm' along the path while also applying the move_to_site plaquette member function (this ensures validity on the grid)
        // We note an assumption is that there are never two Junctions in a row on the grid
        unsigned int through_J = 0;
        for (unsigned int i=1; i<path.size(); i++) {
            if (grid[path[i]] == 'J') {
                through_J = 1;
                continue;
            }
            circuit.push_back(HW_Instruction("Move", path[i-1-through_J], path[i], time, step, 'm', ' ', p.get_shape(), p.get_operator_type()));
            time += TI_ops.at("Move") + through_J*TI_ops.at("Junction");
            if (through_J == 1) {through_J = 0;}
            p.move_to_site('m', path[i]);
        }

    }

    // Move a qubit starting at qsite along a path
    double HardwareModel::move_along_path(unsigned int qsite, unsigned int step, std::vector<HW_Instruction>& circuit, double time,
        const std::vector<unsigned int>& path, GridManager& grid) const {

        // Validity check 
        if (qsite != path[0]) {
            std::cerr << "HardwareModel::move_along_path: current site inconsistent with path given." << std::endl;
            std::cerr << qsite << " " << path[0] << std::endl;
            abort();
        }

        // Construct HW_Instructions to Move along the path 
        // We note an assumption is that there are never two Junctions in a row on the grid
        unsigned int through_J = 0;
        for (unsigned int i=1; i<path.size(); i++) {
            if (grid[path[i]] == 'J') {
                through_J = 1;
                continue;
            }
            circuit.push_back(HW_Instruction("Move", path[i-1-through_J], path[i], time, step, 'X', ' ', 'X', 'X'));
            time += TI_ops.at("Move") + through_J*TI_ops.at("Junction");
            grid.move_qubit(path[i-1-through_J], path[i]);
            if (through_J == 1) {through_J = 0;}
        }

        return time;

    }
    HardwareModel::HardwareModel() : trap_width(420.0) {
        cell_width = 4*trap_width;
        init_TI_ops();
    }

    void HardwareModel::print_TI_ops() const {
        for (auto const &op_pair : TI_ops) {
            std::cout << op_pair.first << ": " << op_pair.second << std::endl;
        }
    }
}