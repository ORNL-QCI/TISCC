#include <TISCC/plaquette.hpp>

namespace TISCC 
{
    // Access qsite occupied by a named qubit
    unsigned int Plaquette::get_qsite(char qubit) const {
        if (qubit == 'a') {return a_;}
        else if (qubit =='b') {return b_;}
        else if (qubit =='c') {return c_;}
        else if (qubit =='d') {return d_;}
        else if (qubit =='m') {return m_;}
        else {std::cerr << "Plaquette::get_qsite: Invalid character given." << qubit << std::endl; abort();}
    }

    // Move qubit 1 to location of qubit 2 (deprecated since we no longer allow qubits to occupy the same site)
    // void Plaquette::apply_move(char q1, char q2) {
    //     mod_qsite(q1) = get_qsite(q2);
    // }

    // Move qubit to a specified site after checking validity and recording change on grid
    void Plaquette::move_to_site(char q, unsigned int site) {

        /* TODO: A qubit cannot sit at a junction
            - Will need to update HardwareModel to skip through junctions
            - Will need to have twice the time for instructions that skip through junctions */
        // if (grid_[site] == 'J') {
        //     std::cerr << "Plaquette::move_to_site: attempted move to a junction." << std::endl;
        //     abort();
        // }

        // For a move to be valid, it must target an unoccupied site (this is covered by grid_.move_qubit())
        // if (grid_.is_occupied(site)) {
        //     std::cerr << q << " " << site << std::endl;
        //     std::cerr << "Plaquette::move_to_site: attempted move to occupied site." << std::endl;
        //     abort();            
        // }
        
        /* TODO: For a move to be valid, it must target an adjacent site (skipping junctions)
            - Must write a function to return the set of allowed moves from a given site */
        // std::set<unsigned int> adjacent = grid_.get_adjacent(get_qsite(q));
        // if (adjacent.find(site)==adjacent.end()) {
        //     std::cerr << "Plaquette::move_to_site: attempted move to non-adjacent site." << std::endl;
        //     abort();
        // }

        // Update the site of the qubit and update the set of occupied sites in the grid
        grid_.move_qubit(get_qsite(q), site);
        mod_qsite(q) = site;
        
    }

    // Private member function that returns an lvalue allowing one to modify the qsite of a qubit
    unsigned int& Plaquette::mod_qsite(char qubit) {
        if (qubit == 'a') {return a_;}
        else if (qubit =='b') {return b_;}
        else if (qubit =='c') {return c_;}
        else if (qubit =='d') {return d_;}
        else if (qubit =='m') {return m_;}
        else {std::cerr << "Plaquette::mod_qsite: Invalid character given." << qubit << std::endl; abort();}
    }

    // Move named qubit to its 'home' qsite on the grid (deprecated since we no longer allow qubits move arbitrary distances)
    // Note: If this functionality is needed in the future, it should be moved into GridManager
    // void Plaquette::move_home(char qubit) {
    //     if (qubit == 'm') {
    //         mod_qsite(qubit) = grid_.index_from_coords(row_, col_, 1);
    //     }
    //     else {
    //         std::cerr << "Plaquette::move_home: Not yet implemented for non-'m' qubits." << std::endl; 
    //         abort();
    //     }
    // }  
}