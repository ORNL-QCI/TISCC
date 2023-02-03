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

    // Private member function that returns an lvalue allowing you to modify the qsite of a qubit
    unsigned int& Plaquette::mod_qsite(char qubit) {
        if (qubit == 'a') {return a_;}
        else if (qubit =='b') {return b_;}
        else if (qubit =='c') {return c_;}
        else if (qubit =='d') {return d_;}
        else if (qubit =='m') {return m_;}
        else {std::cerr << "Plaquette::move_home: Invalid character given." << qubit << std::endl; abort();}
    }

    // Move named qubit to its 'home' qsite on the grid
    // TODO: This function depends on the specifics of the grid and its functionality should be moved into GridManager in the future
    void Plaquette::move_home(char qubit) {
        if (qubit == 'm') {
            mod_qsite(qubit) = grid_.index_from_coords(row_, col_, 1);
        }
        else {
            std::cerr << "Plaquette::move_home: Not yet implemented for non-'m' qubits." << std::endl; 
            abort();
        }
    }  
}