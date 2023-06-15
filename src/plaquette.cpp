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

    // Ask whether the plaquette has support on the site in question and remove if so
    bool Plaquette::remove_supported_qsite(unsigned int site) {
        if ((a_ == site) || (b_ == site) || (c_ == site) || (d_ == site)) {
            return true;
        }
        else if (m_ == site) {
            std::cerr << "Plaquette::remove_supported_qsite: Cannot remove measure qubit of active plaquette." << std::endl;
            abort();
        }
        else {return false;}
    }

    // Move qubit to a specified site after checking validity and recording change on grid
    void Plaquette::move_to_site(char q, unsigned int site) {

        // Update the site of the qubit and update the set of occupied sites in the grid
        grid_->move_qubit(get_qsite(q), site);
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

    // Check whether a qubit is currently at its home location on the grid
    bool Plaquette::is_home(char qubit) const {
        if (qubit == 'm') {
            return m_ == grid_->index_from_coords(row_, col_, 1);
        }
        else if (qubit == 'a') {
            return a_ == grid_->index_from_coords(row_, col_-1, 5);
        }
        else if (qubit == 'b') {
            return b_ == grid_->index_from_coords(row_, col_, 5);
        }
        else if (qubit == 'c') {
            return c_ == grid_->index_from_coords(row_+1, col_-1, 5);
        }
        else if (qubit == 'd') {
            return d_ == grid_->index_from_coords(row_+1, col_, 5);
        }
        else {
            std::cerr << "Plaquette::is_home: Invalid qubit character given." << std::endl;
            abort();
        }
    }

   // Check to see if a given instruction is valid on this plaquette 
    bool Plaquette::is_instr_valid(const Instruction& instr) const {
        bool n_valid = !(shape_ == 'n' && ((instr.get_q1() == 'a') || (instr.get_q1() == 'b') || (instr.get_q2() == 'a') || (instr.get_q2() == 'b')));
        bool s_valid = !(shape_ == 's' && ((instr.get_q1() == 'c') || (instr.get_q1() == 'd') || (instr.get_q2() == 'c') || (instr.get_q2() == 'd')));
        bool e_valid = !(shape_ == 'e' && ((instr.get_q1() == 'b') || (instr.get_q1() == 'd') || (instr.get_q2() == 'b') || (instr.get_q2() == 'd')));
        bool w_valid = !(shape_ == 'w' && ((instr.get_q1() == 'a') || (instr.get_q1() == 'c') || (instr.get_q2() == 'a') || (instr.get_q2() == 'c')));
        bool return_val = n_valid && s_valid && e_valid && w_valid;
        return return_val;
    }
}