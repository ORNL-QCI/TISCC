#ifndef TISCC_INSTRUCTION_HPP
#define TISCC_INSTRUCTION_HPP

#include<string>

namespace TISCC {

// The Instruction class represents a named instruction which acts on the qubit labels (a, b, c, d, or m) of a surface code plaquette
class Instruction {
public:
    explicit Instruction(std::string name, char q1, char q2) : name_(name), q1_(q1), q2_(q2) {};
    const std::string& get_name() const {return name_;}
    const char get_q1() const {return q1_;}
    const char get_q2() const {return q2_;}
private:
    std::string name_;
    char q1_;
    char q2_;
};

// The HW_Instruction class represents a named instruction acting at a time at particular site(s) on the grid
class HW_Instruction {
public:
    // We note that the latter four of these instructions are not strictly necessary but are stored for debugging output
    explicit HW_Instruction(std::string name, unsigned int site1, unsigned int site2, float time, unsigned int step, char q1, char q2, char shape, char type) : 
        name_(name), site1_(site1), site2_(site2), time_(time), step_(step), q1_(q1), q2_(q2), shape_(shape), type_(type) {};
    const std::string& get_name() const {return name_;}
    const unsigned int get_site1() const {return site1_;}
    const unsigned int get_site2() const {return site2_;}
    const float get_time() const {return time_;}
    const unsigned int get_step() const {return step_;}

    // These accessors and associated member variables are only being included now for the purpose of debugging
    const char get_q1() const {return q1_;}
    const char get_q2() const {return q2_;}
    const char get_shape() const {return shape_;}
    const char get_type() const {return type_;}
private:
    std::string name_;
    unsigned int site1_;
    unsigned int site2_;
    float time_;
    unsigned int step_;
    char q1_;
    char q2_;
    char shape_;
    char type_;
};

// Comparison operator for use in sorting hardware instructions before printing them
bool operator<(const HW_Instruction& i1, const HW_Instruction& i2);

}

#endif //TISCC_INSTRUCTION_HPP