#ifndef TISCC_INSTRUCTION_HPP
#define TISCC_INSTRUCTION_HPP

#include<string>

namespace TISCC {

class Instruction {
public:
    explicit Instruction(std::string name, char q1, char q2, float time) : name_(name), q1_(q1), q2_(q2), time_(time) {};
    const std::string& get_name() const {return name_;}
    const char get_q1() const {return q1_;}
    const char get_q2() const {return q2_;}
    const float get_time() const {return time_;}
private:
    std::string name_;
    char q1_;
    char q2_;
    float time_; 
};

class HW_Instruction {
public:
    explicit HW_Instruction(std::string name, unsigned int site1, unsigned int site2, float time, unsigned int step) : 
        name_(name), site1_(site1), site2_(site2), time_(time), step_(step) {};
    const std::string& get_name() const {return name_;}
    const unsigned int get_site1() const {return site1_;}
    const unsigned int get_site2() const {return site2_;}
    const float get_time() const {return time_;}
    const unsigned int get_step() const {return step_;}
private:
    std::string name_;
    unsigned int site1_;
    unsigned int site2_;
    float time_;
    unsigned int step_;
};

// Comparison operator for use in sorting hardware instructions before printing them
bool operator<(const HW_Instruction& i1, const HW_Instruction& i2);
}
#endif //TISCC_INSTRUCTION_HPP