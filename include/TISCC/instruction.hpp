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

}
#endif //TISCC_INSTRUCTION_HPP