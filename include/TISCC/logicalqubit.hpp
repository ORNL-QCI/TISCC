#ifndef TISCC_LOGICALQUBIT_HPP
#define TISCC_LOGICALQUBIT_HPP

#include <TISCC/gridmanager.hpp>

#include <vector>

namespace TISCC {

struct Plaquette {
    std::vector<GridManager::SiteType*> data_qubits;
    std::vector<GridManager::SiteType*> measure_qubits;
};

class LogicalQubit {
public:
    explicit LogicalQubit(unsigned int dx, unsigned int dz, const GridManager& a);
    void idle(unsigned int cycles);

private:
    std::vector<Plaquette> x_plaquettes;
    std::vector<Plaquette> z_plaquettes;
    std::vector<GridManager::SiteType*> data_qubits;
    std::vector<GridManager::SiteType*> measure_qubits;
};
}

#endif //TISCC_LOGICALQUBIT_HPP