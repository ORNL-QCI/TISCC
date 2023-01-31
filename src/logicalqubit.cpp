#include <TISCC/logicalqubit.hpp>

namespace TISCC 
{
    LogicalQubit::LogicalQubit(unsigned int dx, unsigned int dz, const GridManager& a) {
        // TODO: Construct the plaquettes of the surface code

        // Keep track of all data and measure qubits separately so that we can perform collision checks with other LQ

                // // Create an array of data qubits and measure qubits, which are really just pointers to sites on the lattice
        // unsigned int num_data_qubits = dx*dz;
        // unsigned int num_measure_qubits = num_data_qubits-1;
        // SiteType* data_qubits[num_data_qubits];
        // SiteType* measure_qubits[num_measure_qubits];
        // for (unsigned int i=0; i<nrows; i++) {
        //     for (unsigned int j=0; j<ncols; j++) {

        //     }
        // }
    }

    void LogicalQubit::idle(unsigned int cycles) {
        // TODO: Perform operations and print them
    }

}