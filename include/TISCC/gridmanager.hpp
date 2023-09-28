#ifndef TISCC_GRIDMANAGER_HPP
#define TISCC_GRIDMANAGER_HPP

#include <TISCC/plaquette.hpp>
#include <TISCC/instruction.hpp>
#include <TISCC/hardwaremodel.hpp>

#include <iostream>
#include <vector>
#include <set>

namespace TISCC {

// Forward declaration of dependencies
class Plaquette;
class HW_Instruction;

/**
 * @brief Provides access to an array representing a grid-like trapped-ion hardware architecture with a particular repeating pattern of trapping zones.
 * 
 * The repeating unit is {`M', `O', `M', `J', `M', `O', `M'} (see SiteType for meaning of characters).
 * It has two straight segments, one pointed down-ward and one pointed right-ward, connected by a junction. 
 */
class GridManager {
public:
    /**
     * @brief Represents different types of qsites in the hardware architecture.
     */
    enum SiteType : char {
        QSite_Memory = 'M', /**< Memory site ('M') */
        QSite_Memory_and_Ops = 'O', /**< Operation site ('O') */
        Junction = 'J', /**< Junction site ('J') */
        Nothing = ' ' /**< Non-site (not used by default) */
    };
    
    /**
     * @brief Constructor for the GridManager class.
     * @param nrows Number of rows in the hardware grid.
     * @param ncols Number of columns in the hardware grid.
     */
    explicit GridManager(unsigned int nrows, unsigned int ncols);

    /**
     * @brief Destructor for the GridManager class.
     */
    ~GridManager() {delete[] grid_;}

    /**
     * @brief Provides read-only access to the array representing the hardware grid.
     * @param i The grid array index (representing a qsite).
     * @return A const reference to the SiteType at the specified qsite.
     */ 
    const SiteType& operator[] (unsigned int i) const {return grid_[i];}

    /**
     * @brief Provides read-only access to the hardware grid based on coordinates.
     * @param row The row coordinate.
     * @param col The column coordinate.
     * @param idx The index within repeating unit of hardware architecture.
     * @return A const reference to the SiteType at the specified coordinates.
     */
    const SiteType& from_coords(unsigned int row, unsigned int col, unsigned int idx) const; 

    /**
     * @brief Provides grid array index (representing a qsite) based on hardware grid coordinates.
     * @param row The row coordinate.
     * @param col The column coordinate.
     * @param idx The index within repeating hardware unit.
     * @return The grid array index (representing a qsite) at the specified hardware grid coordinates.
     */
    unsigned int index_from_coords(unsigned int row, unsigned int col, unsigned int idx) const;

    /**
     * @brief Provides index within repeating hardware unit from grid array index (representing a qsite).
     * @param i The grid array index (representing a qsite).
     * @return The index within repeating hardware unit.
     */
    unsigned int get_idx(unsigned int i) const {return i%7;}

    /**
     * @brief Provides column coordinate for qsite from grid array index.
     * @param i The grid array index (representing a qsite).
     * @return The column coordinate.
     */
    unsigned int get_col(unsigned int i) const {return ((i-i%7)/7)%ncols_;}

    /**
     * @brief Provides row coordinate for qsite from grid array index.
     * @param i The grid array index (representing a qsite).
     * @return The row coordinate.
     */
    unsigned int get_row(unsigned int i) const {return (((i-i%7)/7)-get_col(i))/ncols_;}

    /**
     * @brief Provides shifted grid array index corresponding to the translation of a qsite on the hardware grid.
     * @param nrows Number of rows to shift east-ward.
     * @param ncols Number of columns to shift south-ward
     * @return The shifted grid array index (representing a qsite).
     */
    unsigned int shift_qsite(unsigned int i, int nrows, int ncols) const;

    /**
     * @brief Returns the number of rows in the hardware grid.
     * @return The number of rows in the hardware grid.
     */
    unsigned int get_nrows() const {return nrows_;}

    /**
     * @brief Returns the number of columns in the hardware grid.
     * @return The number of columns in the hardware grid.
     */
    unsigned int get_ncols() const {return ncols_;}

    /**
     * @brief Returns a set containing the grid array indices (qsites) that are occupied (by default, all 'O' sites).
     * @return A const reference to the set containing the grid array indices (qsites) that are occupied (by default, all 'O' sites).
     */
    const std::set<unsigned int>& get_occ_sites() const {return occupied_sites;}

    /**
     * @brief Provides a path from one 'O' site to the closest site adjacent to a surrounding 'O' site.
     * 
     * Depends on the details of the repeating unit being used by GridManager.
     * @param site1 The first 'O' site.
     * @param site2 The second 'O' site.
     * @return A vector representing the path between the two sites (specified by grid array indices).
     */
    std::vector<unsigned int> get_path(unsigned int site1, unsigned int site2) const;

    /**
     * @brief Returns all sites adjacent to a given qsite depending on its index within the repeating unit.
     * 
     * Depends on the details of the repeating unit being used by GridManager.
     * @return A set containing all sites adjacent to a given qsite depending on its index within the repeating unit.
     */
    std::set<unsigned int> get_adjacent(unsigned int site) const;

    /**
     * @brief Checks occupation of a particular qsite.
     * @param site The qsite being queried.
     * @return A boolean representing the occupation of site.
    */
    bool is_occupied(unsigned int site) const {return (occupied_sites.find(site) != occupied_sites.end());}

    /**
     * @brief Flips occupation state of two sites i.e. "move a qubit'' (relies on there only ever being one qubit per site).
     * 
     * Ensures that the target site is either (a) adjacent to or (b) adjacent to an adjacent junction (and if so, tracks that junction).
     * Checks validity of movement based on the current state of occupied_sites.
     * @param site1 The site occupied by the qubit to be moved.
     * @param site2 The site on which to move the qubit.
     * @return A qsite corresponding to the junction passed through (if applicable).
    */
    unsigned int move_qubit(unsigned int site1, unsigned int site2);

    /** 
     * @brief Provides a plaquette object "pinned" at a given grid point.
     * 
     * Depends on the details of the repeating unit being used by GridManager and the hardware mapping of a surface code plaquette.
     * @param row The row coordinate of the grid point.
     * @param col The column coordinate of the grid point.
     * @param shape Defines either a four-qubit plaquette ('f') or a two-qubit plaquette with edge type {'n','s','e','w'}.
     * @param type Defines whether it is an 'X' or a 'Z' plaquette.
     * @return A Plaquette object "pinned" at a given grid point
    */
    Plaquette get_plaquette(unsigned int row, unsigned int col, char shape, char type);

    /**
     * @brief De-occupies a qsite on the grid
     * @param site The qsite being de-occupied.
    */
    void deoccupy_site(unsigned int site) {occupied_sites.erase(site);}

    /**
     * @brief Enforces the validity of hardware instructions by simulating qubit movements on the hardware grid.
     * @param hw_master A vector of hardware instructions representing the master circuit.
    */
    void enforce_hw_master_validity(std::vector<HW_Instruction>& hw_master);

    /**
     * @brief Prints resource counts for a given hardware circuit using the HardwareModel.
     * 
     * Outputs grid area, computation time, space-time volume, trapping zone-seconds, and active trapping zone-seconds.
     * @param hw_master A vector of hardware instructions representing the master circuit.
    */
    void resource_counter(const std::vector<HW_Instruction>& hw_master) const;

    /**
     * @brief Prints a table of qsites corresponding to hardware coordinates (row, col, repeating unit idx).
    */
    void print_qsite_mapping() const;

    /**
     * @brief Prints the sites occupied on the grid.
    */
    void print_occ_sites() const;

    /**
     * @brief Constructs ASCII visualization of the hardware grid.
     *  
     * Depends on the details of the repeating unit being used by GridManager.
     * @param occ_mode If true, show SiteType labels for only occupied sites.
     * @return A vector of strings that can be readily printed using print_grid.
    */
    std::vector<std::string> ascii_grid(bool occ_mode) const;

    /**
     * @brief Constructs ASCII visualization of the hardware grid including a Pauli operator to visualize.
     * 
     * Depends on the details of the repeating unit being used by GridManager.
     * @param occ_mode If true, show SiteType labels for only occupied sites.
     * @param qsites A vector of pairs corresponding with qsites and corresponding Pauli operators.
     * @return A vector of strings that can be readily printed using print_grid.
    */
    std::vector<std::string> ascii_grid_with_operator(const std::vector<std::pair<unsigned int,char>>& qsites, bool occ_mode) const;

    /**
     * @brief Print grid using the output of ascii_grid or ascii_grid_with_operator.
     * @param ascii_grid Output of ascii_grid or ascii_grid_with_operator.
    */
    void print_grid(std::vector<std::string>& ascii_grid) const;

private:
    SiteType* grid_;
    unsigned int nrows_;
    unsigned int ncols_;
    std::set<unsigned int> occupied_sites;
};
}

#endif //TISCC_GRIDMANAGER_HPP