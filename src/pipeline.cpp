#include <TISCC/pipeline.hpp>

#include <argparse/argparse.h>
#include <iostream>

namespace TISCC 
{
    // Define possible types of sites in the trapped ion micro layout
    enum SiteType : char {
        QSite_Memory = 'M',
        QSite_Memory_and_Ops = 'O',
        Junction = 'J',
        Nothing = 'X'
    };

    // Parses command line arguments and runs through the code pipeline
    int run_tiscc(
        int argc, const char* argv[],
        std::istream& in_stream,
        std::ostream& out_stream,
        std::ostream& err_stream)
    {
        // Setting up command line argument parser
        std::string prog_name{argv[0]};
        argparse::ArgumentParser parser(prog_name, "Trapped-Ion Surface Code Compiler");
        parser.add_argument()
                .names({"-d", "--dummy"})
                .description("Dummy argument to test argparse")
                .required(false);
        parser.enable_help();
        auto err = parser.parse(argc, argv);
        if (err)
        {
            err_stream << err << std::endl;
            parser.print_help();
            return -1;
        }
        if (parser.exists("help"))
        {
            parser.print_help();
            return 0;
        }

        // Use of dummy command line flag
        if (parser.exists("d"))
        {
            std::cout << "Hello world!" << std::endl;
        }

        // Initializing grid on the heap
        // TODO: Create an object called GridManager that has the root pointer as a private member variable
        //  - Construct the grid in a private member function with repeating unit specified in the constructor
        //  - Q: if I have the root as a private member variable and expose it through a const public member function, will the whole array be treated as const?
        unsigned int nrows = 3;
        unsigned int ncols = 3;
        SiteType*** a = new SiteType**[nrows]; 
        for (unsigned int i=0; i<nrows; i++) {
            a[i] = new SiteType*[ncols];
            // The third level of the array contains the repeating unit of the grid
            for (unsigned int j=0; j<ncols; j++) {
                a[i][j] = new SiteType[7];
                a[i][j][0] = SiteType::QSite_Memory;
                a[i][j][1] = SiteType::QSite_Memory_and_Ops;
                a[i][j][2] = SiteType::QSite_Memory;
                a[i][j][3] = SiteType::Junction;
                a[i][j][4] = SiteType::QSite_Memory;
                a[i][j][5] = SiteType::QSite_Memory_and_Ops;
                a[i][j][6] = SiteType::QSite_Memory;
            }
        }

        // Print out grid
        std::cout << a << std::endl;
        for (unsigned int i=0; i<nrows; i++) {
            std::cout << a[i] << std::endl;
            for (unsigned int j=0; j<ncols; j++) {
                std::cout << a[i][j] << std::endl;
                SiteType* k;
                for (SiteType* k=a[i][j]; k<a[i][j]+7; k++) {
                    std::cout << k << std::endl;
                }
            }
        }

        // Create a hash table
        //std::unordered_map<SiteType*, int> qsite_hash

        return 0;
    }
}