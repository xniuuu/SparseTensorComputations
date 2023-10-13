#include "src/helpers.h"
#include <ctype.h>
#include <thread> 

using namespace std;

int main(int argc, char *argv[]){
    
    // Parse command line arguments
    if (argc < 8) {
        cerr << "Error: Not enough arguments provided \n";
        std::cerr << "Usage: " << argv[0] << "<TENSOR> <OFFSET> <SEP> <TENSOR_DIM> <BLOCK_SIZE> <MODE> <TRAVERSAL> [<VECTOR>]\n";
        cerr << " Please check the make file, if all neccessary parameters were provided\n";
        return 1;
    }

    string tns_file = argv[1];
    bool offset = (argv[2][0] != '0') ? 1 : 0;
    char sep = argv[3][0];

    vector<unsigned int> n_d;
    stringstream entry(argv[4]);
    string elem;
    while (getline(entry, elem, ',')) {
        n_d.push_back(stoi(elem));
    }
   
    int block_size = stoi(argv[5]);
    int mode = stoi(argv[6]);
    string traversal_curve = argv[7];
    
    // Add a one-second delay
    this_thread::sleep_for(chrono::seconds(1));
    
    BlockedCOO block_coo(tns_file, offset, sep, n_d, block_size, mode, traversal_curve);
 
    // Construct BlockedBitIf object
    BlockedBitIf b_bitif(block_coo);

    //Compute TVM if vector file is provided
    if (argc >= 9) {
        vector<double> vec;
        string vec_file = argv[8];
        //Load vector from file
        vec = getVecFromFile(vec_file);
        b_bitif.computeTVMBlockedRuntime(vec);

        //block_coo.computeTVM(vec);
        //vector<vector<double> > valTVM = block_coo.getValTVM();
        //b_bitif.validateBlocked(valTVM);
    }
}       