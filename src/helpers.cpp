#include <cmath>
#include <stdexcept>
#include "helpers.h"

using namespace std;

vector<vector<int> > csvToCOOIndices(const string file_name, char sep){
    vector<vector<int> > arr; 
    string line, elem;
 
    fstream file (file_name, ios::in);
    if(file.is_open()){
        while(getline(file, line)){
        stringstream entry(line);
        vector<int> row;
            while(getline(entry, elem, sep)){
                row.push_back(stod(elem));
            }
        row.pop_back(); // Remove the last element (column)
        arr.push_back(row);
        }
    }
    else
        std::cerr << "Error: Failed to open file: " << file_name << std::endl;
    
    return arr;
}


vector<double> csvToCOOValues(const string file_name, char sep){
    vector<double> vals; 
    string line, elem;
 
    fstream file (file_name, ios::in);
    if(file.is_open()){
        while(getline(file, line)){
        stringstream entry(line);
        vector<double> row;
            while(getline(entry, elem, sep)){
                row.push_back(stod(elem));
            }
        vals.push_back(row.back());
        }
    }
    else
        std::cerr << "Error: Failed to open file: " << file_name << std::endl;
    
    return vals;
}

vector<double> getVecFromFile(const string file_name){
    vector<double> vec;
    string line;
    fstream file (file_name, ios::in);
    if (file.is_open()) {
        while (getline(file, line)) {
            if (line[0] != '#') {
                vec.push_back(stod(line));
            }
        }
    } else {
        std::cerr << "Error: Failed to open file: " << file_name << std::endl;
    }
    return vec;
}

void printTensorValues(vector<double> val){
    cout<<"Tensor values:"<<"\n"<< "["; 
    for(int i=0;i<val.size();i++)
        cout<<val[i]<<" ";
    std::cout<<"]"<<"\n";
}

unsigned int* determineIncrements(vector<vector<int>> COO){
    unsigned int nnz = COO.size();
    int dims = COO[0].size(); 
    cout << dims <<endl;
    cout<< nnz <<endl;

    unsigned int *counter = new unsigned int[dims];

    for(int i = 0; i<dims; i++){
        unsigned int count = 0; 
        for(unsigned int row = 0; row <nnz-1; row++){
            
            int delta_j = COO[row+1][i] - COO[row][i];
            if (delta_j !=0)
               count +=1;      
        } 
        counter[i] = count;   
    }
    cout<<endl;
    cout<<"successfully constructed counter *"<<endl;
    return counter; 
}


// -try it out with smaller scale 
// own code: nz*d*log(n_tensor) n_tensor = tensor size
// try find hilbert curve function from someone else
//try z-curve first: use binary to order; convert indices to binary 
        //convert block number to binary and sort them to z-curves
//google the details

