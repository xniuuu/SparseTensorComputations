/*void blockedCOO(const vector<vector<int>> &COO, const vector<int> &n_d, const vector<double>  &val, const int block_size, 
                vector<vector<vector<int> > > &COO_block, vector<vector<double> > &val_block){

    unsigned int nnz = COO.size();
    int dim = COO[0].size(); 
    if (n_d.size() != dim)
        throw std::runtime_error("Dimension and array containing the dimension sizes do not match.");  

    unsigned int num_blocks = 1;  
    for (int i = 0; i < dim; i++) 
        num_blocks *= (n_d[i]%block_size == 0)? n_d[i]/block_size: n_d[i]/block_size + 1;

    COO_block.resize(num_blocks); //overcome memory: use a map instead of vector
    //check with map, if the block numer blick_i is already in map, if not construct block otherwise but it in existing block
    //use std::map or std::unordered map and copy it to std:vector then (to initialize with only non-empty blocks) or copy it in std::vector

    //64 bit integers for very large arrays (use size_t) typedef it to switch between 32 and 64 bit
    //implement sorting (maybe get rid of empty blocks first)
    //out of place sorting c++ 
    val_block.resize(num_blocks);

    unsigned int block_i = 0; 
    for (unsigned int row = 0; row < nnz; row++){
        vector<int> inds(COO[row].begin(), COO[row].begin()+dim);
        block_i = getBlock(inds, block_size, n_d);
        COO_block[block_i].push_back(COO[row]); //overwrite block_i with hilbert curve if it's used
        val_block[block_i].push_back(val[row]);
    }
}*/

/*
unsigned int getHilbertOrder(const std::vector<int>& indices, int hilbert_curve_order) {

    for (int i = 0; i < indices.size(); i++)
    {
        cout << "index: "<<indices[i] <<endl; 
    }
    
    unsigned int hilbert_order = 0;
    int d = indices.size();
    unsigned int pow2n = pow(2, hilbert_curve_order);

    std::vector<int> indices_permuted(d);
    std::vector<int> indices_rotated(d);

    int shift = pow2n / 2;
    int rotations = pow2n / 2;
    int level = 0;
    while (level < d) {
        for (int i = 0; i < d; i++) {
            indices_permuted[i] = indices[(i + level) % d];
            cout << "indices_permuted: "<<indices[i] <<endl;
        }

        for (int i = 0; i < d; i++) {
            indices_rotated[i] = indices_permuted[(i + rotations) % d];
            cout << "indices_rotated: "<<indices_rotated[i] <<endl;
        }

        int index = 0;
        for (int i = 0; i < d; i++) {
            index += indices_rotated[i] * pow(pow2n, i);
            cout << "index 2: "<<index <<endl;
        }

        hilbert_order += index * pow(pow2n, level) / pow2n;
        rotations /= 2;
        shift /= 2;

        if (rotations == 0) {
            rotations = pow2n / 2;
            shift = pow2n / 2;
            level++;
        }
    }

    return hilbert_order;
}
*/



/*COO::COO(const string file_name, char sep){
    string line, elem;
 
    fstream file_tns (file_name, ios::in);
    if(file_tns.is_open()){
        while(getline(file_tns, line)){
        stringstream entry(line);
        vector<int> coo;
        vector<double> value;
            while(getline(entry, elem, sep)){
                coo.push_back(stoi(elem));
                value.push_back(stod(elem));
            }
        val.push_back(value.back());
        coo.pop_back(); // Remove the last element (column)
        indices.push_back(coo);
        }
    }
    else
        cout<<"Could not open the file\n";
    
}*/


/*
BlockedCOO::BlockedCOO(const string file_name_tns, const string file_name_vec, char sep = " ", int b, vector<int> n_d){

    dim = n_d.size();
    tensor_size = n_d; 
    for (int i = 0; i < dim; i++)
    {
        cout << tensor_size[i]<<endl; 
    }
    
    unsigned int num_blocks = 1;  
    for (int i = 0; i < dim; i++) {
        unsigned int n_blocks_i = (tensor_size[i] + block_size - 1) / block_size;
        num_blocks *= n_blocks_i;
    }

    blocked_indices.resize(num_blocks);
    blocked_val.resize(num_blocks);

    unsigned int block_i = 0; 

    string line_tns, elem_tns, line_vec, elem_vec;
 
    fstream file_tns (file_name_tns, ios::in);
    fstream file_vec (file_name_vec, ios::in);

    if(file_tns.is_open() && file_vec.is_open()){
        while(getline(file_tns, line_tns) && getline(file_vec, line_vec)){
            nnz +=1; 

            stringstream entry(line_tns);
            stringstream entry_vec(line_vec);

            vector<int> coo;
            vector<double> value;
            vector<double> vec; 
                while(getline(entry, elem_tns, sep)){
                    coo.push_back(stoi(elem_tns));
                    value.push_back(stod(elem_tns));
                }
                while(getline(entry_vec, elem_vec, sep)){
                    vec.push_back(stod(elem_vec));
                }

        block_i = getBlock(coo, n_d);
        coo.pop_back(); // Remove the last element (column)
        blocked_indices[block_i].push_back(coo);
        blocked_val[block_i].push_back(value.back());
        blocked_vec[block_i].push_back(vec.front());
        }
    }
    else
        cout<<"Error opening files \n";

    for (auto it = blocked_indices.begin(); it != blocked_indices.end(); ) {
        if (it->empty()) {
            it = blocked_indices.erase(it);
        } else {
            ++it;
        }
    }
    for (auto it = blocked_val.begin(); it != blocked_val.end(); ) {
        if (it->empty()) {
            it = blocked_val.erase(it);
        } else {
            ++it;
        }
    }
    for (auto it = blocked_vec.begin(); it != blocked_vec.end(); ) {
        if (it->empty()) {
            it = blocked_vec.erase(it);
        } else {
            ++it;
        }
    }
    
}


BlockedCOO::BlockedCOO(const string file_name_tns, const string file_name_vec, char sep, int b, vector<int> n_d, const string traversal_curve){

    dim = n_d.size();
    tensor_size = n_d; 
    block_size = b; 
    nnz = 0; 
    for (int i = 0; i < dim; i++)
    {
        //cout << tensor_size[i]<<endl; 
    }
    
    vector<unsigned int> n_block(dim, 0); //stores the maximum number of blocks along a dimension

    unsigned int num_blocks = 1;  
    for (int i = 0; i < dim; i++) {
        unsigned int n_blocks_i = (tensor_size[i] + block_size - 1) / block_size;
        num_blocks *= n_blocks_i;
        n_block[i] = n_blocks_i;
    }

    unsigned int block_i = 0; 
    map<unsigned int, vector<vector<int> > > block_map_COO;
    map<unsigned int, vector<double > > block_map_val;
    map<unsigned int, vector<double > > block_map_vec;

    string line_tns, elem_tns, line_vec, elem_vec;
 
    fstream file_tns (file_name_tns, ios::in);
    fstream file_vec (file_name_vec, ios::in);

    if(file_tns.is_open() && file_vec.is_open()){
        while(getline(file_tns, line_tns) && getline(file_vec, line_vec)){
            nnz +=1; 

            stringstream entry(line_tns);
            stringstream entry_vec(line_vec);

            vector<int> coo;
            vector<double> value;
            vector<double> vec; 
                while(getline(entry, elem_tns, sep)){
                    coo.push_back(stoi(elem_tns));
                    value.push_back(stod(elem_tns));
                }
                while(getline(entry_vec, elem_vec, sep)){
                    vec.push_back(stod(elem_vec));
                }
            coo.pop_back(); // Remove the last element (column)
            block_i = getBlock(coo, tensor_size);
            block_map_COO[block_i].push_back(coo);
            //cout << block_i << endl; 
            block_map_val[block_i].push_back(value.back());
            block_map_vec[block_i].push_back(vec.front());
            //cout << value.back() << endl; 
        }

        for (const auto& block : block_map_val){
            for (const auto& entry : block.second) {
                cout << entry << endl;
            }
            cout << endl; 
        }

        for (const auto& block : block_map_vec){
            for (const auto& entry : block.second) {
                cout << entry << endl;
            }
            cout << endl; 
        }

        for (const auto& block : block_map_COO){
            for (const auto& entry : block.second) {
                for (auto& element : entry) {
                cout << element << " ";
            }
            cout << endl; 
            }
            cout << endl; 
        }
    }
    else
        cout<<"Error opening files \n";

    sortBlockedCOOZCurve(dim, n_block, block_map_COO, block_map_val, block_map_vec);

}*/


/*int main(){
    //install gdb!!!!!!
    //valgrind (checker for memory)

    vector<vector<int>> COO = csvToCOOIndices("test_tensor.csv",' ');
    vector<double> COO_val = csvToCOOValues("test_tensor.csv", ' ');

    //vector<int> n_d = {1605, 4198, 1631, 4209, 868131};
    vector<int> n_d = {16, 4, 8};
    int block_size = 4; 

    vector<vector<vector<int> > > COO_block;
    vector<vector<double> > val_block;

    blockedCOO(COO, COO_val, n_d, block_size, COO_block, val_block);

    
    for (const auto &vec1 : COO_block) {
        for (const auto &vec2 : vec1) {
            for (const auto &num : vec2) {
                std::cout << num << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    for (const auto &vec1 : val_block) {
        for (const auto &num : vec1) {
                std::cout << num << " ";
            }
            std::cout << std::endl;
        }
    
    BitIf lbln_bitif(COO, COO_val); 
    
    lbln_bitif.printIncrementArrays();
    lbln_bitif.printBitEncoder();

    int dim_TVM = lbln_bitif.getDim()- 1; 

    unsigned int* num_incr = lbln_bitif.getNumIncrements();
    double vec[16] = {1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0, 
                    1.0, 1.0, 1.0, 1.0,
                 1.0, 1.0, 1.0, 1.0};
    int** delta_tTVM = new int*[dim_TVM];

    for (int i = 0; i < dim_TVM; i++)
        delta_tTVM[i] = new int[num_incr[i]];
    

    vector<double> TVMval; 
    boost::dynamic_bitset<> TVMbit_encoder;
    lbln_bitif.computeTVM(vec, delta_tTVM, TVMval, TVMbit_encoder);

    cout<<"[";
    for (int i = 0; i < TVMval.size(); i++)
    {
        cout<<TVMval[i]<<" ";
    }
    cout<<"]";

    cout<<endl;

    cout<<"[";
    for (int i = 0; i < TVMbit_encoder.size(); i++)
    {
        cout<<TVMbit_encoder[i]<<" ";
    }
    cout<<"]"<<endl;

    for (int i = 0; i < dim_TVM; i++) 
        delete [] delta_tTVM[i];
 
    //delete [] delta_tTVM;
    //delete [] num_incr;
    lbln_bitif.removeBitIf();
    return 0;
}*/


/*
BlockedCOO::BlockedCOO(const string file_name_tns, const string file_name_vec, const char sep, vector<int> n_d, 
                        int b, int mode, const string traversal_curve){
    this->dim = n_d.size();
    this->tensor_size = n_d; 
    this->block_size = b; 
    this->nnz = 0; 

    unsigned int num_blocks = 1;  
    unsigned int block_i = 0; 

    string line_tns, elem_tns, line_vec, elem_vec;
 
    fstream file_tns (file_name_tns, ios::in);
    fstream file_vec (file_name_vec, ios::in);

    if(traversal_curve == "hierarchical"){
        for (int i = 0; i < this->dim; i++) {
            unsigned int n_blocks_i = (this->tensor_size[i] + block_size - 1) / block_size;
            num_blocks *= n_blocks_i;
        }
        
        blocked_indices.resize(num_blocks);
        blocked_val.resize(num_blocks);
        blocked_vec.resize(num_blocks);

        if(file_tns.is_open() && file_vec.is_open()){
            while(getline(file_tns, line_tns) && getline(file_vec, line_vec)){
                nnz +=1; 

                stringstream entry(line_tns);
                stringstream entry_vec(line_vec);

                vector<int> coo;
                vector<double> value;
                vector<double> vec; 
                    while(getline(entry, elem_tns, sep)){
                        coo.push_back(stoi(elem_tns));
                        value.push_back(stod(elem_tns));
                    }
                    while(getline(entry_vec, elem_vec, sep)){
                    
                        vec.push_back(stod(elem_vec));
                    }
 
            coo.pop_back(); // Remove the last element (column)
            block_i = getBlock(coo, this->tensor_size);
            
            if(mode>coo.size() && mode > 0){
                cout << "Invalid mode, the mode cannot be greather than the order of the tensor \n";
                break; 
            }
            if (coo.size() >= 2 && mode < coo.size() && mode > 0) {
                int temp = coo[mode-1];

                for (int i = mode-1; i < coo.size(); i++) {
                    coo[i] = coo[i+1];   
                }
                coo[coo.size() - 1] = temp;
            }

            blocked_indices[block_i].push_back(coo);
            //cout<< block_i<<endl; 
            blocked_val[block_i].push_back(value.back());
            blocked_vec[block_i].push_back(vec.front());
            }
        }
        else
            cout<<"Error opening files \n";

        for (auto it = blocked_indices.begin(); it != blocked_indices.end(); ) {
            if (it->empty()) {
                it = blocked_indices.erase(it);
            } else {
                ++it;
            }
        }
        for (auto it = blocked_val.begin(); it != blocked_val.end(); ) {
            if (it->empty()) {
                it = blocked_val.erase(it);
            } else {
                ++it;
            }
        }
        for (auto it = blocked_vec.begin(); it != blocked_vec.end(); ) {
            if (it->empty()) {
                it = blocked_vec.erase(it);
            } else {
                ++it;
            }
        }

        for (auto& x : blocked_indices) {
        for (auto& y : x) {
            for (auto& z : y) {
                std::cout << z << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    }
    else{
        map<unsigned int, vector<vector<int> > > block_map_COO;
        map<unsigned int, vector<double > > block_map_val;
        map<unsigned int, vector<double > > block_map_vec;

        vector<unsigned int> n_block(this->dim, 0);
        for (int i = 0; i < this->dim; i++) {
            unsigned int n_blocks_i = (this->tensor_size[i] + block_size - 1) / block_size;
            num_blocks *= n_blocks_i;
            n_block[i] = n_blocks_i;
        }

        if(file_tns.is_open() && file_vec.is_open()){
            while(getline(file_tns, line_tns) && getline(file_vec, line_vec)){
                nnz +=1; 

                stringstream entry(line_tns);
                stringstream entry_vec(line_vec);

                vector<int> coo;
                vector<double> value;
                vector<double> vec; 
                    while(getline(entry, elem_tns, sep)){
                        coo.push_back(stoi(elem_tns));
                        value.push_back(stod(elem_tns));
                    }
                    while(getline(entry_vec, elem_vec, sep)){
                        vec.push_back(stod(elem_vec));
                    }
                coo.pop_back(); // Remove the last element (column)
                block_i = getBlock(coo, this->tensor_size);

                if(mode>coo.size() && mode > 0){
                    cout << "Invalid mode, the mode cannot be greather than the order of the tensor \n";
                    break; 
                }
                if (coo.size() >= 2 && mode < coo.size() && mode > 0) {
                    int temp = coo[mode-1];

                    for (int i = mode-1; i < coo.size(); i++) {
                        coo[i] = coo[i+1];   
                    }
                    coo[coo.size() - 1] = temp;
                }

                block_map_COO[block_i].push_back(coo);
                //cout << block_i << endl; 
                block_map_val[block_i].push_back(value.back());
                block_map_vec[block_i].push_back(vec.front());
                //cout << value.back() << endl; 
            }

            for (const auto& block : block_map_val){
                for (const auto& entry : block.second) {
                    cout << entry << endl;
                }
                cout << endl; 
            }

            for (const auto& block : block_map_vec){
                for (const auto& entry : block.second) {
                    cout << entry << endl;
                }
                cout << endl; 
            }

            for (const auto& block : block_map_COO){
                for (const auto& entry : block.second) {
                    for (auto& element : entry) {
                    cout << element << " ";
                }
                cout << endl; 
                }
                cout << endl; 
            }
        }
        else{
            cout<<"Error opening files \n";
        }

        if(traversal_curve == "zcurve"){
            sortBlockedCOOZCurve(n_block, block_map_COO, block_map_val, block_map_vec);
        }
        else if(traversal_curve == "hilbert"){
            sortBlockedCOOHilbert(n_block, block_map_COO, block_map_val, block_map_vec);
        }
        else{
             cout<< "Invalid traversing curve\n";   
        }

    }

}

*/

/*
void BlockedCOO::sortBlockedCOOZCurve(vector<unsigned int> n_block, map<unsigned int, vector<vector<int>>> &block_map_COO, 
                                    map<unsigned int, vector<double>> &block_map_val, map<unsigned int, vector<double> >& block_map_vec){
    bool are_keys_equal = true;
    for (const auto& kv : block_map_COO) {
        if (block_map_val.find(kv.first) == block_map_val.end()) {
            are_keys_equal = false;
            break;
        }
    }
    if (are_keys_equal) {
        for (const auto& kv : block_map_val) {
            if (block_map_COO.find(kv.first) == block_map_COO.end()) {
                are_keys_equal = false;
                break;
            }
        }
    }

    if (!are_keys_equal) {
        throw std::runtime_error("Blocked COO storage and associated blocked value storage does not match.");  
    }

    unsigned int nneb = block_map_COO.size(); //nneb = number of non-empy blocks
    blocked_indices.resize(nneb);
    blocked_val.resize(nneb);
    blocked_vec.resize(nneb);

    unsigned int z_number = 0; 
    map<unsigned int, unsigned int > z_order;
    unsigned int z_count = 0; 

    for (const auto& block : block_map_val){
        unsigned int block_i = block.first;
        z_number = returnZOrder(n_block, block_i);
        z_order[z_number] = block_i;
        //cout << "z: "<< z_number<< endl; 
    }

    for (const auto& num : z_order){
        //cout << "z order: "<< num.first << " " << num.second << endl; 
    }

    if (!(nneb == z_order.size())) {
        throw std::runtime_error("Ordering map of z curve does not match coo block map");  
    }

    unsigned int count = 0; 
    for (const auto& num : z_order){
        unsigned int block_i = num.second; 
        blocked_indices_2d.insert(blocked_indices_2d.end(), block_map_COO[block_i].begin(), block_map_COO[block_i].end());
        blocked_val_2d.insert(blocked_val_2d.end(), block_map_val[block_i].begin(), block_map_val[block_i].end());

        blocked_indices[count] = block_map_COO[block_i];
        blocked_val[count] = block_map_val[block_i];
        blocked_vec[count] = block_map_vec[block_i];
        count +=1; 
    }
    
        // Print the modified 3D vector
    for (auto& dim1 : blocked_indices) {
        for (auto& dim2 : dim1) {
            for (auto& element : dim2) {
                cout << element << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    for (auto& dim1 : blocked_val) {
        for (auto& element : dim1) {
                cout << element << " ";
            }
            cout << endl;
        }
    for (auto& dim1 : blocked_vec) {
        for (auto& element : dim1) {
                cout << element << " ";
            }
            cout << endl;
        }
}
*/

/*
void BlockedCOO::sortBlockedCOOHilbert(vector<unsigned int> n_block, map<unsigned int, vector<vector<int>>> &block_map_COO, map<unsigned int, vector<double>> &block_map_val, map<unsigned int, vector<double>> &block_map_vec){
    bool are_keys_equal = true;
    for (const auto& kv : block_map_COO) {
        if (block_map_val.find(kv.first) == block_map_val.end()) {
            are_keys_equal = false;
            break;
        }
    }
    if (are_keys_equal) {
        for (const auto& kv : block_map_val) {
            if (block_map_COO.find(kv.first) == block_map_COO.end()) {
                are_keys_equal = false;
                break;
            }
        }
    }

    if (!are_keys_equal) {
        throw std::runtime_error("Blocked COO storage and associated blocked value storage does not match.");  
    }

    unsigned int nneb = block_map_COO.size(); //nneb = number of non-empy blocks
    blocked_indices.resize(nneb);
    blocked_val.resize(nneb);
    blocked_vec.resize(nneb);

    unsigned int hilbert_number = 0; 
    map<unsigned int, unsigned int > hilbert_order;

    unsigned int hilbert_curve_order = 3; 

    for (const auto& block : block_map_val){
        unsigned int block_i = block.first;
        hilbert_number = returnHilbertNumber(n_block, block_i, hilbert_curve_order);
        hilbert_order[hilbert_number] = block_i;
        //cout << "z: "<< z_number<< endl; 
    }

    for (const auto& num : hilbert_order){
        //cout << "z order: "<< num.first << " " << num.second << endl; 
    }

    if (!(nneb == hilbert_order.size())) {
        throw std::runtime_error("Ordering map of hilbert does not match coo block map");  
    }

    unsigned int count = 0; 
    for (const auto& num : hilbert_order){
        unsigned int block_i = num.second; 
        blocked_indices_2d.insert(blocked_indices_2d.end(), block_map_COO[block_i].begin(), block_map_COO[block_i].end());
        blocked_val_2d.insert(blocked_val_2d.end(), block_map_val[block_i].begin(), block_map_val[block_i].end());

        blocked_indices[count] = block_map_COO[block_i];
        blocked_val[count] = block_map_val[block_i];
        blocked_vec[count] = block_map_vec[block_i];
        count +=1; 
    }

            // Print the modified 3D vector
    for (auto& dim1 : blocked_indices) {
        for (auto& dim2 : dim1) {
            for (auto& element : dim2) {
                cout << element << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    for (auto& dim1 : blocked_val) {
        for (auto& element : dim1) {
                cout << element << " ";
            }
            cout << endl;
        }
    for (auto& dim1 : blocked_vec) {
        for (auto& element : dim1) {
                cout << element << " ";
            }
            cout << endl;
        }

}
*/


/*
//old from blockedCOO constructor

    if(traversal_curve == "hierarchical"){
        //cout << num_blocks<<endl; 
        blocked_indices.resize(num_blocks);
        blocked_val.resize(num_blocks);

        if(file_tns.is_open()){
            while(getline(file_tns, line_tns)){
                nnz +=1; 

                stringstream entry(line_tns);

                vector<int> coo;
                vector<double> value;
                    while(getline(entry, elem_tns, sep)){
                        coo.push_back(stoi(elem_tns));
                        value.push_back(stod(elem_tns));
                    }

            coo.pop_back(); // Remove the last element (column)
            block_i = getBlock(coo, this->tensor_size);
            
            if(mode>coo.size() && mode > 0){
                cout << "Invalid mode, the mode cannot be greather than the order of the tensor \n";
                break; 
            }
            if (coo.size() >= 2 && mode < coo.size() && mode > 0) {
                int temp = coo[mode-1];

                for (int i = mode-1; i < coo.size(); i++) {
                    coo[i] = coo[i+1];   
                }
                coo[coo.size() - 1] = temp;
            }

            blocked_indices[block_i].push_back(coo);
            //cout<< block_i<<endl; 
            blocked_val[block_i].push_back(value.back());
            }
        }
        else
            cout<<"Error opening files \n";

        for (auto it = blocked_indices.begin(); it != blocked_indices.end(); ) {
            if (it->empty()) {
                it = blocked_indices.erase(it);
            } else {
                ++it;
            }
        }
        for (auto it = blocked_val.begin(); it != blocked_val.end(); ) {
            if (it->empty()) {
                it = blocked_val.erase(it);
            } else {
                ++it;
            }
        }

        for (auto& x : blocked_indices) {
        for (auto& y : x) {
            for (auto& z : y) {
                std::cout << z << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    }
    else {

    //cout << num_blocks<<endl; 
        if(file_tns.is_open()){
            while(getline(file_tns, line_tns)){
                nnz +=1; 

                stringstream entry(line_tns);

                vector<int> coo;
                vector<double> value;
                    while(getline(entry, elem_tns, sep)){
                        coo.push_back(stoi(elem_tns));
                        value.push_back(stod(elem_tns));
                    }
                coo.pop_back(); // Remove the last element (column)
                block_i = getBlock(coo, this->tensor_size);

                if(mode>coo.size() && mode > 0){
                    cout << "Invalid mode, the mode cannot be greather than the order of the tensor \n";
                    break; 
                }
                if (coo.size() >= 2 && mode < coo.size() && mode > 0) {
                    int temp = coo[mode-1];

                    for (int i = mode-1; i < coo.size(); i++) {
                        coo[i] = coo[i+1];   
                    }
                    coo[coo.size() - 1] = temp;
                }

                block_map_COO[block_i].push_back(coo);
                //cout << block_i << endl; 
                block_map_val[block_i].push_back(value.back());
                //cout << value.back() << endl; 
            }
            cout << block_map_val.size()<<endl;
            for (const auto& block : block_map_val){
                for (const auto& entry : block.second) {
                    cout << entry << endl;
                }
                cout << endl; 
            }

            for (const auto& block : block_map_COO){
                for (const auto& entry : block.second) {
                    for (auto& element : entry) {
                    cout << element << " ";
                }
                cout << endl; 
                }
                cout << endl; 
            }
        }
        else{
            cout<<"Error opening files \n";
        }
*/


/*
COO::COO(const string file_name, bool offset, const char sep, vector<unsigned int> n_d, int mode, const string traversal_curve){

    this->tensor_size = n_d;                        
    if (n_d.size() >= 2 && mode < n_d.size() && mode > 0) {
        int temp = tensor_size[mode-1];
        for (int i = mode-1; i < n_d.size(); i++) {
            this->tensor_size[i] = this->tensor_size[i+1];   
        }
        this->tensor_size[n_d.size() - 1] = temp;
    }
    //for (int i = 0; i < n_d.size(); i++) 
    //    cout<<tensor_size[i]<<" ";
    //cout<<endl;

    unsigned int order = 0;
    map<unsigned int, vector<vector<int> > > map_COO;
    map<unsigned int, vector<double > > map_val;

    string line, elem;
    this->nnz = 0; 
    fstream file_tns (file_name, ios::in);
    if(file_tns.is_open()){        
        while(getline(file_tns, line)){
        this->nnz+=1;

        stringstream entry(line);
        vector<int> coo;
        vector<double> value;
            while(getline(entry, elem, sep)){
                coo.push_back(stoi(elem)-offset);
                value.push_back(stod(elem));
            }
        
        coo.pop_back(); // Remove the last element (column)
        //put value at mode k at the end:
        if(mode>coo.size() && mode > 0){
            cout << "Invalid mode, the mode cannot be greather than the order of the tensor \n";
            break; 
        }
        if (coo.size() >= 2 && mode < coo.size() && mode > 0) {
            int temp = coo[mode-1];

            for (int i = mode-1; i < coo.size(); i++) {
                coo[i] = coo[i+1];   
            }
            coo[coo.size() - 1] = temp;
        }

        
        if(traversal_curve == "hierarchical"){
            val.push_back(value.back());
            indices.push_back(coo);
        }

        else if(traversal_curve == "zcurve"){
            order = returnZOrder(coo); 
            map_COO[order].push_back(coo);
            map_val[order].push_back(value.back());
            //cout << "hilbert_number: "<< hilbert_number<< endl; 
        
             
        } else if(traversal_curve == "hilbert"){
            order = returnHilbertNumber(coo, 3);
            map_COO[order].push_back(coo);
            map_val[order].push_back(value.back());
        } else{
            cout<< "Invalid traversing curve\n";   
        }

    }

    unsigned int total_entries = 0;
    for (const auto& entry : map_COO) {
        total_entries += entry.second.size();
    }
    if (!(nnz == total_entries)) {
        throw std::runtime_error("Numbers of entries in map does not match number of non zeros");  
    }

    if(traversal_curve != "hierarchical"){
        for (const auto& entries : map_COO){
            const vector< vector<int> >& current_entries = entries.second;
            for (auto entry : current_entries){
                indices.push_back(entry);
            }
        }
        for (const auto& entries : map_val){
            const vector< double >& current_entries = entries.second;
            for (auto entry : current_entries){
                val.push_back(entry);
            }
        }

    }
        for (const auto& row : indices) {
             for (const auto& elem : row) {
                cout << elem << " ";
        }
        cout << endl;
        }  
        //cout << this->nnz<<endl;
    }
    else
        cout<<"Could not open the file\n";


}

*/


/*
BitIf::BitIf(const std::vector<std::vector<int> > COO, const std::vector<double> COO_values) {
    
    this->nnz = COO_values.size();
    this->dim = COO[0].size();
    this->delta_d = new int*[dim]();
    this->val = new double[nnz]();
    this->bit_encoder.resize(dim*nnz);
    unsigned int *incr_counter = new unsigned int[dim]();
    this->num_increments = new unsigned int[dim]();

    for(int i = 0; i<dim; i++){
        unsigned int count = 1; 
        for(unsigned int row = 0; row <nnz-1; row++){
            
            int delta_j = COO[row+1][i] - COO[row][i];
            if (delta_j !=0)
               count +=1;      
        } 
        this->num_increments[i] = count;   
    }


    this->val[0] = COO_values[0];
    for (int j = 0; j < dim; j++){
        this->delta_d[j] = new int[num_increments[j]]();
        this->delta_d[j][0] = COO[0][j];
        incr_counter[j] = 1;
        bit_encoder[j] = 1;
    }

    for(int row = 1; row<nnz; row++){
        this->val[row] = COO_values[row];
        for(int j = 0; j< dim ; j++){
            int delta_j = COO[row][j] - COO[row-1][j];
            if (delta_j !=0){
                this->delta_d[j][incr_counter[j]] = delta_j;
                incr_counter[j]++;
                this->bit_encoder[dim*row+j] = 1;
            }
        }
    }
}
*/

/*
BlockedBitIf::BlockedBitIf(const std::vector< std::vector<std::vector<int> > > blocked_indices, 
                            const std::vector<std::vector<double> > blocked_vals){

    this->num_blocks = blocked_indices.size(); 
    this->dim = blocked_indices[0][0].size();
    this->delta_d = new int**[num_blocks]();
    this->val = new double*[num_blocks]();
    this->bit_encoder.resize(num_blocks);

    //Allocate memory for blocked increment arrays and blocked tensor value arrays
    for (unsigned int b = 0; b < num_blocks; b++) {
        int block_size = blocked_indices[b].size();
        if (block_size == 0)
            continue;       
        
        //store the unique indices up to mode k-1 (for arbitrary traversal TVM):
        for (int i = 0; i < block_size; i++) {
            assert(blocked_indices[b][i].begin() != blocked_indices[b][i].end());
            vector<int> current_indices(blocked_indices[b][i].begin(), blocked_indices[b][i].end() - 1);
            //current_indices.push_back(0); 
            if (find(unique_COO.begin(), unique_COO.end(), current_indices) == unique_COO.end()) {
                unique_COO.push_back(move(current_indices)); //double check
            }
        }

        //allocate memory for the class members
        this->delta_d[b] = new int*[dim]();
        this->val[b] = new double[block_size]();
        unsigned int *incr_counter_b = new unsigned int[dim]();
        unsigned int *num_increments = new unsigned int[dim]();
        this->bit_encoder[b].resize(dim*block_size);

        for(int i = 0; i<dim; i++){
            unsigned int count = 1; 
            for(unsigned int row = 0; row <block_size-1; row++){
                int delta_j = blocked_indices[b][row+1][i] - blocked_indices[b][row][i];
                if (delta_j !=0)
                    count +=1;      
                } 
            num_increments[i] = count;   
        }

        this->val[b][0] = blocked_vals[b][0];
        for (int j = 0; j < dim; j++){
            this->delta_d[b][j] = new int[num_increments[j]]();
            this->delta_d[b][j][0] = blocked_indices[b][0][j];
            incr_counter_b[j] = 1;
            this->bit_encoder[b][j] = 1;
        }

        for(int row = 1; row<block_size; row++){
            this->val[b][row] = blocked_vals[b][row];
            for(int j = 0; j< dim ; j++){
                int delta_j = blocked_indices[b][row][j] - blocked_indices[b][row-1][j];
                if (delta_j !=0){
                    this->delta_d[b][j][incr_counter_b[j]] = delta_j;
                    incr_counter_b[j]++;
                    this->bit_encoder[b][dim*row+j] = 1;
                }
            }
        }
        
        if(incr_counter_b !=nullptr)
            delete [] incr_counter_b;
        if(num_increments !=nullptr)
            delete [] num_increments;
    
    }
}
*/



//////////////////////////////old BITIF code//////////////////////////

    /*if(TVMnum_incr !=nullptr){
        for (unsigned int b = 0; b < this->TVMnum_blocks; b++) {
            delete [] TVMnum_incr[b];
        }
        delete [] TVMnum_incr;
    }
    
    /*this-> TVMdim = this->dim - 1; 
    vector<unsigned int> incr_counter(TVMdim);
    unsigned int TVMtotal_increments = 0; 
    unsigned int TVMnnz = 0; 
    vector<unsigned int> TVMincr_offset(TVMnum_blocks*TVMdim);
    //>cout<<TVMnum_blocks<<endl;
    for (unsigned int b = 0; b<TVMnum_blocks; ++b){
        unsigned int TVMsize_b = TVM_COO[b].size();
        //cout<<TVMsize_b<<endl;
        TVMnnz += TVMsize_b;
        for(int i = 0; i<TVMdim; i++){
            //unsigned int count = 1; 
            TVMincr_offset[b*TVMdim + i] = TVMtotal_increments;
            for(unsigned int row = 0; row <TVMsize_b-1; row++){
                int delta_j = TVM_COO[b][row+1][i] - TVM_COO[b][row][i];
                if (delta_j !=0){
                    //count +=1;   
                    TVMtotal_increments += 1; 
                }   
            } 
        }
    }
    //cout<<TVMtotal_increments<<endl;
    boost::dynamic_bitset<> TVMbit_encoder_blocked(TVMnnz*TVMdim);
    vector<int> TVMdelta_d_blocked(TVMtotal_increments);

    unsigned int offset = 0; 
    unsigned int offset_b = 0; 
    unsigned int offset_d = 0; 
    for (unsigned int b = 0; b<TVMnum_blocks; ++b){
        unsigned int TVMsize_b = TVM_COO[b].size();
        for (int j = 0; j < TVMdim; j++){
            offset = TVMincr_offset[b*TVMdim + j];
            TVMdelta_d_blocked[offset] = TVM_COO[b][0][j];
            incr_counter[j] = 1;
            TVMbit_encoder_blocked[offset_b+j] = 1; 
        }
          
        for(int row = 1; row < TVMsize_b; row++){
            for(int j = 0; j< TVMdim ; j++){
                int delta_j = TVM_COO[b][row][j] - TVM_COO[b][row-1][j];
                if (delta_j !=0){
                    offset_d = TVMincr_offset[b*TVMdim + j] + incr_counter[j];
                    //cout<<TVMincr_offset[b*TVMdim + j]<<endl;
                    TVMdelta_d_blocked[offset_d] = delta_j;
                    incr_counter[j]++;
                    TVMbit_encoder_blocked[offset_b+ TVMdim*row+j] = 1;
                }
            }
        }
        offset_b += TVMsize_b*TVMdim;
    }

    TVMend = std::chrono::high_resolution_clock::now();
    this->TVMruntime = chrono::duration_cast<chrono::milliseconds>(TVMend - TVMstart).count();
    cout<<this->TVMruntime <<endl; 
*/