#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <chrono>

// TODO: Check which among map and unordered map gives the best result?
// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

double laplace_const = 0.0035;
/* In the given starter code accessing a node given the name or its index seems to return and iterator and takes
time since it has to traverse through entire graph list -> why not make it index based?*/

class network{
public:
    map<string,int> Name_ind; //Mapping for the name of node to index
    map<int,string> ind_Name; //Mapping for index to Name of node (incase ever needed)
    // vector<string> Nodes; //Names of the nodes
    vector<map<string,int>> cat_val; // Contains the category-value pairs that a particular node can take -> will be useful while trying to update CPT (for index in CPT)
    vector<vector<string>> Values; //Contains just the possible categories a variable can take
    vector<int> nValues; //number of values a node can take -> basically Values[i].size()
    vector<vector<int>> Children; //index of all children of a particular node
    vector<vector<int>> Parents; //index of all parents of a particular node
    vector<vector<double>> CPT; //The entire cpt table ordered according to index of node

    network() {
        Name_ind.clear();
        cat_val.resize(37);
        Values.resize(37); // Set the outer vector to have 37 elements
        nValues.resize(37); // Set the vector of nValues to have 37 elements
        Children.resize(37); // Set the vector of vectors for Children to have 37 elements
        Parents.resize(37); // Set the vector of vectors for Parents to have 37 elements
        CPT.resize(37); // Set the vector of vectors for CPT to have 37 elements
      // Now all the vectors will have 37 elements initially
    }

    // * reads the .bif file and intitializes Name_ind, parents, children, cpt
    void read_network(const std::string& filename){
        ifstream file(filename);
        string line;
        int find=0;
        string temp;
        string name;
        map<string,int> category_val;
        vector<string> values;
        int index=0;
        if (file.is_open())
        {
            while (! file.eof() )
            {
                stringstream ss;
                getline (file,line);
                
                
                ss.str(line);
                ss>>temp;
                
                //Updating the mapping for variablename to Index (i.e. Name_ind) and The categories/values a variable can take
                if(temp.compare("variable")==0)
                {
                        
                        ss>>name;
                        getline (file,line);
                    
                        stringstream ss2;
                        ss2.str(line);
                        for(int i=0;i<4;i++)
                        {
                            
                            ss2>>temp;
                            
                            
                        }
                        category_val.clear();
                        values.clear();
                        int val = 0;
                        while(temp.compare("};")!=0)
                        {
                            category_val.emplace(temp,val);
                            values.push_back(temp);
                            val++;
                            ss2>>temp;
                        }
                        Name_ind.emplace(name,index);
                        ind_Name.emplace(index,name);
                        cat_val[index] = category_val;
                        Values[index] = values;
                        nValues[index] = val;
                        index+=1;
                        // Graph_Node new_node(name,values.size(),values);
                        // int pos=Alarm.addNode(new_node);

                        
                }
                else if(temp.compare("probability")==0)
                {
                        
                        ss>>temp;
                        ss>>temp;
                        
                        // temp now contains the current node we are at whose req index is accessible using map[temp]
                        int node_index = Name_ind[temp];
                        // list<Graph_Node>::iterator listIt;
                        // list<Graph_Node>::iterator listIt1;
                        // listIt=Alarm.search_node(temp);
                        // int index=Alarm.get_index(temp);

                        ss>>temp; // going to next space-delimited string
                        //getting all parents
                        vector<int> temp_parents;
                        while(temp.compare(")")!=0) 
                        {
                            temp_parents.push_back(Name_ind[temp]); //Adding parent index for current node into a temp vector
                            // listIt1=Alarm.search_node(temp);
                            // listIt1->add_child(index);
                            Children[Name_ind[temp]].push_back(node_index); //Adding current node index as child for the parent
                            // values.push_back(temp);
                            ss>>temp;

                        }
                        Parents[node_index] = temp_parents;
                        // listIt->set_Parents(values);
                        getline (file,line); //table line
                        stringstream ss2;
                        
                        ss2.str(line);
                        ss2>> temp; //the word "table" (acc to alarm.bif)
                        
                        ss2>> temp; //cpt val one by one
                        
                        vector<double> curr_CPT;
                        string::size_type sz;
                        while(temp.compare(";")!=0)
                        {
                            curr_CPT.push_back(atof(temp.c_str()));
                            ss2>>temp;
                        }
                        CPT[node_index] = curr_CPT;
                        // listIt->set_CPT(curr_CPT);
                }
                else
                {   
                }   
            }
            if(find==1)
            file.close();
        }
    }

    // * For writing into solved_alarm.bif
	void write_network(const std::string &filename)
	{
		ifstream alarm(filename);
		ofstream solved_alarm("solved_alarm.bif");
		string temp;
		if (alarm.is_open() && solved_alarm.is_open()) {
			string line;
			getline(alarm,line);
			while (!alarm.eof()) {
				solved_alarm << line << endl;
				stringstream ss;
				ss.str(line);
				ss>>temp;
				if (temp.compare("probability")==0) {
					ss>>temp; //contains the "("
                    ss>>temp; //contains the present node name
					int node_index = Name_ind[temp]; //contains index to nodename
					getline(alarm,line);
					stringstream ss1;
					ss1.str(line);
					ss1>>temp; //the word "table"
					solved_alarm << "\t" << temp << " ";
					vector<double> var_cpt = CPT[node_index];
					for (int i=0;i<var_cpt.size();i++){
						if (var_cpt[i] < 0.0001){
							solved_alarm << 0.0001 << " ";
						}
						else{
							solved_alarm << std::fixed << setprecision(5) <<var_cpt[i] << " ";
						}
					}
					solved_alarm << ";" << endl;
				}
				getline(alarm,line);
			}
			alarm.close();
			solved_alarm.close();
        cout << "File has been modified and saved as " << "solved_alarm.bif" << endl;
		} else {
			cerr << "Error opening files!" << endl;
		}
	}

    // * For checking the network intialized
    void view_network(){
        std::cout << "Name_ind:\n";
        for (const auto &entry : Name_ind) {
            std::cout << entry.first << ": " << entry.second << std::endl;
        }
        
        std::cout << "\nValues:\n";
        for (const auto &vec : Values) {
            for (const auto &num : vec) {
                std::cout << num << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "\nCat_values:\n";
        for (const auto &vec : cat_val) {
            for (const auto &entry : vec) {
                std::cout << entry.first << ": " << entry.second << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "\nnValues:\n";
        for (const auto &val : nValues) {
            std::cout << val << " ";
        }
        std::cout << std::endl;

        std::cout << "\nChildren:\n";
        for (const auto &vec : Children) {
            for (const auto &num : vec) {
                std::cout << num << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "\nParents:\n";
        for (const auto &vec : Parents) {
            for (const auto &num : vec) {
                std::cout << num << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "\nCPT:\n";
        for (const auto &vec : CPT) {
            for (const auto &num : vec) {
                std::cout << num << " ";
            }
            std::cout << std::endl;
        }
    }
};

/*
  TODO: Need to make a class for reading .dat file and having functions to update corresponding CPT and data values
  * This class should read the .dat file and store it into a vector of vector of strings ( (~10000) lines * 37 variables)
  * Possible other variables - position of '?' in the data file, ---what else?
  * Parameters to class would prbly be .dat file and pointer to network class instance 
  * Class should also have function for expectation, maximisation
  ? Do we use CPT of network class, or make a new variable that maps node_name (or ind) to the CPT table like thing, and reading each value from data we just add 1 to that particular element?
*/

class compute_CPT{

public:
    network* network_obj;
    vector<vector<int>> data_val; // stores the corresponding values of the categories using the network's Values variable
    vector<vector<string>> data_cat; //stores the corresponding categories from .datfile
    vector<int> q_ind; //stores index of question mark within each vector of data
    vector<double> prob_weight;

    compute_CPT(network* ptr){ // Constructor
        network_obj = ptr;
    }
    
    // Reading the data file
    void read_data(const std::string& filename){
        ifstream file(filename);
        string line;
        vector<int> vals(37); //corersponding values for categories of current line
        vector<string> cats(37); //corresponding categories of current line
        int find=0;
        string temp;
        string name;
        map<string,int> category_val;
        if (file.is_open())
        {
            while (! file.eof() )
            {
                stringstream ss;
                getline (file,line); // Each line in the .dat file
                ss.str(line);
                int index = 0;
                int q_index = -1;
                while(ss>>temp){ // Each word in the particular line
                    if (temp.compare("\"?\"")==0){
                        vals[index] = -1;
                        cats[index] = "?";
                        q_index = index;
                    }
                    else{
                        vals[index] = network_obj->cat_val[index][temp]; //Add the value of the category of a node to vals
                        cats[index] = temp;
                        // q_ind.push_back(-1);
                    }
                    index ++;
                }
                q_ind.push_back(q_index);
                data_val.push_back(vals);
                data_cat.push_back(cats);
            }
            prob_weight.resize(q_ind.size());
        }
        if(find==1)
        file.close();
    }


    /*
     * Possible Intializations:
     ? Random uniform distribution for P(Ai|P) for all i 
     ? 1/(number of possible values) for P(Ai|P) for all i 
     ? 1 for P(Ai|P) and 0 for P(Aj|P) (j!=i) where Ai,Aj's are possible values of A 
     ? Randomly chose one value from possible values and substitute "?" with it in data

     Note - First 3 deal with initalizing of CPT first whereas last one deals with data initialization first
     Need to test out everything to check which performs better -> ig easier to make an init_CPT function rather than modifying expectation everytime 
    */

    void intialize_CPT(){
        double no_of_vals;
        for (int i=0;i<37;i++){ //for all probability tables in CPT
            no_of_vals = network_obj->nValues[i]; //Gives the total number of values for the variable in question
            for (auto &prob: network_obj->CPT[i]){
                prob = 1/no_of_vals;
            }
        }
    }

    // * Update data_val and/or data_cat
    void expectation(){
        for (int i=0;i<q_ind.size();i++){
            int missing_idx = q_ind[i];
            if (missing_idx == -1){
                prob_weight[i] = -1;
                continue;
            }
            else{
                vector<int> data = data_val[i]; //the data of the line
                int no_of_val = network_obj->nValues[missing_idx]; //possible values that the ? variable can take
                vector<int> parents_idx = network_obj->Parents[missing_idx]; //parents of the ? variable

                //calculating values of all parents from data
                int par_size = parents_idx.size();
                vector<int> parents_val(par_size);
                for (int j=0;j<par_size;j++){
                    parents_val[j] = data[parents_idx[j]]; //Value of parents
                }

                //index to be accessed from cpt for P(var|parents)
                int var_cpt_idx = 0;
                int mult = 1;
                if (par_size > 0){
                    for (int j=par_size-1;j>=0;j--){
                        var_cpt_idx += mult*parents_val[j];
                        mult *= network_obj->nValues[parents_idx[j]];
                    }
                }

                //computing P(var|MB(var)) of child for all possible values of variable
                vector<int> children_idx = network_obj->Children[missing_idx]; //children of the ? variable
                vector<double> possible_prob(no_of_val,0.0);
                double prob = 1;
                double max_prob = 1;
                int prediction_val = 0;
                for (int val=0; val<no_of_val; val++){
                    data[missing_idx] = val;
                    var_cpt_idx += mult*val; //index of P(var=val|Parents)
                    
                    prob *= network_obj->CPT[missing_idx][var_cpt_idx]; //P(var=val|Parents)
                
                    //computing P(var's children given children's parents)
                    for (int child=0;child<children_idx.size();child++){ //for every child of variable
                        int child_idx = children_idx[child]; //index of child
                        int val_child = data[child_idx]; //value of child
                        
                        vector<int> child_parents_idx = network_obj->Parents[child_idx]; //index of parents of child
                        int child_par_size = child_parents_idx.size();
                        vector<int> child_parents_val(child_par_size);
                        for (int j=0;j<child_par_size;j++){
                            child_parents_val[j] = data[child_parents_idx[j]]; //Value of parents
                        }

                        //index to be accessed from cpt for child|parents
                        int child_cpt_idx = 0;
                        int mult = 1;
                        if (child_par_size > 0){
                            for (int j=child_par_size-1;j>=0;j--){
                                child_cpt_idx += mult*child_parents_val[j];
                                mult *= network_obj->nValues[child_parents_idx[j]];
                            }
                        }
                        child_cpt_idx += mult*val_child;
                        prob*=network_obj->CPT[child_idx][child_cpt_idx];
                    }

                    possible_prob[val] = prob;
                    if (max_prob < prob){
                        prediction_val = val;
                        max_prob = prob;
                    }
                    prob=1;
                }
                
                //Use the below commented out part in case we want to use maximimum probability for prediction
                // data_val[i][missing_idx] = prediction_val;  
                // prob_weight[i] = possible_prob[prediction_val];

                //Use the below one for expected value method
                double total_prob = 0;
                for (int j = 0;j<no_of_val;j++){
                    total_prob+=possible_prob[j];
                }

                //expectation calculation
                double expected_val=0;
                for (int var_val = 0; var_val<no_of_val;var_val++){
                    expected_val += ((double) possible_prob[var_val]* (double) var_val) / ((double) total_prob);
                }

                data_val[i][missing_idx] = round(expected_val); //update the value 
                prob_weight[i] = possible_prob[round(expected_val)];
            }
        }
    }

    // * Update cpt -> Assumes that expectation is called before maximization (only in the intial case-else there would be ? in data)
    void maximization(){

        // below part is for making it all 0's and then doing the necessary calculations
        for (int i=0;i<37;i++){ //for all probability tables in CPT
            network_obj->CPT[i] = vector<double>(network_obj->CPT[i].size(),laplace_const);
        }

        //CPT numerator updation -> 1 pass through data
        for (int line=0;line<data_val.size();line++){ //for each line
            // cout << "Entered loop iteration "<<line<<endl;
            vector<int> line_data = data_val[line]; // ith line
            int missing_idx = q_ind[line];
            for (int var=0;var<37;var++){ //for every variable
                int var_val = line_data[var]; //value of variable
                vector<int> parents_idx = network_obj->Parents[var]; //index of parents of the variable
                int par_size = parents_idx.size();
                vector<int> parents_val(par_size);
                for (int i=0;i<par_size;i++){
                    parents_val[i] = line_data[parents_idx[i]]; //Value of parents
                }

                //index in CPT to be updated would be given by below code:
                int var_cpt_idx = 0;
                int mult = 1;
                if (par_size > 0){
                    for (int i=par_size-1;i>=0;i--){
                        var_cpt_idx += mult*parents_val[i];
                        mult *= network_obj->nValues[parents_idx[i]];
                    }
                }
                var_cpt_idx += mult*var_val;
                // cout<<var_cpt_idx<<endl;
                //Now CPT of variable is network_obj->CPT[var]
                //Now add 1 to the var_cpt_idx index in CPT[var] to get numerator
                if (var == missing_idx){
                    network_obj->CPT[var][var_cpt_idx] += prob_weight[line];
                }
                else{
                    network_obj->CPT[var][var_cpt_idx] += 1;
                }
            }
        }

        //CPT denominator updation -> 1 pass through CPT
        for (int table_idx=0; table_idx<37; table_idx++){ //ith table in CPT
            vector<int> parents_idx = network_obj->Parents[table_idx]; //index of parents of the variable
            int offset=1; //offset for each value of variable in cpt table
            for (auto idx: parents_idx){
                offset*=network_obj->nValues[idx];
            }
            
            int total_var_val = network_obj->nValues[table_idx];
            for (int idx=0;idx<offset;idx++){ 
                double denom = total_var_val*laplace_const;
                for (int var_val=0;var_val<total_var_val;var_val++){
                    denom+=network_obj->CPT[table_idx][var_val*offset+idx]; //obtaining the required denominator
                }
                for (int var_val=0;var_val<total_var_val;var_val++){
                    network_obj->CPT[table_idx][var_val*offset+idx]/=denom; //dividing the values by required denominator
                }
            }
        }
    }

    void view_data(){
        cout << "q_ind\n" << endl;
        for (const auto &num: q_ind){
            cout << num << " ";
        }

        cout << "\ndata_val\n" << endl;
        for (const auto &vec: data_val){
            for(const auto &val: vec){
                cout << val << " ";
            }
            cout << endl;
        }

        cout << "\ndata_cat\n" << endl;
        for (const auto &vec: data_cat){
            for(const auto &val: vec){
                cout << val << " ";
            }
            cout << endl;
        }
    }
};

int main(int argc,char *argv[]){
    if (argc != 3){
        cerr<<"Missing .bif file and/or datafile"<<endl;
    }
    const int time_limit_seconds = 115;
    auto start_time = std::chrono::high_resolution_clock::now();
    network medical;
    medical.read_network(argv[1]);
    // medical.view_network();

    compute_CPT compute(&medical);
    compute.read_data(argv[2]);
    // compute.view_data();

    compute.intialize_CPT();
    // medical.view_network();
    int j = 0;
    while (true)
    {
        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time);

        if (elapsed_seconds.count() >= time_limit_seconds)
        {
            // Time limit exceeded, exit the loop
            break;
        }
        compute.expectation();
        compute.maximization();
        j++;
        // if(j==20) break;
    }
    cout << "Number of Iterations " << j << endl;
    medical.write_network(argv[1]);
    // medical.view_network();
    
}