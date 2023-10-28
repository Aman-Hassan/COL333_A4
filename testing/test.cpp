#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <cstdlib>


// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

/* In the given starter code accessing a node given the name or its index seems to return and iterator and takes
time since it has to traverse through entire graph list -> why not make it index based?*/

class network{
public:
    map<string, int> Name_ind; // Mapping for the name of node to index
	unordered_map<int, string> ind_Name; // Mapping for index to Name of node (incase ever needed)
	// vector<string> Nodes; //Names of the nodes
	vector<map<string, int>> Values; // Contains the category-value pairs that a particular node can take -> will be useful while trying to update CPT (for index in CPT)
	vector<int> nValues;			 // number of values a node can take -> basically Values[i].size()
	vector<vector<int>> Children;	 // index of all children of a particular node
	vector<vector<int>> Parents;	 // index of all parents of a particular node
	vector<vector<double>> CPT;		 // The entire cpt table ordered according to index of node
	vector<double> probabilities;
	vector<vector<int>> data; //
	vector<vector<int>> all_possible_data;
	vector<int> missing_idx; 


    network() {
        Name_ind.clear();
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
                        int val = 0;
                        while(temp.compare("};")!=0)
                        {
                            category_val.emplace(temp,val);
                            val++;
                            ss2>>temp;
                        }
                        Name_ind.emplace(name,index);
                        ind_Name.emplace(index,name);
                        Values[index] = category_val;
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
                        // node_values.clear();
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
        
        return;
    }

    // * For checking the network intialized
    void view_network(){
        std::cout << "Name_ind:\n";
        for (const auto &entry : Name_ind) {
            std::cout << entry.first << ": " << entry.second << std::endl;
        }

        std::cout << "\nValues:\n";
        for (const auto &vec : Values) {
            for (const auto &entry : vec) {
                std::cout << " " << entry.first << ": " << entry.second;
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


// TODO: Need to make a class for reading .dat file and having functions to update corresponding CPT and data values
// ? This class should read the .dat file and store it into a vector of vector of strings ( (~10000) lines * 37 variables)
// ? Possible other variables - position of '?' in the data file, ---what else?
// ? Parameters to class would prbly be .dat file and pointer to network class instance 
// ? Class should also have function for expectation, maximisation

int main(int argc,char *argv[]){
    if (argc != 3){
        cerr<<"Missing .bif file and/or datafile"<<endl;
    }
    network medical;
    medical.read_network(argv[1]);
    medical.view_network();

}
