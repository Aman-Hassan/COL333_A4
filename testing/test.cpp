#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unordered_map>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

double sampling_constant = 0.0001;

/* In the given starter code accessing a node given the name or its index seems to return and iterator and takes
time since it has to traverse through entire graph list -> why not make it index based?*/

class network
{
public:
	unordered_map<string,int> Name_ind;				 // Mapping for the name of node to index
	unordered_map<int, string> ind_Name; // Mapping for index to Name of node (incase ever needed)
	// vector<string> Nodes; //Names of the nodes
    vector<map<string,int>> cat_val; // Contains the category-value pairs that a particular node can take -> will be useful while trying to update CPT (for index in CPT)
	vector<vector<string>> Values; // Contains the category-value pairs that a particular node can take -> will be useful while trying to update CPT (for index in CPT)
	vector<int> nValues;		   // number of values a node can take -> basically Values[i].size()
	vector<vector<int>> Children;  // index of all children of a particular node
	vector<vector<int>> Parents;   // index of all parents of a particular node
	vector<vector<double>> CPT;	   // The entire cpt table ordered according to index of node
	vector<double> probabilities;
	vector<vector<int>> sample; //
	vector<vector<int>> all_possible_data;
	vector<int> missing_idx_list;

	network()
	{
		Name_ind.clear();
		Values.resize(37);	 // Set the outer vector to have 37 elements
		nValues.resize(37);	 // Set the vector of nValues to have 37 elements
		Children.resize(37); // Set the vector of vectors for Children to have 37 elements
		Parents.resize(37);	 // Set the vector of vectors for Parents to have 37 elements
		CPT.resize(37);		 // Set the vector of vectors for CPT to have 37 elements
							 // Now all the vectors will have 37 elements initially
	}

	// * reads the .bif file and intitializes Name_ind, parents, children, cpt
	void read_network(const std::string &filename)
	{
		ifstream file(filename);
		string line;
		int find = 0;
		string temp;
		string name;
		map<string, int> category_val;
		vector<string> vals;
		int index = 0;
		if (file.is_open())
		{
			while (!file.eof())
			{
				stringstream ss;
				getline(file, line);

				ss.str(line);
				ss >> temp;

				// Updating the mapping for variablename to Index (i.e. Name_ind) and The categories/values a variable can take
				if (temp.compare("variable") == 0)
				{

					ss >> name;
					getline(file, line);

					stringstream ss2;
					ss2.str(line);
					for (int i = 0; i < 4; i++)
					{

						ss2 >> temp;
					}
					category_val.clear();
					vals.clear();
					int val = 0;
					while (temp.compare("};") != 0)
					{
						category_val.emplace(temp, val);
						vals.push_back(temp);
						val++;
						ss2 >> temp;
					}
					Name_ind.emplace(name, index);
					ind_Name.emplace(index, name);
                    cat_val[index] = category_val;
					Values[index] = vals;
					nValues[index] = val;
					index += 1;
					// Graph_Node new_node(name,values.size(),values);
					// int pos=Alarm.addNode(new_node);
				}
				else if (temp.compare("probability") == 0)
				{

					ss >> temp;
					ss >> temp;

					// temp now contains the current node we are at whose req index is accessible using map[temp]
					int node_index = Name_ind[temp];
					// list<Graph_Node>::iterator listIt;
					// list<Graph_Node>::iterator listIt1;
					// listIt=Alarm.search_node(temp);
					// int index=Alarm.get_index(temp);

					ss >> temp; // going to next space-delimited string
					// node_values.clear();
					// getting all parents
					vector<int> temp_parents;
					while (temp.compare(")") != 0)
					{
						temp_parents.push_back(Name_ind[temp]); // Adding parent index for current node into a temp vector
						// listIt1=Alarm.search_node(temp);
						// listIt1->add_child(index);
						Children[Name_ind[temp]].push_back(node_index); // Adding current node index as child for the parent
						// values.push_back(temp);
						ss >> temp;
					}
					Parents[node_index] = temp_parents;
					// listIt->set_Parents(values);
					getline(file, line); // table line
					stringstream ss2;

					ss2.str(line);
					ss2 >> temp; // the word "table" (acc to alarm.bif)

					ss2 >> temp; // cpt val one by one

					vector<double> curr_CPT;
					string::size_type sz;
					while (temp.compare(";") != 0)
					{
						curr_CPT.push_back(atof(temp.c_str()));
						ss2 >> temp;
					}
					CPT[node_index] = curr_CPT;
					// listIt->set_CPT(curr_CPT);
				}
				else
				{
				}
			}
			if (find == 1)
				file.close();
		}

		return;
	}

	// * For reading the .dat file
	void read_data(const std::string& filename){
        ifstream file(filename);
        string line;
        vector<int> vals(37); //corersponding values for categories of current line
        // vector<string> cats(37); //corresponding categories of current line
        int find=0;
        string temp;
        string name;
        map<string,int> category_val;
        if (file.is_open())
        {
            while (! file.eof() )
            {
                vals.clear();
                // cats.clear();
                stringstream ss;
                getline (file,line); // Each line in the .dat file
                ss.str(line);
                int index = 0; //store the value where the missing value occurs in a line
                int missing_idx = -1;
                while(ss>>temp){ // Each word in the particular line
                    if (temp.compare("\"?\"")==0){
                        vals[index] = -1;
                        // cats[index] = "?";
                        missing_idx = index;
                    }
                    else{
                        vals[index] = cat_val[index][temp]; //Add the value of the category of a node to vals
                        // cats[index] = temp;
                        // missing_idx_list.push_back(-1);
                    }
                    index ++;
                }
                missing_idx_list.push_back(missing_idx);
                all_possible_data.push_back(vals);
				if (missing_idx!=-1){
					for (int i=0;i<nValues[missing_idx];i++){
						vals[missing_idx]=i;
						sample.push_back(vals);
					}
				}
				else{
					sample.push_back(vals);
				}
                // data_cat.push_back(cats);
            }
        }
        if(find==1)
        file.close();
    }

	// * For checking the network intialized
	void view_network()
	{
		std::cout << "Name_ind:\n";
		for (const auto &entry : Name_ind)
		{
			std::cout << entry.first << ": " << entry.second << std::endl;
		}

		std::cout << "\nValues:\n";
		for (const auto &vec : Values)
		{
			for (const auto &entry : vec)
			{
				std::cout << entry << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "\nnValues:\n";
		for (const auto &val : nValues)
		{
			std::cout << val << " ";
		}
		std::cout << std::endl;

		std::cout << "\nChildren:\n";
		for (const auto &vec : Children)
		{
			for (const auto &num : vec)
			{
				std::cout << num << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "\nParents:\n";
		for (const auto &vec : Parents)
		{
			for (const auto &num : vec)
			{
				std::cout << num << " ";
			}
			std::cout << std::endl;
		}

		std::cout << "\nCPT:\n";
		for (const auto &vec : CPT)
		{
			for (const auto &num : vec)
			{
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

int prob_index(vector<int> val_vec, vector<int> size_vec)
{
	int a = 0;
	int b = 1;
	for (int i = size_vec.size() - 1; i >= 0; i--)
	{
		a += b * val_vec[i];
		b *= size_vec[i];
	}
	return a;
}

void expectation(network &medical)
{
	medical.probabilities.clear();
	for (int i = 0; i < medical.sample.size(); i++)
	{
		int missing_idx = medical.missing_idx_list[i];
		if (missing_idx == -1)
		{
			medical.probabilities.push_back(1);
		}
		else
		{
			int N_missing = medical.nValues[missing_idx];
			double num = 0.0;
			double den = 0.0;
			vector<double> all_poss_prob;
			for (int s = 0; s < N_missing; s++)
			{
				num = 1.0;
				vector<int> current_sample(medical.sample.begin(), medical.sample.end());
				current_sample[missing_idx] = s;
				vector<int> val_vec, size_vec;
				for (int j = 0; j < medical.Children[missing_idx].size(); j++)
				{
					val_vec.push_back(current_sample[medical.Children[missing_idx][j]]);
					size_vec.push_back(medical.nValues[medical.Children[missing_idx][j]]);
					for (int p = 0; p < medical.Parents[medical.Children[missing_idx][j]].size(); p++)
					{
						val_vec.push_back(current_sample[medical.Parents[medical.Children[missing_idx][j]][p]]);
						size_vec.push_back(medical.nValues[medical.Parents[medical.Children[missing_idx][j]][p]]);
					}
					int a = 0;
					int b = 0;
					num *= medical.CPT[medical.Children[missing_idx][j]][prob_index(val_vec, size_vec)];
				}
				den += num;
				val_vec.clear();
				size_vec.clear();
				val_vec.push_back(s);
				size_vec.push_back(N_missing);
				for (int j = 0; j < medical.Parents[missing_idx].size(); j++)
				{
					val_vec.push_back(current_sample[medical.Parents[missing_idx][j]]);
					size_vec.push_back(medical.nValues[medical.Parents[missing_idx][j]]);
				}
				num *= medical.CPT[missing_idx][prob_index(val_vec, size_vec)];
				all_poss_prob.push_back(num);
			}
			for (int j = 0; j < all_poss_prob.size(); j++)
			{
				medical.probabilities.push_back(all_poss_prob[j] / den);
			}
		}
	}
}

void maximization(network &medical)
{
	for (int i = 0; i < medical.CPT.size(); i++)
	{
		vector<int> val_vec;
		vector<int> size_vec;
		int sz = medical.CPT[i].size() / medical.nValues[i];
		vector<double> denom(sz, 0.0), numer(medical.CPT[i].size(), 0.0);
		size_vec.push_back(medical.nValues[i]);
		for (int j = 0; j < medical.Parents[i].size(); j++)
		{
			size_vec.push_back(medical.nValues[medical.Parents[i][j]]);
		}
		for (int j = 0; j < medical.all_possible_data.size(); j++)
		{
			val_vec.clear();
			val_vec.push_back(i);
			for (int k = 0; k < medical.Parents[i].size(); k++)
			{
				val_vec.push_back(medical.all_possible_data[j][medical.Parents[i][k]]);
			}
			int idx = prob_index(val_vec, size_vec);
			denom[idx % sz] += medical.probabilities[j];
			numer[idx] += medical.probabilities[j];
		}
		double probab;
		for (int j = 0; j < medical.CPT[i].size(); j++)
		{
			probab = (numer[j] + sampling_constant) / (denom[j % sz] + sampling_constant * (medical.nValues[i]));
			if (probab < sampling_constant)
				probab = sampling_constant;
			medical.CPT[i][j] = probab;
		}
	}
}

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		cerr << "Missing .bif file and/or datafile" << endl;
	}
	network medical;
	medical.read_network(argv[1]);
	medical.read_data(argv[2]);
	// medical.view_network();
}
