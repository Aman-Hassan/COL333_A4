#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unordered_map>
#include <iomanip>
#include <algorithm>
#include <vector>
// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

double sampling_constant = 0.0035;

/* In the given starter code accessing a node given the name or its index seems to return and iterator and takes
time since it has to traverse through entire graph list -> why not make it index based?*/

class network
{
public:
	unordered_map<string, int> Name_ind; // Mapping for the name of node to index
	unordered_map<int, string> ind_Name; // Mapping for index to Name of node (incase ever needed)
										 // vector<string> Nodes; //Names of the nodes
	vector<map<string, int>> cat_val;	 // Contains the category-value pairs that a particular node can take -> will be useful while trying to update CPT (for index in CPT)
	vector<vector<string>> Values;		 // Contains the category-value pairs that a particular node can take -> will be useful while trying to update CPT (for index in CPT)
	vector<int> nValues;				 // number of values a node can take -> basically Values[i].size()
	vector<vector<int>> Children;		 // index of all children of a particular node
	vector<vector<int>> Parents;		 // index of all parents of a particular node
	vector<vector<double>> CPT;			 // The entire cpt table ordered according to index of node
	vector<double> probabilities;
	vector<vector<int>> sample; //
	vector<vector<int>> all_possible_data;
	vector<int> missing_idx_list;
	vector<vector<int>> dummy_sample;
	int number_nodes;

	network()
	{
		Name_ind.clear();
		cat_val.resize(37);
		Values.resize(37);	 // Set the outer vector to have 37 elements
		nValues.resize(37);	 // Set the vector of nValues to have 37 elements
		Children.resize(37); // Set the vector of vectors for Children to have 37 elements
		Parents.resize(37);	 // Set the vector of vectors for Parents to have 37 elements
		CPT.resize(37);		 // Set the vector of vectors for CPT to have 37 elements
							 // Now all the vectors will have 37 elements initially
		number_nodes = 0;
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
		number_nodes = index;
		return;
	}

	// * For reading the .dat file
	void read_data(const std::string &filename)
	{
		ifstream file(filename);
		string line;
		vector<int> vals(37); // corersponding values for categories of current line
		// vector<string> cats(37); //corresponding categories of current line
		int find = 0;
		string temp;
		string name;
		map<string, int> category_val;
		if (file.is_open())
		{
			while (!file.eof())
			{
				// vals.clear();
				// vals.resize(37);
				// cats.clear();
				stringstream ss;
				getline(file, line); // Each line in the .dat file
				ss.str(line);
				int index = 0; // store the value where the missing value occurs in a line
				int missing_idx = -1;
				while (ss >> temp)
				{ // Each word in the particular line
					if (temp.compare("\"?\"") == 0)
					{
						vals[index] = -1;
						// cats[index] = "?";
						missing_idx = index;
					}
					else
					{
						// cout<<index<<" "<<vals.size()<<" "<<endl;
						vals[index] = cat_val[index][temp]; // Add the value of the category of a node to vals
															// cats[index] = temp;
															// missing_idx_list.push_back(-1);
					}
					index++;
				}
				missing_idx_list.push_back(missing_idx);
				sample.push_back(vals);
				if (missing_idx != -1)
				{
					vals[missing_idx] = 0;
					dummy_sample.push_back(vals);
					for (int i = 0; i < nValues[missing_idx]; i++)
					{
						vals[missing_idx] = i;
						all_possible_data.push_back(vals);
						// cout<<sample[sample.size()-1][missing_idx]<<endl;
						// cout<<dummy_sample[dummy_sample.size()-1][missing_idx]<<endl;
					}
				}
				else
				{
					all_possible_data.push_back(vals);
				}
				// data_cat.push_back(cats);
			}
		}
		if (find == 1)
			file.close();
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
							solved_alarm << std::fixed << setprecision(4) <<var_cpt[i] << " ";
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

double probab(int idx,int node,int value,network & medical)
{
	vector<int> size_vec;
	vector<int> val_vec;
	val_vec.push_back(value);
	size_vec.push_back(medical.nValues[node]);
	for(int i = 0;i<medical.Parents[node].size();i++)
	{
		val_vec.push_back(medical.sample[idx][medical.Parents[node][i]]);
		size_vec.push_back(medical.nValues[medical.Parents[node][i]]);
	}
	return medical.CPT[node][prob_index(val_vec,size_vec)];
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
			double total = 0.0;
			vector<double> all_poss_prob;
			for (int s = 0; s < N_missing; s++)
			{
				// cout<<missing_idx<<endl;
				double p1 = probab(i,missing_idx,s,medical);
				for (int j = 0; j < medical.Children[missing_idx].size(); j++)
				{
					medical.sample[i][missing_idx] = s;
					double p2 = probab(i,medical.Children[missing_idx][j],medical.sample[i][medical.Children[missing_idx][j]],medical);
					p1 *= p2;
				}
				total += p1 ;
				all_poss_prob.push_back(p1);
			}
			int req_idx = std::distance(all_poss_prob.begin(), std::max_element(all_poss_prob.begin(), all_poss_prob.end()));
			medical.sample[i][missing_idx] = req_idx;
			for (int j = 0; j < all_poss_prob.size(); j++)
			{
				medical.probabilities.push_back(all_poss_prob[j] / total);
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
		medical.CPT[i] = vector<double>(medical.CPT[i].size(),sampling_constant);
		size_vec.push_back(medical.nValues[i]);
		for (int j = 0; j < medical.Parents[i].size(); j++)
		{
			size_vec.push_back(medical.nValues[medical.Parents[i][j]]);
		}
		vector<double> denom(sz, medical.nValues[i]* sampling_constant);
		for (int j = 0; j < medical.all_possible_data.size(); j++)
		{
			val_vec.clear();
			val_vec.push_back(medical.all_possible_data[j][i]);
			for (int k = 0; k < medical.Parents[i].size(); k++)
			{
				val_vec.push_back(medical.all_possible_data[j][medical.Parents[i][k]]);
			}
			int idx = prob_index(val_vec, size_vec);
			// denom[idx % sz] += medical.probabilities[j];
			// numer[idx] += medical.probabilities[j];
			medical.CPT[i][idx] += medical.probabilities[j];
			denom[idx %sz] += medical.probabilities[j];
		}
		double probab;
		for (int j = 0; j < medical.CPT[i].size(); j++)
		{
			if(denom[j%sz] == 0) medical.CPT[i][j] = 0.0;
			else medical.CPT[i][j] /= denom[j%sz];
		}
	}
}

void init_CPT(network &medical)
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
		for (int j = 0; j < medical.dummy_sample.size(); j++)
		{
			val_vec.clear();
			val_vec.push_back(medical.dummy_sample[j][i]);
			for (int k = 0; k < medical.Parents[i].size(); k++)
			{
				val_vec.push_back(medical.dummy_sample[j][medical.Parents[i][k]]);
			}
			int idx = prob_index(val_vec, size_vec);
			
			denom[idx % sz] += 1;
			numer[idx] += 1;
		}
		double probab;
		for (int j = 0; j < medical.CPT[i].size(); j++)
		{
			probab = (numer[j] + sampling_constant) / (denom[j % sz] + sampling_constant * (medical.nValues[i]));
			if (probab ==0)
				probab = 0.0;
			medical.CPT[i][j] = probab;
		}
	}
	// for(int i = 0; i < medical.number_nodes; i++)
	// 	{
	// 		for(int j = 0; j < medical.CPT[i].size(); j++)
	// 		{
	// 			medical.CPT[i][j] = 1.0/medical.nValues[i];
	// 		}
	// 	}
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
	init_CPT(medical);
	int j = 10;
	while (j--)
	{
		expectation(medical);
		cout << j << endl;
		maximization(medical);
	}
	medical.write_network(argv[1]);
	// medical.view_network();
}
