#include <vector>
#include <fstream>
#include "util.h"
using namespace std;


class LIBSVMFileParser {
	public:
		static void  parseSVM(char* fileName, int& numFea, vector<SparseVec*>& data, vector<int>& labels){
			// open input file stream
			ifstream fin(fileName);
			// quit if fail to establish stream
			if( fin.fail() ){
				cerr << "cannot find file." << endl;
				exit(0);
			}
			// initialize local variables
			char* line = new char[MAX_LINE];
			vector<string> tokens;
			numFea=0;
			while( !fin.eof() ){
				fin.getline(line, MAX_LINE);
				string str = string(line);
				tokens = split(str," ");
				// parse label
				labels.push_back( atoi(tokens[0].c_str()) );
				// parse feature
				SparseVec* ins = new SparseVec();
				for(int i=1;i<tokens.size();i++){
					if (tokens[i].find(':') == string::npos) 
						continue;
					vector<string> pv = split(tokens[i],":");
					pair<int,double> pair;
					pair.first = atoi(pv[0].c_str());
					pair.second = atof(pv[1].c_str());
					ins->push_back(pair);
					// cout << "i=" << i << ", "<< pair.first << ":" << pair.second << endl;
				}
				//cerr << "fea="<< ins->fea.back().second << endl;
				//cerr << data->size() << ", " << ins->fea.size() <<  endl;
				if( ins->size() > 0 && ins->back().first > numFea )
					numFea = ins->back().first;
				data.push_back(ins);
			}

			data.pop_back();
			
			delete[] line;
		}

		static vector<string> split(string str, string pattern){
			vector<string> str_split;
			size_t i=0;
			size_t index=0;
			while( index != string::npos ){

				index = str.find(pattern,i);
				str_split.push_back(str.substr(i,index-i));

				i = index+1;
			}
			if ( str_split.back()=="" )
				str_split.pop_back();

			return str_split;
		}
};
