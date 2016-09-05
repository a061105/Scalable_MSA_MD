#include "LIBSVMFileParser.h"
#include "Lowrank.h"
#include "util.h"

class Param{
	public:
	char* dataFile;
	double lambda;
	double eta;
	int max_iter;
	int cache_size;
	int do_tightening;
	char* exemplar_cost_file;
	
	Param(){
		dataFile = NULL;
		lambda = -1.0;
		eta = 0.01;
		max_iter = 10000;
		cache_size = 10000;
		exemplar_cost_file = NULL;
		do_tightening = 0;
	}
};
