#include <stdlib.h>
#include <iostream>
#include "cache.h"

int main(int argc, char** argv){
	
	ArrayCache* cache = new ArrayCache(atoi(argv[1]));
	
	int index_range = 1000;
	int iter = 1000;
	int size = 1000000;
	int hit = 0;
	int create_count = 0;
	for(int i=0; i<iter; i++){
		
		int ind = rand()%index_range;
		
		//RETRIEVAL
		if( cache->get(ind) == NULL ){
			double* arr = new double[size];
			for(int j=0;j<size;j++){
				arr[j] = 0.0;
			}
			//INSERTION
			cache->insert( ind, arr );
			create_count++;
		}else{
			hit++;
		}
	}
	cerr << "done, hit rate=" << (double)hit/iter << ", create_count=" << create_count << endl;
	cache->dumpInfo(cerr);

	while(1);
}
