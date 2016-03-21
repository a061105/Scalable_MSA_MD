#include <vector>
#include <map>
#include <iostream>
#include <stdlib.h>
#include "util.h"
#include "cache.h"

/** Define interface for a Low-rank matrix
 */

/** matrix of the form: M=UV^T
 */
class LowrankMat{
	public:
	int R;
	int C;	
	LowrankMat(SparseMat* _U, SparseMat* _V, int _D):U(_U),V(_V),D(_D){
		R = U->size();
		C = V->size();
	}

	void sumRows(vector<int>& row_index, double* sums){
		
		double* eTU = new double[D];
		for(int i=0;i<D;i++)
			eTU[i] = 0.0;
		
		//(e^TU)
		for(vector<int>::iterator it=row_index.begin(); it!=row_index.end(); it++){
			int r = *it;
			SparseVec* ui = U->at(r);
			for(SparseVec::iterator it=ui->begin(); it!=ui->end(); it++){
				eTU[ it->first ] += it->second;
			}
		}
		
		//(e^TU)V^T
		for(int j=0;j<C;j++){
			SparseVec* vj = V->at(j);
			double sum = 0.0;
			for(SparseVec::iterator it=vj->begin(); it!=vj->end(); it++){
				sum += eTU[ it->first ] * it->second;
			}
			sums[j] = sum;
		}
		
		delete[] eTU;
	}

	double eval(int i, int j){
		
		return dot(U->at(i), V->at(j));
	}
	
	private:
	int D; // rank
	SparseMat* U;
	SparseMat* V;
};

/** matrix of the form: \sum_i c_i M_i,
 *  where M_i is low-rank matrix.
 */
class CompositeLowrankMat{
	
	public:
	int cache_size;
	ArrayCache* col_cache;
	
	CompositeLowrankMat(int _R, int _C, int _cache_size):R(_R),C(_C),cache_size(_cache_size){

		col_cache = new ArrayCache(cache_size);
	}
	
	void addMat(double c_i, LowrankMat* mat_i){
		
		if( mat_i->R != R || mat_i->C != C ){
			cerr << "mis-match matrices in Composite Lowrank Matrix." << endl;
			exit(0);
		}
		c.push_back(c_i);
		mats.push_back(mat_i);
	}
	
	void getColumns(vector<int>& col_index, vector<double*>& columns){
		
		columns.resize(col_index.size());
		for(int i=0;i<col_index.size();i++){
			
			if( (columns[i]=col_cache->get( col_index[i] )) == NULL ){
				
				columns[i] = new double[R];
				for(int j=0;j<R;j++){
					columns[i][j] = eval( j, col_index[i] );
				}
				col_cache->insert( col_index[i], columns[i] );
			}
		}
	}
	
	
	double* getColumn(int col_index){
		
		double* column;
		if( (column=col_cache->get( col_index )) == NULL ){
			//cerr << "col_index="  << col_index << endl;
			column = new double[R];
			for(int j=0;j<R;j++){
				column[j] = eval( j, col_index );
			}
			col_cache->insert( col_index, column );
		}
		return column;
	}

	void sumRows(vector<int>& row_index, double* sums){
		
		for(int j=0;j<C;j++)
			sums[j] = 0.0;
		
		double* tmp_sums = new double[C];
		for(int i=0;i<mats.size();i++){
			mats[i]->sumRows(row_index, tmp_sums);
			vadd( sums, c[i], tmp_sums, sums, C );
		}

		delete[] tmp_sums;
	}
	
	double eval(int i, int j){
		
		double sum = 0.0;
		for(int k=0;k<mats.size();k++){
			sum += c[k] * mats[k]->eval(i,j);
		}
		return sum;
	}
	
	private:
	int R;
	int C;
	vector<double> c;
	vector<LowrankMat*> mats;
	
	
	/*void cache_column( int col_index, double* col ){
		cache.insert( make_pair(col_index, col) );
	}
	double* get_column_from_cache( int col_index ){
		map<int,double*>::iterator it;
		if( (it=cache.find( col_index )) == cache.end() )
			return NULL;
		else
			return it->second;
	}*/
};
