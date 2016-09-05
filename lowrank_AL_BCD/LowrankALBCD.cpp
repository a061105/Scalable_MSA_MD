#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <omp.h>
#include "LowrankALBCD.h"
using namespace std;

double PRIMAL_INF_TOL = 1e-5;
//double PRIMAL_INF_TOL = -1;
double EPS_PHASE_TRANS = 1e-3;
//double EPS_PHASE_TRANS = -1;
double EPS_INNER = 1e-0;
int NUM_REFER  = 10;
int DEFAULT_GEN_SIZE = 10;
int STOP_COUNT = 4;

double dissim_cost( SparseMat2& WT, CompositeLowrankMat* dissim_mat ){
	
	double dist_cost = 0.0;
	double* dissim_col;
	for(SparseMat2::iterator it=WT.begin();it!=WT.end();it++){
		int col_index = it->first;
		SparseVec* Wj = &(it->second);
		dissim_col = dissim_mat->getColumn(col_index);
		for(SparseVec::iterator it=Wj->begin(); it!=Wj->end(); it++){
			dist_cost += dissim_col[it->first] * it->second;
		}
	}
	return dist_cost;
}

double group_norm( SparseMat2& WT, double* lambda ){
	
	double sum = 0.0;
	double max_val;
	for(SparseMat2::iterator it=WT.begin();it!=WT.end();it++){
		SparseVec* Wj = &(it->second);
		max_val = sv_max( Wj );
		sum += lambda[it->first]*max_val;
	}
	return sum;
}

double objective( SparseMat2& WT, CompositeLowrankMat* dissim_mat, double* lambda ){
	
	return dissim_cost(WT, dissim_mat) + group_norm( WT, lambda);
}

/** Primal infeasibility \|W1-1\|_{\infty}
 */
double primal_inf(SparseMat2& WT, int N){
	
	double* W1_1 = new double[N];
	for(int i=0;i<N;i++)
		W1_1[i] = -1.0;
	for(SparseMat2::iterator it=WT.begin();it!=WT.end();it++){
		SparseVec* Wj = &(it->second);
		for(SparseVec::iterator it2=Wj->begin(); it2!=Wj->end(); it2++)
			W1_1[it2->first] += it2->second;
	}
	double inf_max = -1e300;
	double inf;
	for(int i=0;i<N;i++){
		inf = abs(W1_1[i]);
		if( inf > inf_max )
			inf_max = inf;
	}

	delete[] W1_1;
	return inf_max;
}

/** proximal operator of ||.||_{\infty,1} together with W_ij >= 0
 */
void prox( double* Wj, int N, double lambda ){
	
	//prox of non-negative constraints
	for(int i=0;i<N;i++)
		if( Wj[i] < 0.0 )
			Wj[i] = 0.0;
	
	//find break point 
	vector<int> index;
	for(int i=0;i<N;i++)
		index.push_back(i);
	sort( index.begin(), index.end(), ScoreComp(Wj) );
	
	vector<int>::iterator it = index.begin();
	double sum = 0.0;
	int count = 0;
	double prox_avg;
	while(1){
		int i = *it;
		double val = Wj[i];
		sum += val;
		count++;
		prox_avg = (sum-lambda)/count;
		
		it++;
		if( it==index.end() || prox_avg >= Wj[*it] )
			break;
	}
	if( it==index.end() && prox_avg <= 0.0){ // prox results in zero column
		for(int i=0;i<N;i++)
			Wj[i] = 0.0;
		return ;
	}else{
		for(int s=0;s<count;s++) //i in the prox_avg are crunched to val prox_avg
			Wj[ index[s] ] = prox_avg;
		
		//other i remain unchanged
		return ;
	}
}

/** r = W1-1+\alpha^t
 */
int block_RCD( CompositeLowrankMat* DT, int N, double* lambda, vector<int>& act_col_index, double eta, int max_iter, SparseMat2& WT, double* r, int& act_size){
	
	vector<int> index = act_col_index;
	act_size = act_col_index.size();
	
	int iter = 0;
	double PG_1norm, PG_1norm_max; // proximal gradient of column
	double* Wj_new = new double[N];
	double* Wj_change = new double[N];
	double* Dj;
	while( iter < max_iter ){
		
		PG_1norm_max = -1e300;
		random_shuffle( index.begin(), index.end() );
		
		for(int s=0;s<act_size;s++){
			
			int j = index[s];
			double* Dj = DT->getColumn(j);
			//compute quadratic part closed-form solution for column j
			SparseVec* Wj = &(WT[j]);
			sparseToDense(Wj, Wj_new, N);
			for(int i=0;i<N;i++)
				Wj_new[i] = Wj_new[i] - (Dj[i]+eta*r[i])/eta; // W_ij - g_ij / H_ij;
			
			//apply proximal operator of ||.||_{\infty,1} together with W_ij >= 0
			prox( Wj_new, N, lambda[j]/eta );
			
			//compute proximal gradient (for stopping criteria)
			PG_1norm = 0.0;
			for(int i=0;i<N;i++)
				Wj_change[i] = Wj_new[i];
			for(SparseVec::iterator it=Wj->begin(); it!=Wj->end(); it++){
				Wj_change[it->first] -= it->second;
			}
			for(int i=0;i<N;i++)
				PG_1norm += fabs(Wj_change[i]);
			
			if( PG_1norm > PG_1norm_max )
				PG_1norm_max = PG_1norm;
			
			if( PG_1norm < 0.0 ){ //shrinking
				/*act_size--;
				swap(index[s], index[act_size]);
				s--;*/
				continue; //no need to update
			}
			//update Wj and residual
			Wj->clear();
			for(int i=0;i<N;i++){
				if( Wj_new[i] > 1e-3 ){
					Wj->push_back(make_pair(i,Wj_new[i]));
				}
				r[i] += Wj_change[i];
			}
			if( Wj->size() == 0 ){
				WT.erase(j);
			}
		}
		iter++;
		if( iter % 10 == 0 )
			cerr << ".";
		//checking stopping criteria
		if( PG_1norm_max < EPS_INNER ){
			cerr << "*";
			break;
		}
	}

	delete[] Wj_new;
	delete[] Wj_change;
	return iter;
}

void sample_from_residual(double* r, int size, int num_draw, vector<int>& sample_index){
	
	double* tmp = new double[size];
	for(int i=0; i<size; i++){
		if( r[i] < 0.0 )
			tmp[i] = -r[i];
		else
			tmp[i] = 0.0;
	}
	//cumulative
	for(int i=1;i<size;i++){
		tmp[i] = tmp[i-1] + tmp[i];
	}
	double sum = tmp[size-1];
	//sample
	sample_index.clear();
	for(int i=0;i<num_draw;i++){
		double u = sum*((double)rand()/RAND_MAX);
		int j;
		for(j=0; j<size && tmp[j]<u; j++);
		
		sample_index.push_back(j);
	}
	
	delete[] tmp;
}

void gen_greedy_columns(CompositeLowrankMat* dissim_mat, int N, double* r, double eta, double* lambda, int max_gen_size, vector<int>& gen_col_index, SparseMat2& WT){
	
	int* ind = new int[N];
	double* score = new double[N]; //the smaller the better
	for(int i=0;i<N;i++){
		ind[i] = i;
		score[i] = -1e300;
	}
	
	double* trunc_col_grad_sum = new double[N];
	vector<int> ref_index;
	vector<int> row_index;
		
	//sample point according to residual
	sample_from_residual(r, N, NUM_REFER, ref_index);
	//random_uniform(N, NUM_REFER, ref_index);
	
	//random_uniform(N, MAX_GEN_SIZE, gen_col_index);
	//return ;
	
	//compute greedy scores
	for(vector<int>::iterator it=ref_index.begin(); it!=ref_index.end(); it++){
		//compute active row index based on reference point
		int j = *it;
		double* dissim_col = dissim_mat->getColumn(j);
		row_index.clear();
		for(int i=0;i<N;i++){
			if( dissim_col[i] + eta*r[i] < 0.0 )
				row_index.push_back(i);
		}

		//sum active rows to get reference score for each column
		dissim_mat->sumRows( row_index,  trunc_col_grad_sum );
		double sum_r = 0.0;
		for(vector<int>::iterator it=row_index.begin(); it!=row_index.end(); it++)
			sum_r += r[*it];
		for(int j=0;j<N;j++)
			trunc_col_grad_sum[j] += eta*sum_r;
		
		//get min cost for each column over all reference
		for(int j=0;j<N;j++){
			if( -trunc_col_grad_sum[j] > score[j] )
				score[j] = -trunc_col_grad_sum[j];
		}
	}

	//adding exemplar_cost
	for(int j=0;j<N;j++)
		score[j] -= lambda[j];
	
	//sort and retrieve candidates
	sort(ind, ind+N, ScoreComp(score));
	gen_col_index.clear();
	for(int s=0;s<N && gen_col_index.size()<max_gen_size;s++){
		int j = ind[s];
		if( WT.find(j) != WT.end() ) //already in active set
			continue;
		
		gen_col_index.push_back(j);
	}
	
	delete[] ind;
	delete[] score;
	delete[] trunc_col_grad_sum;
}

double AL_BCD( CompositeLowrankMat* dissim_mat, int N, double* lambda, double eta, int max_iter, SparseMat2& WT ){
	
	WT.clear();
	double* r = new double[N]; // r = W1-1+alpha
	for(int i=0;i<N;i++)
		r[i] = -1.0;
	
	vector<int> gen_col_index;
	vector<int> act_col_index;
	
	int max_gen_size = DEFAULT_GEN_SIZE;
	int phase = 0;
	double start = omp_get_wtime();
	int iter = 0;
	int inner_iter;
	int inner_max_iter = 1;
	int active_size = N;
	double p_inf_last = 1e300;
	vector<int> tmp;
	int stable_count = 0;
	int stable_count2 = 0;
	while(iter < max_iter){
		
		//active column generation
		gen_greedy_columns( dissim_mat, N, r, eta, lambda, max_gen_size, gen_col_index, WT);
		for(vector<int>::iterator it=gen_col_index.begin(); it!=gen_col_index.end(); it++)
			act_col_index.push_back(*it);

		//solves subproblem
		//inner_max_iter = (active_size!=0)?(N/active_size):N;
		inner_iter = block_RCD( dissim_mat, N, lambda, act_col_index, eta, inner_max_iter,  WT, r, active_size );
		
		//dual update
		for(int i=0;i<N;i++)
			r[i] += -1.0;
		for(SparseMat2::iterator it=WT.begin(); it!=WT.end(); it++){
			SparseVec* Wj = &(it->second);
			for(SparseVec::iterator it2=Wj->begin(); it2!=Wj->end(); it2++){
				r[it2->first] += it2->second;
			}
		}
		
		//remove inactive columns
		tmp.clear();
		for(vector<int>::iterator it=act_col_index.begin(); it!=act_col_index.end(); it++)
			if( WT.find(*it) != WT.end()  )
				tmp.push_back(*it);
		act_col_index = tmp;

		//dump info
		double obj = objective( WT, dissim_mat, lambda);
		double p_inf = primal_inf( WT, N );
		int nnz = num_nonzeros(WT);
		int nz_col = num_nonzero_rows(WT);
		double end = omp_get_wtime();
		cerr << "iter=" << iter << ", obj=" << obj << ", p_inf=" << p_inf << ", nnz=" << nnz << ", nz_col=" << nz_col << ", #gen_col=" << gen_col_index.size() << ", #inner=" << inner_iter << ", eta=" << eta << ", time=" << end-start <<  ", |cache|=" << dissim_mat->col_cache->size() << endl;
		
		//set column generation size to be no less than half of nnz_col
		max_gen_size = max(DEFAULT_GEN_SIZE, nz_col/3);

		//consider phase transition
		if( phase==0 && p_inf < EPS_PHASE_TRANS ){
			stable_count2++;
		}else{
			stable_count2=0;
		}
		if( stable_count2 >= STOP_COUNT ){
			inner_max_iter = 10000;
			phase = 1;
		}
		
		if( phase==1 ){
			eta *= 2.0;
			for(int i=0;i<N;i++)
				r[i] /= 2.0;
		}
		
		//stropping criteria
		if( p_inf < PRIMAL_INF_TOL ){
			stable_count++;
			if( stable_count >= STOP_COUNT )
				break;
		}else{
			stable_count=0;
		}

		p_inf_last = p_inf;
		iter++;
	}
	
	delete[] r;
	return objective( WT, dissim_mat, lambda );
}

void runALBCD( vector<SparseVec*>& data, double* exemplar_cost, int N, int D, double eta, int max_iter, int cache_size, SparseMat2& WT ){
	
	//truncate lambda_j to zero since y_j (exemplar indicator)=1 when lambda_j <= 0, 
	//which then becomes a constant to W.
	for(int j=0;j<N;j++)
		exemplar_cost[j] = max(exemplar_cost[j],0.0);
	
	//normalize data
	for(int i=0;i<N;i++)
		normalize(data[i]);
	
	//construct low rank dissim matrix
	SparseMat* ones_vec = ones(N,1);
	int rank = 1;
	
	LowrankMat rank_one_11T(ones_vec, ones_vec, rank);
	rank = D;
	LowrankMat sim_mat( &data, &data, rank );
	CompositeLowrankMat dissim_mat(N,N, cache_size);
	
	dissim_mat.addMat(2.0, &rank_one_11T);
	dissim_mat.addMat(-2.0, &sim_mat);
	
	double obj = AL_BCD( &dissim_mat, N, exemplar_cost, eta, max_iter, WT );
	cout << "obj=" << obj << endl;
	cout << "dissim_cost=" << dissim_cost(WT, &dissim_mat) << endl;
}

/** Output: WT (W transpose)
 *  WT row corresponds to exemplar
 *  WT column corredsponds to sample
 **/
void py_ALBCDexemplar(double** X, double* exemplar_cost, int N, int D, double eta, int max_iter, int cache_size, 
		int& W_nnz, int* W_rind, int* W_cind, double* W_val){
	
	vector<SparseVec*> data;
	for(int i=0;i<N;i++){
		SparseVec* sv = new SparseVec();
		for(int j=0;j<D;j++){
			sv->push_back(make_pair(j,X[i][j]));
		}
		data.push_back(sv);
	}
	
	SparseMat2 WT;
	Param* param = new Param();
	
	runALBCD(data, exemplar_cost, N, D, eta, max_iter, cache_size, WT);
	
	//convert WT to W in primary format
	W_nnz = 0;
	for(SparseMat2::iterator it=WT.begin(); it!=WT.end(); it++){
		SparseVec* sv = &(it->second);
		W_nnz += sv->size();
	}
	W_rind = new int[W_nnz];
	W_cind = new int[W_nnz];
	W_val = new double[W_nnz];
	int count=0;
	for(SparseMat2::iterator it=WT.begin(); it!=WT.end(); it++){

		int j = it->first;
		SparseVec* sv = &(it->second);
		for(SparseVec::iterator it2=sv->begin(); it2!=sv->end(); it2++, count++){
			W_rind[count] = it2->first;
			W_cind[count] = j;
			W_val[count] = it2->second;
		}
	}
}

void exit_with_help(){

	cerr << "Usage: LowrankALBCD [data]" << endl;
	cerr << "output: \"clusters\"" << endl;
	cerr << "	-l lambda: parameter controling number of clusters (default #samples/20)" << endl;
	cerr << "	-t tightening: incrementally rounding to integer solution (0 or 1, default 0)" << endl;
	cerr << "	-e eta: Augmented Lagrangian parameter (default 0.01)" << endl;
	cerr << "	-m max_iter: maximum number of iterations (default 10000)" << endl;
	cerr << "	-c cache_size: maximum number of columns of dissimilarity matrix stored (default 10000)" << endl;
	cerr << "	(p.s. Normalized square Euclidean distance d(xi,xj)=2-<xi,xj> is used by default.)" << endl;
	exit(0);
}

void parse_cmd_line(int argc, char** argv, Param* param){

	int i;
	for(i=1;i<argc;i++){
		if( argv[i][0] != '-' )
			break;
		if( ++i >= argc )
			exit_with_help();

		switch(argv[i-1][1]){
			
			case 'l': param->lambda = atof(argv[i]);
				  break;
			case 't': param->do_tightening = atoi(argv[i]);
				  break;
			case 'e': param->eta = atof(argv[i]);
				  break;
			case 'm': param->max_iter = atoi(argv[i]);
				  break;
			case 'c': param->cache_size = atoi(argv[i]);
				  break;
			default:
				  cerr << "unknown option: -" << argv[i-1][1] << endl;
				  exit(0);
		}
	}
	
	if(i>=argc)
		exit_with_help();
	
	param->dataFile = argv[i++];
}

int main(int argc, char** argv){
	
	if( argc < 1+1 )
		exit_with_help();
	
	// parse arguments
	Param* param = new Param();
	parse_cmd_line(argc,argv, param);

	srand(time(NULL));
	//read data (libsvm format)
	int D = -1;
	int N;
	vector<SparseVec*> data;
	vector<int> labels;
	LIBSVMFileParser::parseSVM(param->dataFile, D, data, labels);
	D += 1;//adding bias, feature started at index 1
	N = data.size();
	if( param->lambda < 0 )
		param->lambda = N/20.0;
	cerr << "D=" << D << ", N=" << N << ", lambda=" << param->lambda << endl;
	double* exemplar_cost = new double[N];
	if( param->exemplar_cost_file != NULL ){
		ifstream fin(param->exemplar_cost_file);
		readVec(fin, exemplar_cost, N);
		fin.close();
	}else{
		for(int i=0;i<N;i++)
			exemplar_cost[i] = 0.0;
	}
	for(int i=0;i<N;i++)
		exemplar_cost[i] += param->lambda;
	
	SparseMat2 WT;
	runALBCD(data, exemplar_cost, N, D, param->eta, param->max_iter, param->cache_size, WT);
	
	ofstream fout("clusters");
	printSparseMat(WT, fout);
	fout.close();
}
