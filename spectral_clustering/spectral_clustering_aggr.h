/**************************************************************************
*	 
*
*	 Copyright (c) 2018 Camila P.S. Tautenhain, Mariá C.V. Nascimento
*
*	 This file is part of MOSpecG/SpecG software.
*
* 	 MOSpecG/SpecG software is free software: you can redistribute it and/or modify
*    it under the terms of the GNU Affero General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    MOSpecG/SpecG software is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU Affero General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*****************************************************************************/

#ifndef Spectral_Clustering_H
#define Spectral_Clustering_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string.h>
#include <list>

#include "../defines.h"
#include "util_arrays.h"
#include "flow_matrix.h"
#include "util_eigen.h"


#define DEF_KNUMBER 4
#define DEF_NPOP 5
#define DEF_NPARETO 11

#define DEF_PERC_OFFS 0.4

#define DEF_PERC_EIGEN 0.10//0.25
#define DEF_NGEN 50

#define DEF_IS_NEWMAN TRUE
using namespace std;

class Spectral_Clustering{

private:
	bool directed;

	list<int> l_eig_val_pos;
	list<int> l_eig_val_neg;	

	//Non backtracking matrix
	Flow_Matrix Flow_matrix;

	int max_it;
	int p_leading;
	double alpha;
	double lambda;
	char graph_file[FILE_SIZE];
	char results_folder[FILE_SIZE];
	char time_folder[FILE_SIZE];
	char pareto_file[FILE_SIZE];
	char details_file_eigenvalues[FILE_SIZE];
	char details_file_base_name[FILE_SIZE];
	char details_file_base[FILE_SIZE];
	char details_file_partitions[FILE_SIZE];
	char details_file_eigenvectors[FILE_SIZE];

	int NVert;
	double MEdges;
	int NSampling;
	double MSampling;

	int NPop;
	int NPopOff;
	int ind_offspring;
	int KNumber;
	int PEigen;
	int NGen;
	
	double PercFromOff;
	double PercEigen;

	int PEigen_pos;
	int PEigen_neg;

	double alpha_sc;

	double **adj;
	double *deg;
	double *Bcsc;//Upper triangular
	int *irow; 
	int *pcol;
	int nnz_A;
	int nnz_Bcsc;

	double *edge_csc;//Upper triangular
	int *edge_csc_irow; 
	int *edge_csc_pcol;

	//Move to class population, ind
	double* pop_Q;
	int** pop_part;
	double*** pop_RGroup_pos;
	double*** pop_RGroup_neg;

	double* offs_Q;
	int** offs_part;
	double*** offs_RGroup_pos;
	double*** offs_RGroup_neg;	

	double** ri_vert_pos;		
	double** ri_vert_neg;

	//Pareto 
	double mo_weight_Qin;
	double mo_weight_Qnull;		

	int num_pareto;
	int** pareto_parts;

	//local variable
	int* map_sampl_vert;
	int* map_group_vert;

	int alloc_pareto();

public:
	//Spectral_Clustering();
	Spectral_Clustering(int param_max_it, double param_pLeading, double param_alpha, double param_lambda, char* param_graph_file, char* param_results_folder, char* param_time_folder, char* param_pareto_file, int param_num_pareto);
	~Spectral_Clustering();	
	void initialize_params(int param_max_it, double param_pLeading, double param_alpha, double param_lambda, char* param_graph_file, char* param_results_folder, char* param_time_folder, char* param_pareto_file, int param_num_pareto);

	Spectral_Clustering(int param_max_it, double param_pLeading, double param_alpha, double param_lambda, char* param_graph_file, char* param_results_folder, char* param_time_folder, char* param_pareto_file, int param_num_pareto, int param_n_gen, int param_n_pop, double param_p_offs);
	void initialize_params(int param_max_it, double param_pLeading, double param_alpha, double param_lambda, char* param_graph_file, char* param_results_folder, char* param_time_folder, char* param_pareto_file, int param_num_pareto,int param_n_gen, int param_n_pop, double param_p_offs);

	int set_details_file_base(char* file_details);
	int set_details_file_extra(const char* extra);

	int destroy_data();
	int* select_random_vertices();

	/** Eigenvalues and Eigenvectors **/
	double calc_ri_vall_pos(int l, int i, int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);	
	double calc_ri_vall_neg(int l, int i, int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);	
	double calc_ri_vert_pos(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);	
	double calc_ri_vert_neg(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);	

	//double calc_RGroup(double* RGroup, int* part, int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec);
	double calc_RGroup_pos(int s, double* RGroupS_pos, int* part);
	double calc_RGroup_neg(int s, double* RGroupS_neg, int* part);

	double calc_Q_total(double** RGroup_pos,double** RGroup_neg, int* part);
	double calc_Q_pos(double** RGroup, int* part);
	double calc_Q_neg(double** RGroup, int* part);
	double calc_QModularity_classical(int* part);
	double calc_QModularity(int* part);
	double calc_QModularity(int* part, double *storeObj);
	double calc_QModularity(int* part, int* map_sampl_vert);
	int adjust_positive_PEigen(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);
	int adjust_ortogonal_eig_vec(int* nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);
	int largest_eig_val_eig_vec(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);	
	int index_eig_vec_ri(int lambda,int i);
	void calc_PEigen_pos_neg(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);

	int define_KNumber_square(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec, int NSampling, double multiply);
	int define_KNumber_square(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec, int NSampling);
	/** Population functions **/
	int alloc_pop_arrays();
	int destroy_pop_arrays();

	void print_ind(int ind);
	void print_offs_ind(int ind);
	void print_population();
	void print_population_Q();
	void print_offspring_Q();
	int* best_from_population_Q();

	void copy_rgroup_pop(int ind_offs, int* offs_part, double* offs_Q, double** offs_RGroup_pos, double** offs_RGroup_neg, int ind_dest, int* dest_part, double* dest_Q, double** dest_RGroup_pos, double** dest_RGroup_neg);

	void genetic_algorithm(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec);
	void crossover(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec);//int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec);
	void crossover_offspring(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec, int g);//int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec);
	void mutation(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec, int g);
	void local_search(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec, int g);
	int construct_initial_population(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec);

	double calc_deltaQ_move(int i, double* RGroupSource_pos,double* RGroupSource_neg,  double* RGroupDest_pos, double* RGroupDest_neg);
	double calc_deltaQ_move_penalty(int i, double* RGroupSource_pos,double* RGroupSource_neg,  double* RGroupDest_pos, double* RGroupDest_neg);
	double move_vert(int i,   double* RGroupSource_pos, double* RGroupSource_neg, int* part, int lDest, double* RGroupDest_pos, double* RGroupDest_neg);
	double move_vert(int i, double* RGroupSource_pos, double* RGroupSource_neg, int* part, int lDest, double* RGroupDest_pos, double* RGroupDest_neg, double deltaQ);

	double move_vert_ls(int i,   double* RGroupSource_pos, double* RGroupSource_neg, int* part, int lDest, double* RGroupDest_pos, double* RGroupDest_neg);
	double move_vert_ls(int i, double* RGroupSource_pos, double* RGroupSource_neg, int* part, int lDest, double* RGroupDest_pos, double* RGroupDest_neg, double deltaQ);

	double construct_part_sampling_vert(double** RGroup_pos, double** RGroup_neg, int* part, int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec);
	double construct_part_random(double** RGroup_pos, double** RGroup_neg, int* part, int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec);

	int select_crossover_dest_cluster(double* RGroup_pos_pSourceQ2_ls, double* RGroup_neg_pSourceQ2_ls,double**RGroup_pos_pDestQ1,double** RGroup_neg_pDestQ1);
	int select_crossover_dest_cluster_acos(double* RGroup_pos_pSourceQ2_ls, double* RGroup_neg_pSourceQ2_ls,double**RGroup_pos_pDestQ1,double** RGroup_neg_pDestQ1);

	int update_population_offspring(double perc);

	int** get_pareto_parts();

	/** Read functions **/
	int readpajek_adj(char* graphFile);
	double readpajek_geral(FILE *f, double **B,int *N);	
	int readpajek_geral_csc(char* graphFile);
	int readpajek_geral_sampling_csc(char* graphFile, int* map_sampl_vert);
	int alloc_Bcsc();
	int get_NVert();

	/** Sampling **/
	int sampling_original_vert();//double perc_sampling);
	int sampling_random_vert(double perc_sampling);
	int sampling_modularity(double perc_sampling);
	double recover_sampling_modularity(int *final_part, int* sampl_part, double* storeObj);
	
	/** Spectral Clustering **/
	void print_part_file(char* results_folder, int it, int* part, int size);
	void print_part_file_ensemble(char* results_folder, int* part, int size);
	void print_part_file_mod(char* results_folder, int* part, int size);
	int calc_eigen_values_vectors(double* eig_val_R, double* eig_val_I,	 double* eig_vec,int N_size, int nnz, double* matrix_csc, int* matrix_irow, int* matrix_pcol);
	int** start_spectral_clustering_MO(double perc_sampling);
	int start_spectral_clustering_single(double perc_sampling);	
	int* start_spectral_clustering_ensemble_mod(int k_consensus, double** matrix_m, int mult_matrix_m);
	int** start_spectral_clustering_ensemble_MO(int k_consensus, double** matrix_m, int mult_matrix_m);
	int start_spectral_clustering_matrix(int k_consensus, int nnz, double** matrix_m,double* matrix, int* matrix_irow, int* matrix_pcol);

	/** Details method **/

	void print_details_pop_file(const char* extra);
	void print_details_offspring_file(const char* extra);
	void print_details_part_row(FILE* file, int* part, int size);
	void print_details_crossover(int pSourceQ2,int pDest,int ls,int ld,int ind_offs);
	void print_details_mutation(const char* extra,int ind,int lold,int lnew);

	void print_details_spectral_groups(const char* extra,const char* name_vector, int ind, double** mat, int n1, int n2);
	void print_details_vertex_vectors();
	void print_details_pop_all_group_vectors(const char* extra);
	void print_details_all_group_vectors(const char* extra,int ind);
	void print_details_offspring_all_group_vectors(const char* extra);
	void print_details_offspring_group_vectors(const char* extra,int ind);

};

#endif
