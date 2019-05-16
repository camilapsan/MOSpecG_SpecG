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

#include "spectral_clustering_aggr.h"
/*
Spectral_Clustering::Spectral_Clustering(){

}
*/
Spectral_Clustering::Spectral_Clustering(int param_max_it, double param_pLeading, double param_alpha, double param_lambda, char* param_graph_file, char* param_results_folder, char* param_time_folder, char* param_pareto_file, int param_num_pareto){
	initialize_params(param_max_it, param_pLeading, param_alpha, param_lambda, param_graph_file, param_results_folder, param_time_folder, param_pareto_file, param_num_pareto);
}

Spectral_Clustering::Spectral_Clustering(int param_max_it, double param_pLeading, double param_alpha, double param_lambda, char* param_graph_file, char* param_results_folder, char* param_time_folder, char* param_pareto_file, int param_num_pareto, int param_n_gen, int param_n_pop, double param_p_offs){
	initialize_params(param_max_it, param_pLeading, param_alpha, param_lambda, param_graph_file, param_results_folder, param_time_folder, param_pareto_file, param_num_pareto,param_n_gen, param_n_pop, param_p_offs);

}

Spectral_Clustering::~Spectral_Clustering(){	
	destroy_data();		
}

void Spectral_Clustering::initialize_params(int param_max_it, double param_p_leading, double param_alpha, double param_lambda, char* param_graph_file, char* param_results_folder, char* param_time_folder, char* param_pareto_file, int param_num_pareto){

	max_it = param_max_it;
	alpha = param_alpha;
	lambda = param_lambda;	

	PercEigen=DEF_PERC_EIGEN;

	NVert=0;
	MEdges=0;

	adj=NULL;
	deg=NULL;
	Bcsc=NULL;

	irow=NULL;
	pcol=NULL;
	
	pop_Q=NULL;
	pop_part=NULL;

	ri_vert_pos=NULL;
	pop_RGroup_pos=NULL;
	
	ri_vert_neg=NULL;
	pop_RGroup_neg=NULL;	

	map_sampl_vert=NULL;
	
	edge_csc=NULL;
	edge_csc_irow=NULL;
	edge_csc_pcol=NULL;
	
	KNumber=DEF_KNUMBER;
	
	NPop=DEF_NPOP;
	NPopOff=NPop+1;
	NGen=DEF_NGEN;

	sprintf(graph_file,"%s",param_graph_file);
	sprintf(results_folder,"%s",param_results_folder);
	sprintf(time_folder,"%s",param_time_folder);
	sprintf(pareto_file,"%s",param_pareto_file);

	mo_weight_Qin=1;
	mo_weight_Qnull=1;

	pareto_parts=NULL;

	num_pareto = param_num_pareto;

	srand(time(NULL));
}

void Spectral_Clustering::initialize_params(int param_max_it, double param_p_leading, double param_alpha, double param_lambda, char* param_graph_file, char* param_results_folder, char* param_time_folder, char* param_pareto_file, int param_num_pareto, int param_n_gen, int param_n_pop, double param_p_offs){

	max_it = param_max_it;
	alpha = param_alpha;
	lambda = param_lambda;	
	KNumber=DEF_KNUMBER;
	
	NPop=param_n_pop;
	NPopOff=NPop+1;
	NGen=param_n_gen;
	PercFromOff=param_p_offs;

	PercEigen=param_p_leading;

	NVert=0;
	MEdges=0;

	adj=NULL;
	deg=NULL;
	Bcsc=NULL;

	irow=NULL;
	pcol=NULL;
	
	pop_Q=NULL;
	pop_part=NULL;

	ri_vert_pos=NULL;
	pop_RGroup_pos=NULL;
	
	ri_vert_neg=NULL;
	pop_RGroup_neg=NULL;	

	map_sampl_vert=NULL;
	
	edge_csc=NULL;
	edge_csc_irow=NULL;
	edge_csc_pcol=NULL;
	

	sprintf(graph_file,"%s",param_graph_file);

	sprintf(results_folder,"%s",param_results_folder);
	sprintf(time_folder,"%s",param_time_folder);
	sprintf(pareto_file,"%s",param_pareto_file);

	mo_weight_Qin=1;
	mo_weight_Qnull=1;

	pareto_parts=NULL;

	num_pareto = param_num_pareto;

	srand(time(NULL));
}



int Spectral_Clustering::set_details_file_base(char* file_details){	
	strcpy(details_file_base,file_details);
} 

int Spectral_Clustering::set_details_file_extra(const char* extra){
	sprintf(details_file_base_name,"%s_%s",details_file_base,extra);
	sprintf(details_file_partitions,"%s_partitions_%s.txt",details_file_base,extra);
	sprintf(details_file_eigenvalues,"%s_eigenvalues_%s.txt",details_file_base,extra);	
	sprintf(details_file_eigenvectors,"%s_eigenvectors_%s.txt",details_file_base,extra);	

	#if SAVE_DETAILS == TRUE
		FILE* file=fopen(details_file_partitions,"w");
		fclose(file);
		file=NULL;
	#endif
} 

int Spectral_Clustering::destroy_data(){
	cout << "destroy data 1" << endl;

	///adj can be destroyed before
	adj=Util_Arrays::destroy_2d(adj,NVert,NVert);	
	deg=Util_Arrays::destroy_1d(deg,NVert);

	//can be destroyed after calculating eigenvalues and vectors
	Bcsc=Util_Arrays::destroy_1d(Bcsc,nnz_Bcsc);	

	irow=Util_Arrays::destroy_1d(irow,nnz_Bcsc);	
	pcol=Util_Arrays::destroy_1d(pcol,NSampling+1);
		
	ri_vert_pos= Util_Arrays::destroy_2d(ri_vert_pos, NSampling, PEigen_pos);
	ri_vert_neg= Util_Arrays::destroy_2d(ri_vert_neg, NSampling, PEigen_neg);	

	if(pop_part!=NULL)	pop_part = Util_Arrays::destroy_2d(pop_part, NPopOff, NVert);
	if(pop_Q!=NULL) pop_Q = Util_Arrays::destroy_1d(pop_Q,NPopOff);
	if(pop_RGroup_pos!=NULL) pop_RGroup_pos= Util_Arrays::destroy_3d(pop_RGroup_pos, NPopOff, KNumber, PEigen_pos);
	if(pop_RGroup_neg!=NULL) pop_RGroup_neg= Util_Arrays::destroy_3d(pop_RGroup_neg, NPopOff, KNumber, PEigen_neg);
	
	map_sampl_vert = Util_Arrays::destroy_1d(map_sampl_vert, NSampling);	
	map_group_vert = Util_Arrays::destroy_1d(map_group_vert, NVert);

		
	if(offs_part!=NULL) offs_part = Util_Arrays::destroy_2d(offs_part, NPopOff, NVert);
	if(offs_Q!=NULL) offs_Q = Util_Arrays::destroy_1d(offs_Q,NPopOff);
	if(offs_RGroup_pos!=NULL) offs_RGroup_pos= Util_Arrays::destroy_3d(offs_RGroup_pos, NPopOff, KNumber, PEigen_pos);
	if(offs_RGroup_neg!=NULL) offs_RGroup_neg= Util_Arrays::destroy_3d(offs_RGroup_neg, NPopOff, KNumber, PEigen_neg);

	if(pareto_parts!=NULL){
		pareto_parts = Util_Arrays::destroy_2d(pareto_parts,num_pareto,NVert);
	}

 	return 0;
}

/************************************
 Eigenvalues and Eigenvectors
 ************************************/

int  Spectral_Clustering::index_eig_vec_ri(int lambda,int i){
	return (lambda*NSampling)+i;
}

double Spectral_Clustering::calc_ri_vall_pos(int l, int i, int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec){
	//[r_i]_l = \sqrt{lambda_l}  U_{il}		
	return sqrt(eig_val_R[l])*eig_vec[index_eig_vec_ri(l,i)]; //symetric
	}

double Spectral_Clustering::calc_ri_vall_neg(int l, int i, int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec){
	//[r_i]_l = \sqrt{lambda_l}  U_{il}		
	return sqrt(-eig_val_R[l])*eig_vec[index_eig_vec_ri(l,i)]; //symetric
	}

double Spectral_Clustering::calc_ri_vert_pos(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec){
	int l=0,ind=0;
		
	ri_vert_pos = Util_Arrays::alloc_2d_double(NSampling, PEigen_pos);	
	for(int i=0; i<NSampling;i++){
		//[r_i]_l = \sqrt{lambda_l}  U_{il}		
		ind=0;						
		for (list<int>::iterator it = l_eig_val_pos.begin(); it != l_eig_val_pos.end(); it++,ind++){
			l=*(it);
			assert(eig_val_R[l]>=0);
			ri_vert_pos[i][ind] = calc_ri_vall_pos(l,i,nconv,eig_val_R,eig_val_I,eig_vec);			
		}	
	}
	#if COUT_DEBUG == true
		#endif
}


double Spectral_Clustering::calc_ri_vert_neg(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec){
	int l=0,ind=0;	
	ri_vert_neg = Util_Arrays::alloc_2d_double(NSampling, PEigen_neg);
	for(int i=0; i<NSampling;i++){
		//[r_i]_l = \sqrt{lambda_l}  U_{il}		
		ind=0;
		for (list<int>::iterator it = l_eig_val_neg.begin(); it != l_eig_val_neg.end(); it++,ind++){
			l=(*it);
			assert(eig_val_R[l]<0);
			ri_vert_neg[i][ind] = calc_ri_vall_neg(l,i,nconv,eig_val_R,eig_val_I,eig_vec);			
		}	
	}
	#if COUT_DEBUG == true
	#endif
}

//pop_RGroup ***
//RGroup **
//RGroupS *
//\sum_{i \in s } ri
double Spectral_Clustering::calc_RGroup_pos(int S, double* RGroupS_pos, int* part){
	double contrQ=0;
	//ri is already calculated? yes 

	for(int i=0; i<NSampling;i++){				
		//lambda_0: [0..NVert-1], lambda_1: [NVert..2NVert-1], lambda_2: [2NVert..3NVert-1],...,lambda_l: [l*NVert...(l+1)*NVert-1]
		//real eigenvector of vertex i associated with lambda_l U_{il}			
		if(part[i]==S){
			Util_Arrays::Util_Arrays::add_array(RGroupS_pos, ri_vert_pos[i], PEigen_pos);
		}
	}	
	contrQ = Util_Arrays::norm2_1d(RGroupS_pos,PEigen_pos);
	contrQ = contrQ/(2.0*MEdges);
	return contrQ;
}

double Spectral_Clustering::calc_RGroup_neg(int S, double* RGroupS_neg, int* part){
	double contrQ=0;
	//ri is already calculated? yes 

	for(int i=0; i<NSampling;i++){				
		//lambda_0: [0..NVert-1], lambda_1: [NVert..2NVert-1], lambda_2: [2NVert..3NVert-1],...,lambda_l: [l*NVert...(l+1)*NVert-1]
		//real eigenvector of vertex i associated with lambda_l U_{il}			
		if(part[i]==S){
			Util_Arrays::add_array(RGroupS_neg, ri_vert_neg[i], PEigen_neg);
		}
	}	
	contrQ = Util_Arrays::norm2_1d(RGroupS_neg,PEigen_neg);
	contrQ = contrQ/(2.0*MEdges);
	return contrQ;
}

/*
	Assuming that ri was already calculated!
*/
double Spectral_Clustering::calc_Q_total(double** RGroup_pos, double** RGroup_neg, int* part){
	double valQ=0;
	int s;

	return calc_Q_pos(RGroup_pos,part) - calc_Q_neg(RGroup_neg,part);
}

double Spectral_Clustering::calc_Q_pos(double** RGroup_pos, int* part){
	double valQ=0;
	int s;


	for(s=0;s<KNumber;s++){
		valQ += Util_Arrays::norm2_1d(RGroup_pos[s], PEigen_pos);
	}

	 valQ = valQ/(2.0*MEdges);
	 return valQ;
}

double Spectral_Clustering::calc_Q_neg(double** RGroup_neg, int* part){
	double valQ=0;
	int s;

	for(s=0;s<KNumber;s++){
		valQ += Util_Arrays::norm2_1d(RGroup_neg[s], PEigen_neg);
	}

	valQ = valQ/(2.0*MEdges);
	return valQ;
}

/*
	Original vertices
*/
double Spectral_Clustering::calc_QModularity_classical(int* part){
	double valQ=0;
	int s;

	for(int i=0; i<NVert;i++){				
		for(int j=0; j<NVert;j++){	
			if(part[i]==part[j]){
				valQ += adj[i][j] - (deg[i]*deg[j])/(2.0*MEdges);
			}		
			
		}
	}	

	valQ = valQ/(2.0*MEdges);
	return valQ;
}

double Spectral_Clustering::calc_QModularity(int* part){
	double valQ=0;
	int s;

	for(int i=0; i<NVert;i++){				
		for(int j=0; j<NVert;j++){	
			if(part[i]==part[j]){
				valQ += mo_weight_Qin * adj[i][j] - mo_weight_Qnull * (deg[i]*deg[j])/(2.0*MEdges);
			}
		}
	}	

	valQ = valQ/(2.0*MEdges);
	return valQ;
}


double Spectral_Clustering::calc_QModularity(int* part, double *storeObj){
	double valQ=0;
	int s;
	storeObj[0] =0;
	storeObj[1] =0;

	for(int i=0; i<NVert;i++){				
		for(int j=0; j<NVert;j++){	
			if(part[i]==part[j]){
				storeObj[0] += mo_weight_Qin * adj[i][j];
				storeObj[1] += mo_weight_Qnull * (deg[i]*deg[j])/(2.0*MEdges);
				valQ += mo_weight_Qin * adj[i][j] - mo_weight_Qnull * (deg[i]*deg[j])/(2.0*MEdges);
			}
		}
	}	

	valQ = valQ/(2.0*MEdges);
	storeObj[0]= storeObj[0]/(2.0*MEdges);
	storeObj[1]= storeObj[1]/(2.0*MEdges);
	return valQ;
}


/*
	Considering sampling
*/
double Spectral_Clustering::calc_QModularity(int* part, int* map_sampl_vert) {
	double valQ=0;
	int s, sampl_i, sampl_j;

	for(int i=0; i<NSampling; i++){				
		sampl_i=map_sampl_vert[i];
		for(int j=0; j<NSampling; j++){
			sampl_j=map_sampl_vert[j];
			if(part[sampl_i]==part[sampl_j]){
				valQ += mo_weight_Qin * adj[sampl_i][sampl_j] - mo_weight_Qnull* (deg[sampl_i]*deg[sampl_j])/(2.0*MEdges);
			}		
		}
	}

	valQ = valQ/(2.0*MEdges);
	return valQ;
}

int Spectral_Clustering::adjust_positive_PEigen(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec){
	int countPos=0;
	for (int i=0; i<nconv; i++) {
		if(eig_val_R[i]>0){
			countPos++;
		}
	}
	return countPos;
}

/*eig_val_R is ordered??*/
void Spectral_Clustering::calc_PEigen_pos_neg(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec){	
	PEigen_pos=0;
	PEigen_neg=0;

	for(int l=0;l<PEigen;l++){
		if(eig_val_R[l] >= 0){			
			PEigen_pos+=1;
			l_eig_val_pos.push_back(l);
		}else{			
			PEigen_neg+=1;
			l_eig_val_neg.push_back(l);
		}
	}

	#if COUT_DEBUG == true
		cout << "PEigen=" << PEigen << ", pos=" << PEigen_pos << ", neg="  << PEigen_neg << endl;
	#endif
}


int Spectral_Clustering::adjust_ortogonal_eig_vec(int *nconv, double* eig_val_R, double* eig_val_I, double* eig_vec){
	int newLast=*nconv;	
	int indexVec=0;
	//new last eigenvalue must be 0 associated with the uniform vector 1=(1,1,...,1)		
	PEigen+=1;
	*nconv=*nconv+=1;
	eig_val_R[newLast]=0;
	indexVec=index_eig_vec_ri(newLast, 0);
  	
  	//cout << "..P Eigenvalues:" << endl;
  	for (int i=0; i<NVert; i++,indexVec++){
  		eig_vec[indexVec]=1;
  	}	
}


int Spectral_Clustering::largest_eig_val_eig_vec(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec){
	int indexVal=0,l,indexVec=0,minVal=0;
	double min=eig_val_R[0];
	//smaller eigenvalue must be 1 associated with the uniform vector 1=(1,1,...,1)	

	//find the smaller eigenvalue
	for(int l=0;l<nconv;l++){
		if(eig_val_R[l]<min){
			min=eig_val_R[l];
			minVal=l;
		}
	}
	
	cout << "min=" << min << endl;
	eig_val_R[minVal]=0;		

	indexVec=index_eig_vec_ri(minVal, 0);
  	
  	for (int i=0; i<NVert; i++,indexVec++){
  		eig_vec[indexVec]=1;
  	}
}
/************************************
 Population functions + Genetic Algorithm (GA)
 ************************************/

int Spectral_Clustering::alloc_pop_arrays(){	
	pop_part = Util_Arrays::alloc_2d_int(NPopOff, NVert);
	pop_Q = Util_Arrays::alloc_1d_double(NPopOff);
	
	pop_RGroup_pos = Util_Arrays::alloc_3d_double(NPopOff, KNumber,PEigen_pos);
	pop_RGroup_neg = Util_Arrays::alloc_3d_double(NPopOff, KNumber,PEigen_neg);
	
	offs_part = Util_Arrays::alloc_2d_int(NPopOff, NVert);
	offs_Q = Util_Arrays::alloc_1d_double(NPopOff);

	offs_RGroup_pos = Util_Arrays::alloc_3d_double(NPopOff, KNumber,PEigen_pos);
	offs_RGroup_neg = Util_Arrays::alloc_3d_double(NPopOff, KNumber,PEigen_neg);
	
	return 0; //
}


int Spectral_Clustering::destroy_pop_arrays(){	

	if(pop_part!=NULL)	pop_part = Util_Arrays::destroy_2d(pop_part, NPopOff, NVert);
	if(pop_Q!=NULL) pop_Q = Util_Arrays::destroy_1d(pop_Q,NPopOff);
	if(pop_RGroup_pos!=NULL) pop_RGroup_pos= Util_Arrays::destroy_3d(pop_RGroup_pos, NPopOff, KNumber, PEigen_pos);
	if(pop_RGroup_neg!=NULL) pop_RGroup_neg= Util_Arrays::destroy_3d(pop_RGroup_neg, NPopOff, KNumber, PEigen_neg);
		
	if(offs_part!=NULL) offs_part = Util_Arrays::destroy_2d(offs_part, NPopOff, NVert);
	if(offs_Q!=NULL) offs_Q = Util_Arrays::destroy_1d(offs_Q,NPopOff);
	if(offs_RGroup_pos!=NULL) offs_RGroup_pos= Util_Arrays::destroy_3d(offs_RGroup_pos, NPopOff, KNumber, PEigen_pos);
	if(offs_RGroup_neg!=NULL) offs_RGroup_neg= Util_Arrays::destroy_3d(offs_RGroup_neg, NPopOff, KNumber, PEigen_neg);

	return 0; //
}


void Spectral_Clustering::print_ind(int ind){	
	double conf_Q=calc_QModularity(pop_part[ind],map_sampl_vert);

	cout << "pop_part[" << ind << "]: ";
	Util_Arrays::print_1d(pop_part[ind],NSampling);		
	cout << endl << "pop_Q[" << ind << "]=" << pop_Q[ind] << endl;		
	cout << "conf_Q(mod)=" << conf_Q << endl;
	cout << "conf_Q(vector)=" << calc_Q_pos(pop_RGroup_pos[ind],pop_part[ind]) << "-" << calc_Q_neg(pop_RGroup_neg[ind],pop_part[ind]) << "=" << calc_Q_total(pop_RGroup_pos[ind],pop_RGroup_neg[ind],pop_part[ind]) << endl;
}


void Spectral_Clustering::print_offs_ind(int ind){	
	double conf_Q=calc_QModularity(offs_part[ind],map_sampl_vert);

	cout << "offs_part[" << ind << "]: ";
	Util_Arrays::print_1d(offs_part[ind],NSampling);		
	cout << endl << "offs_Q[" << ind << "]=" << offs_Q[ind] << endl;		
	cout << "conf_Q(mod)=" << conf_Q << endl;
	cout << "conf_Q(vector)=" << calc_Q_pos(offs_RGroup_pos[ind],offs_part[ind]) << "-" << calc_Q_neg(offs_RGroup_neg[ind],offs_part[ind]) << "=" << calc_Q_total(offs_RGroup_pos[ind],offs_RGroup_neg[ind],offs_part[ind]) << endl;
}

void Spectral_Clustering::print_population(){
	for(int ind=0;ind<NPop;ind++){
		cout << "...P Population Ind " << ind << endl;
		print_ind(ind);
		cout << "...................." << endl;
		cout << endl;		
	}
}

void Spectral_Clustering::print_details_pop_file(const char* extra){		
	FILE* file;

	file=fopen(details_file_partitions,"a");	

	for(int ind=0;ind<NPop;ind++){
		fprintf(file,"%s; %d;\tQ=%lf;\t",extra,ind,pop_Q[ind]);
		print_details_part_row(file,pop_part[ind],NSampling);		
	}

	fprintf(file,"\n");
	fclose(file);
	file=NULL;
}


void Spectral_Clustering::print_details_offspring_file(const char* extra){		
	FILE* file;

	file=fopen(details_file_partitions,"a");	

	for(int ind=0;ind<NPop;ind++){
		fprintf(file,"%s; %d;\tQ=%lf;\t",extra,ind,offs_Q[ind]);
		print_details_part_row(file,offs_part[ind],NSampling);		
	}

	fprintf(file,"\n");
	fclose(file);
	file=NULL;
}

void Spectral_Clustering::print_details_part_row(FILE* file, int* part, int size){
	
	for(int i=0;i<size;i++){
		fprintf(file,"%d ",part[i]);		
	}

	fprintf(file,"\n");
	
}

void Spectral_Clustering::print_population_Q(){
	int ind=0;
	double conf_Q=0,conf_S=0;
	
	cout << "...................." << endl;

	for(int ind=0;ind<NPop;ind++){		
		cout << "...P Population Ind " << ind << endl;
		cout << "pop_Q[" << ind << "]=" << pop_Q[ind] << endl;		
		conf_Q=calc_QModularity(pop_part[ind]);
		cout << "conf_Q(mod)=" << conf_Q << endl;	
		conf_S=calc_QModularity(pop_part[ind],map_sampl_vert);
		cout << "conf_Q_Sampl(mod)=" << conf_S << endl;	
		cout << "conf_Q(vector)=" << calc_Q_pos(pop_RGroup_pos[ind],pop_part[ind]) << "-" << calc_Q_neg(pop_RGroup_neg[ind],pop_part[ind]) << "=" << calc_Q_total(pop_RGroup_pos[ind],pop_RGroup_neg[ind],pop_part[ind]) << endl;

		cout << "...................." << endl;		
	}
	cout << endl;		
}


void Spectral_Clustering::print_offspring_Q(){
	int ind=0;
	double conf_Q=0,conf_S=0;
	
	cout << "...................." << endl;

	for(int ind=0;ind<NPop;ind++){		
		cout << "...P OFFS-Population Ind " << ind << endl;
		cout << "pop_Q[" << ind << "]=" << offs_Q[ind] << endl;		
		conf_Q=calc_QModularity(offs_part[ind]);
		cout << "conf_Q(mod)=" << conf_Q << endl;	
		conf_S=calc_QModularity(offs_part[ind],map_sampl_vert);
		cout << "conf_Q_Sampl(mod)=" << conf_S << endl;	
		cout << "conf_Q(vector)=" << calc_Q_pos(offs_RGroup_pos[ind], offs_part[ind]) << "-" << calc_Q_neg(offs_RGroup_neg[ind],offs_part[ind]) << "=" << calc_Q_total(offs_RGroup_pos[ind],offs_RGroup_neg[ind],offs_part[ind]) << endl;

		cout << "...................." << endl;		
	}
	cout << endl;		
}

int* Spectral_Clustering::best_from_population_Q(){
	int ind=0,best_ind=-1;
	double conf_Q=0;
	double best_Q=-2;

	cout << "best from population " << endl;
	for(int ind=0;ind<NPop;ind++){
		if(pop_Q[ind] > best_Q){
			best_Q = pop_Q[ind];
			best_ind = ind;
		}
	}
	assert(best_ind>=0);
	 cout << "..N Sampling=" << NSampling << endl;
	#if COUT_DEBUG == true
		cout << "..P best ind = " << best_ind << endl;
	#endif
	return pop_part[best_ind];
}

/*
	Print spectral groups
*/
void Spectral_Clustering::print_details_spectral_groups(const char* extra,const char* name_vector, int ind, double** mat, int n1, int n2){
	char file_name[FILE_SIZE];
	sprintf(file_name,"%s_%s_%s_%d.txt",details_file_base_name,extra,name_vector,ind);	
	Util_Arrays::print_2d_file(file_name, mat, n1, n2);
}

void Spectral_Clustering::print_details_vertex_vectors(){	
	char file_name[FILE_SIZE];
	sprintf(file_name,"%s_%s.txt",details_file_base_name,"rp");	
	Util_Arrays::print_2d_file(file_name, ri_vert_pos, NSampling, PEigen_pos);

	sprintf(file_name,"%s_%s.txt",details_file_base_name,"rn");	
	Util_Arrays::print_2d_file(file_name, ri_vert_neg, NSampling, PEigen_neg);
}

void Spectral_Clustering::print_details_pop_all_group_vectors(const char* extra){	
	for(int ind=0;ind<NPop;ind++){		
		print_details_spectral_groups("Rp",extra,ind,pop_RGroup_pos[ind],KNumber,PEigen_pos);
		print_details_spectral_groups("Rn",extra,ind,pop_RGroup_neg[ind],KNumber,PEigen_neg);		
	}
}

void Spectral_Clustering::print_details_all_group_vectors(const char* extra,int ind){
	print_details_spectral_groups("Rp","initial",ind,pop_RGroup_pos[ind],KNumber,PEigen_pos);
	print_details_spectral_groups("Rn","initial",ind,pop_RGroup_neg[ind],KNumber,PEigen_neg);		
}

void Spectral_Clustering::print_details_offspring_all_group_vectors(const char* extra){	
	for(int ind=0;ind<NPop;ind++){		
		print_details_spectral_groups("Rp",extra,ind,offs_RGroup_pos[ind],KNumber,PEigen_pos);
		print_details_spectral_groups("Rn",extra,ind,offs_RGroup_neg[ind],KNumber,PEigen_neg);		
	}
}

void Spectral_Clustering::print_details_offspring_group_vectors(const char* extra,int ind){
	print_details_spectral_groups("Rp","initial",ind,offs_RGroup_pos[ind],KNumber,PEigen_pos);
	print_details_spectral_groups("Rn","initial",ind,offs_RGroup_neg[ind],KNumber,PEigen_neg);		
}

/*
	Must calculate ri first
*/
void Spectral_Clustering::genetic_algorithm(int nconv, double* eig_val_R, double* eig_val_I, 	double* eig_vec){	
	//popArrays
	#if COUT_DEBUG == true
		// cout << "*** Construct initial population...**" << endl;
	#endif

	construct_initial_population(nconv, eig_val_R, eig_val_I,eig_vec);
	#if SAVE_DETAILS == TRUE 
		print_details_pop_file("initial");
		print_details_pop_all_group_vectors("initial");
		print_details_vertex_vectors();
	#endif

	#if COUT_POP == true
	print_population_Q();
	print_population();
	#endif


	for(int g=0;g<NGen;g++){
		#if COUT_DEBUG == true
			 cout << "*** GA.. generation:" << g << " ...**" << "\t";
		#endif


		#if SAVE_DETAILS==true
			 FILE* file=fopen(details_file_partitions,"a");
			 fprintf(file,"GEN; %d\n",g);
			 fclose(file);
			 file=NULL;
		#endif
		crossover_offspring(nconv, eig_val_R, eig_val_I,eig_vec,g);		
		mutation(nconv, eig_val_R, eig_val_I,eig_vec,g);
		local_search(nconv, eig_val_R, eig_val_I,eig_vec,g);	
		update_population_offspring(PercFromOff);
	}
	
}

/*
Update population with offspring
*/
int Spectral_Clustering::update_population_offspring(double perc){
	//Up to the perc of the highest fitness individuals of offspring replace the lowest fitness individuals of the old population
	int nhighest = floor(perc*(double)(NPop));
	int ind;
	//find the highest nreplace individuals from offspring
    // cout << "*** Update population... **	" << endl;
    double vmax[nhighest+1];
    int indmax[nhighest+1];
    int a,jins;
    for(a=0;a<nhighest;a++){
        vmax[a]=-NVert;	
        indmax[a]=-1;
    }
    
    for(ind=0; ind<NPop; ind++) {        	
    	for(jins=0; offs_Q[ind] < vmax[jins] && jins < nhighest; jins++){                        
         }
         if(jins < nhighest){
             for(a=nhighest; a>=jins; a--){            
                 vmax[a] = vmax[a-1];
                 indmax[a] = indmax[a-1];
             }                        
                
             vmax[jins]=offs_Q[ind];
             indmax[jins]=ind;     
         }               
    }

	// print_population_Q();
    double vmin[nhighest+1];
    int indmin[nhighest+1];    
    for(a=0;a<nhighest;a++){
        vmin[a]=NVert;	
        indmin[a]=-1;
    }

    for(ind=0; ind<NPop; ind++) {        	
    	for(jins=0; pop_Q[ind] > vmin[jins] && jins < nhighest; jins++){                        
         }
              if(jins < nhighest){
             for(a=nhighest; a>=jins; a--){            
                 vmin[a] = vmin[a-1];
                 indmin[a] = indmin[a-1];
             }                        
                
             vmin[jins]=offs_Q[ind];
             indmin[jins]=ind;     
         }               
    }

    int ind_offs, ind_pop;
    for(a=0;a<nhighest;a++){         
         ind_pop=indmin[a];
         ind_offs=indmax[a];

         copy_rgroup_pop(ind_pop, pop_part[ind_pop], pop_Q, pop_RGroup_pos[ind_pop], pop_RGroup_neg[ind_pop], ind_offs, offs_part[ind_offs], offs_Q, offs_RGroup_pos[ind_offs], offs_RGroup_neg[ind_offs]);         
    }
}

/*
Construct initial population
*/
int Spectral_Clustering::construct_initial_population(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec){	

	for(int ind=0;ind<(NPop);ind++){	
		pop_Q[ind]=construct_part_sampling_vert(pop_RGroup_pos[ind],pop_RGroup_neg[ind], pop_part[ind], nconv, eig_val_R, eig_val_I, eig_vec );		
	}

	return 0;
}

/*
	Return: modularity (vector part. problem) of the partition
*/
double Spectral_Clustering::construct_part_random(double** RGroup_pos, double** RGroup_neg, int* part, int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec){
	int indexeig_vec=0;
	double valQ=0;
	int s,v,i,ind=0;

	//each vertex in a random cluster (0...k)
	//Ignore removed vertices

	for(int i=0; i<NSampling;i++){
		v=map_sampl_vert[i];
		part[v]=Util_Arrays::gen_rand_int(0,KNumber-1);				
		s=part[v];
		Util_Arrays::add_array(RGroup_pos[s],ri_vert_pos[i], PEigen_pos);
		Util_Arrays::add_array(RGroup_neg[s],ri_vert_neg[i], PEigen_neg);
	}
	for(s=0;s<KNumber;s++){		
		valQ += Util_Arrays::norm2_1d(RGroup_pos[s], PEigen_pos);
		valQ -= Util_Arrays::norm2_1d(RGroup_neg[s], PEigen_neg);
	}

	valQ=valQ/(2*MEdges);

	return valQ;	
}



/*	Construct initial partition by randomly sampling KNumber vector arrays ri_vert
	Return: modularity (vector part. problem) of the partition
*/
double Spectral_Clustering::construct_part_sampling_vert(double** RGroup_pos, double** RGroup_neg, int* part, int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec){	
	double valQ=0,prod=0,maxProd;
	int l,i,a,s,sbest=0,ind=0,v;
	
	//sort KNumber of vertices and choose their vertex points as the initial partition
	int* shuffleVert = Util_Arrays::alloc_1d_int(NSampling);
	Util_Arrays::shuffle_array_fisher_yates(shuffleVert,NSampling);	

	for(i=0;i<NSampling;i++){	
		v=map_sampl_vert[i];
		part[v]=-1;
	}

	//select the first KNumber vertices in shuffleVert as the initial RGroup
	for(s=0;s<KNumber;s++){		
		i = shuffleVert[s];
		v=map_sampl_vert[i];
		part[v]=s;
		Util_Arrays::add_array(RGroup_pos[s],ri_vert_pos[i], PEigen_pos);
		Util_Arrays::add_array(RGroup_neg[s],ri_vert_neg[i], PEigen_neg);
	}

	for(i=0;i<NSampling;i++){	
		//define partition of i as the closest R
		v=map_sampl_vert[i];
		if(part[v]==-1){
			//find the closest community:
			maxProd=-2.0;
			for(s=0;s<KNumber;s++){
				prod = Util_Arrays::inner_product_1d(RGroup_pos[s],ri_vert_pos[i],PEigen_pos);
				prod -= Util_Arrays::inner_product_1d(RGroup_neg[s],ri_vert_neg[i],PEigen_neg); //ERROR POS 
				if(prod>maxProd){
					maxProd=prod;
					sbest=s;
				}
			}
			//add i to that community
			part[v]=sbest;
			Util_Arrays::add_array(RGroup_pos[sbest],ri_vert_pos[i],PEigen_pos);		
			Util_Arrays::add_array(RGroup_neg[sbest],ri_vert_neg[i],PEigen_neg);
		}
	}	
	shuffleVert = Util_Arrays::destroy_1d(shuffleVert,NSampling);

	//after calculating 
	for(s=0;s<KNumber;s++){
		valQ += Util_Arrays::norm2_1d(RGroup_pos[s],PEigen_pos);
		valQ -= Util_Arrays::norm2_1d(RGroup_neg[s],PEigen_neg);
	}

	valQ=valQ/(2*MEdges);

	return valQ;	
}

double Spectral_Clustering::calc_deltaQ_move_penalty(int i, double* RGroupSource_pos, double* RGroupSource_neg, double* RGroupDest_pos, double* RGroupDest_neg){
	double inProdT=0,inProdS=0,inProdI=0;
	double inProdT_neg=0,inProdS_neg=0,inProdI_neg=0;

	inProdT= Util_Arrays::inner_product_1d(RGroupDest_pos,ri_vert_pos[i],PEigen_pos);
	inProdS= Util_Arrays::inner_product_1d(RGroupSource_pos,ri_vert_pos[i],PEigen_pos);
	inProdI= Util_Arrays::inner_product_1d(ri_vert_pos[i],ri_vert_pos[i],PEigen_pos);

//- 
	inProdT_neg= Util_Arrays::inner_product_1d(RGroupDest_neg,ri_vert_neg[i],PEigen_neg);
	inProdS_neg= Util_Arrays::inner_product_1d(RGroupSource_neg,ri_vert_neg[i],PEigen_neg);
	inProdI_neg= Util_Arrays::inner_product_1d(ri_vert_neg[i],ri_vert_neg[i],PEigen_neg);

	double deltaQ= (((inProdT- inProdS + inProdI) - (inProdT_neg - inProdS_neg + inProdI_neg ))/(double(MEdges)));

	if(Util_Arrays::compare_equal_array(RGroupSource_pos,ri_vert_pos[i],PEigen_pos)==true && Util_Arrays::compare_equal_array(RGroupSource_neg,ri_vert_neg[i],PEigen_neg)==true){
		deltaQ*=0;
		printf("************ DELTA 0 *********** \n");
	}
	return deltaQ;
}

double Spectral_Clustering::calc_deltaQ_move(int i, double* RGroupSource_pos, double* RGroupSource_neg, double* RGroupDest_pos, double* RGroupDest_neg){
	double inProdT=0,inProdS=0,inProdI=0;
	double inProdT_neg=0,inProdS_neg=0,inProdI_neg=0;

	inProdT= Util_Arrays::inner_product_1d(RGroupDest_pos,ri_vert_pos[i],PEigen_pos);
	inProdS= Util_Arrays::inner_product_1d(RGroupSource_pos,ri_vert_pos[i],PEigen_pos);
	inProdI= Util_Arrays::inner_product_1d(ri_vert_pos[i],ri_vert_pos[i],PEigen_pos);

//-
	inProdT_neg= Util_Arrays::inner_product_1d(RGroupDest_neg,ri_vert_neg[i],PEigen_neg);
	inProdS_neg= Util_Arrays::inner_product_1d(RGroupSource_neg,ri_vert_neg[i],PEigen_neg);
	inProdI_neg= Util_Arrays::inner_product_1d(ri_vert_neg[i],ri_vert_neg[i],PEigen_neg);

	return (((inProdT- inProdS + inProdI) - (inProdT_neg - inProdS_neg + inProdI_neg ))/(double(MEdges)));
}

double Spectral_Clustering::move_vert(int i,  double* RGroupSource_pos, double* RGroupSource_neg, int* part, int lDest, double* RGroupDest_pos, double* RGroupDest_neg){
	double deltaQ=0;

	//Check if RGroupSource_pos is different from ri_vert_pos -> if so, do not move, because the group only has vertex i
	deltaQ=calc_deltaQ_move(i,RGroupSource_pos, RGroupSource_neg, RGroupDest_pos, RGroupDest_neg);
	assert(part[map_sampl_vert[i]]!=lDest);

	part[map_sampl_vert[i]] = lDest;

	Util_Arrays::sub_array(RGroupSource_pos,ri_vert_pos[i],PEigen_pos);
	Util_Arrays::add_array(RGroupDest_pos,ri_vert_pos[i],PEigen_pos);

	Util_Arrays::sub_array(RGroupSource_neg,ri_vert_neg[i],PEigen_neg);
	Util_Arrays::add_array(RGroupDest_neg,ri_vert_neg[i],PEigen_neg);
	
	return deltaQ;	
}


double Spectral_Clustering::move_vert_ls(int i,  double* RGroupSource_pos, double* RGroupSource_neg, int* part, int lDest, double* RGroupDest_pos, double* RGroupDest_neg){
	double deltaQ=0;

	//Check if RGroupSource_pos is different from ri_vert_pos -> if so, do not move, because the group only has vertex i
	deltaQ=calc_deltaQ_move(i,RGroupSource_pos, RGroupSource_neg, RGroupDest_pos, RGroupDest_neg);
	assert(part[map_sampl_vert[i]]!=lDest);

	part[map_sampl_vert[i]] = lDest;

	Util_Arrays::sub_array(RGroupSource_pos,ri_vert_pos[i],PEigen_pos);
	Util_Arrays::add_array(RGroupDest_pos,ri_vert_pos[i],PEigen_pos);
	Util_Arrays::sub_array(RGroupSource_neg,ri_vert_neg[i],PEigen_neg);
	Util_Arrays::add_array(RGroupDest_neg,ri_vert_neg[i],PEigen_neg);

	return deltaQ;	
}


double Spectral_Clustering::move_vert(int i, double* RGroupSource_pos, double* RGroupSource_neg, int* part, int lDest, double* RGroupDest_pos, double* RGroupDest_neg, double deltaQ){	
	assert(part[map_sampl_vert[i]]!=lDest);
	part[map_sampl_vert[i]] = lDest;
	Util_Arrays::sub_array(RGroupSource_pos,ri_vert_pos[i],PEigen_pos);
	Util_Arrays::add_array(RGroupDest_pos,ri_vert_pos[i],PEigen_pos);

	Util_Arrays::sub_array(RGroupSource_neg,ri_vert_neg[i],PEigen_neg);
	Util_Arrays::add_array(RGroupDest_neg,ri_vert_neg[i],PEigen_neg);

	return deltaQ;	
}



double Spectral_Clustering::move_vert_ls(int i, double* RGroupSource_pos, double* RGroupSource_neg, int* part, int lDest, double* RGroupDest_pos, double* RGroupDest_neg, double deltaQ){	
	assert(part[map_sampl_vert[i]]!=lDest);
	part[map_sampl_vert[i]] = lDest;
	Util_Arrays::sub_array(RGroupSource_pos,ri_vert_pos[i],PEigen_pos);
	Util_Arrays::add_array(RGroupDest_pos,ri_vert_pos[i],PEigen_pos);

	Util_Arrays::sub_array(RGroupSource_neg,ri_vert_neg[i],PEigen_neg);
	Util_Arrays::add_array(RGroupDest_neg,ri_vert_neg[i],PEigen_neg);

	return deltaQ;	
}

void Spectral_Clustering::local_search(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec, int g){
	int it=0,lold=0,ind=0,i=0,s=0,sBest=0,count_i;
	double prod=0,maxProd=-2.0;
	
	int* shuffle = Util_Arrays::alloc_1d_int(NSampling);

	#if SAVE_DETAILS == TRUE 
		print_details_offspring_file("local_search");
		char extra[30];
		sprintf(extra,"local_search_%d",g);
		if(g==0){
			print_details_offspring_all_group_vectors(extra);
		}
	#endif

	for(ind=0;ind<NPop;ind++){

		#if COUT_POP == true
			cout << "*** Local Search..;**" <<  endl; 
			cout << "offs_Qbefore=" << offs_Q[ind] << endl;
		#endif 

		for(it=0;it<max_it/10;it++){
			Util_Arrays::shuffle_array_fisher_yates(shuffle,NSampling);	//ESWA: random order 	

			for(count_i=0;count_i<NSampling;count_i++){	
				i = shuffle[count_i];
				lold=offs_part[ind][map_sampl_vert[i]];
	
				//define partition of i as the closest R
				maxProd=-2.0;
				//find the closest community:
						
				for(s=0;s<KNumber;s++){
					if(s!=lold){
						prod = calc_deltaQ_move(i,offs_RGroup_pos[ind][lold],offs_RGroup_neg[ind][lold], offs_RGroup_pos[ind][s],offs_RGroup_neg[ind][s]);

						if(prod>maxProd){
							maxProd=prod;
							sBest=s;
						}
					}
				}
				
				if(maxProd>0 && lold != sBest){
					offs_Q[ind] += move_vert(i, offs_RGroup_pos[ind][lold], offs_RGroup_neg[ind][lold], offs_part[ind], sBest, offs_RGroup_pos[ind][sBest], offs_RGroup_neg[ind][sBest], maxProd);		
				}
			}
		}
	
		#if COUT_POP == true
			cout << "..P AFTER LOCAL SEARCH" << endl;
			print_ind(ind);
		#endif
	}

	#if SAVE_DETAILS == TRUE 
		print_details_offspring_file("local_search");		
		sprintf(extra,"local_search_%d",g);
		if(g==0){
			print_details_offspring_all_group_vectors(extra);
		}
	#endif

}

void Spectral_Clustering::crossover(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec){
	int p1, p2, pDestQ1, pSourceQ2,pAux,vi,lold; //lnew
	int ls,ld;
	double sumQ=0;
	double deltaQ;
	int* partDest=NULL;
	double ind_prob[NPop];
	double offset = 0.0;
	double rnd=0.0;
	
	#if COUT_POP == true
		cout << "*** Crossover (one-way)...**" << endl;
	#endif
	
	//select the two with the highest fitness		
	//roulete whell	
	for(int ind=0;ind<NPop;ind++){		
		sumQ += pop_Q[ind] + 1;
		ind_prob[ind] = pop_Q[ind] + 1;
	}
	for(int ind=0;ind<NPop;ind++){
		ind_prob[ind] = ind_prob[ind]/sumQ;

	}

	//pSourceQ2
	rnd = Util_Arrays::gen_rand_double(0,1);
	offset=0;
	for(int ind=0;ind<NPop;ind++){
	    offset += ind_prob[ind];
	    if (rnd < offset) {
	        pSourceQ2 = ind;
	        break;
	    }
	}
	
	for(int ind=0;ind<NPop;ind++){
		if(ind!=pSourceQ2){
			sumQ += pop_Q[ind] + 1;
			ind_prob[ind] = pop_Q[ind] + 1;
		}		
	}

	for(int ind=0;ind<NPop;ind++){
		if(ind!=pSourceQ2){
			ind_prob[ind] = ind_prob[ind]/sumQ;
		}	
	}
	
	rnd = Util_Arrays::gen_rand_double(0,1);
	offset=0;
	for(int ind=0;ind<NPop;ind++){
		if(ind!=pSourceQ2){
			offset += ind_prob[ind];
			if (rnd < offset) {
				pDestQ1 = ind;
		    	break;
			}
		}		
	}

	//pSource must be the parent with highest fitness
	if(pop_Q[pDestQ1] > pop_Q[pSourceQ2]){
		pAux=pDestQ1;
		pDestQ1=pSourceQ2;
		pSourceQ2=pAux;
	}

	//update only if new dest is better than source
	vi=Util_Arrays::gen_rand_int(0,NSampling-1);
	ls=pop_part[pSourceQ2][map_sampl_vert[vi]];	
	ld=select_crossover_dest_cluster(pop_RGroup_pos[pSourceQ2][ls],pop_RGroup_neg[pSourceQ2][ls],pop_RGroup_pos[pDestQ1],pop_RGroup_neg[pDestQ1]);

		#if COUT_POP == true
			cout << "ParentSource=" << pSourceQ2 << " , ParentDest=" << pDestQ1 << " , l[" << vi << "," << map_sampl_vert[vi] << "]=" << pop_part[pDestQ1][map_sampl_vert[vi]] << " to " <<  ls <<","<< ld  << endl; 
		#endif 

		for(int i=0;i<NSampling;i++){		
			if(ls == pop_part[pSourceQ2][map_sampl_vert[i]]){	
				if(ld!=pop_part[pDestQ1][map_sampl_vert[i]]){		
					lold=pop_part[pDestQ1][map_sampl_vert[i]];
					pop_Q[pDestQ1] += move_vert(i, pop_RGroup_pos[pDestQ1][lold], pop_RGroup_neg[pDestQ1][lold], pop_part[pDestQ1], ld, pop_RGroup_pos[pDestQ1][ld],pop_RGroup_neg[pDestQ1][ld]);
				}
			}
		}
		#if COUT_POP == true
			print_ind(pDestQ1);	
		#endif
}

void Spectral_Clustering::copy_rgroup_pop(int ind_offs, int* offs_part, double* offs_Q, double** offs_RGroup_pos, double** offs_RGroup_neg, int ind_dest, int* dest_part, double* dest_Q, double** dest_RGroup_pos, double** dest_RGroup_neg){
	int ind,k,p;
	offs_Q[ind_offs] = dest_Q[ind_dest];
	for(int i=0;i<NVert; i++){
		offs_part[i] = dest_part[i];
	}
	for(k=0; k<KNumber; k++){
		for(p=0; p<PEigen_pos; p++){
			offs_RGroup_pos[k][p] = dest_RGroup_pos[k][p];
		}
		for(p=0; p<PEigen_neg; p++){
			offs_RGroup_neg[k][p] = dest_RGroup_neg[k][p];
		}
	}
}

int Spectral_Clustering::select_crossover_dest_cluster(double* RGroup_pos_pSourceQ2_ls, double* RGroup_neg_pSourceQ2_ls,double** RGroup_pos_pDestQ1,double** RGroup_neg_pDestQ1){
	double cos_pos, cos_neg, val_pos, val_neg, val, val_best;
	int sbest,s;
	
	#if COUT_POP == true
		printf("... Select crossover dest cluster ld\n");
	#endif

	val_best=-2.0;
	sbest=-1;

	//for each cluster in pop_RGroup_pos and pop_RGroup_neg:
	for(s=0;s<KNumber;s++){
		val_pos = Util_Arrays::inner_product_1d(RGroup_pos_pSourceQ2_ls,RGroup_pos_pDestQ1[s],PEigen_pos);				
		val_neg = Util_Arrays::inner_product_1d(RGroup_neg_pSourceQ2_ls,RGroup_neg_pDestQ1[s],PEigen_neg); //ERROR POS ???		

		val = val_pos + val_neg;

		if(sbest==-1 || val > val_best){
			val_best=val;
			sbest=s;
		}
	}

	return sbest;
}

int Spectral_Clustering::select_crossover_dest_cluster_acos(double* RGroup_pos_pSourceQ2_ls, double* RGroup_neg_pSourceQ2_ls,double** RGroup_pos_pDestQ1,double** RGroup_neg_pDestQ1){
	double cos_pos, cos_neg, arc_cos_pos, arc_cos_neg, val, val_best;
	int sbest,s;
	
	
	val_best=-2.0;
	sbest=-1;

	//for each cluster in pop_RGroup_pos and pop_RGroup_neg:
	for(s=0;s<KNumber;s++){
		cos_pos = Util_Arrays::inner_product_1d(RGroup_pos_pSourceQ2_ls,RGroup_pos_pDestQ1[s],PEigen_pos);
		cos_pos = cos_pos / ( sqrt(Util_Arrays::norm2_1d(RGroup_pos_pSourceQ2_ls,PEigen_pos))*sqrt(Util_Arrays::norm2_1d(RGroup_pos_pDestQ1[s],PEigen_pos)) );		
		arc_cos_pos = acos(cos_pos);
		

		cos_neg = Util_Arrays::inner_product_1d(RGroup_neg_pSourceQ2_ls,RGroup_neg_pDestQ1[s],PEigen_neg); //ERROR POS ???
		cos_neg = cos_neg / ( sqrt(Util_Arrays::norm2_1d(RGroup_neg_pSourceQ2_ls,PEigen_neg))*sqrt(Util_Arrays::norm2_1d(RGroup_neg_pDestQ1[s],PEigen_neg)) );		
		arc_cos_neg = acos(cos_neg);		

		val = arc_cos_pos + arc_cos_neg;

		if(sbest==-1 || val < val_best){
			val_best=val;
			sbest=s;
		}
	}

	//calc inner product
	#if COUT_POP == true
		printf("... sbest=%d\n",sbest);
	#endif

	return sbest;
}
void Spectral_Clustering::crossover_offspring(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec, int g){
	int p1, p2, pDestQ1, pSourceQ2,pAux,vi,lold,ls,ld;
	double sumQ=0;
	double deltaQ;
	int* partDest=NULL;
	double ind_prob[NPop];
	double offset = 0.0;
	double rnd=0.0;

	
	//select the two with the highest fitness
	for(int indoffs=0; indoffs<NPop; indoffs++){
		for(int ind=0;ind<NPop;ind++){		
			sumQ += pop_Q[ind] + 1;
			ind_prob[ind] = pop_Q[ind] + 1;
		}
		
		for(int ind=0;ind<NPop;ind++){
			ind_prob[ind] = ind_prob[ind]/sumQ;
		}

		//pSourceQ2		
		rnd = Util_Arrays::gen_rand_double(0,1);
		offset=0;
		for(int ind=0;ind<NPop;ind++){
		    offset += ind_prob[ind];
		    if (rnd < offset) {
		        pSourceQ2 = ind;	 
		        break;
		    }
		}
		
		sumQ=0;
		//remove pSourceQ2 from possibilities		
		for(int ind=0;ind<NPop;ind++){
			if(ind!=pSourceQ2){
				sumQ += pop_Q[ind] + 1;
				ind_prob[ind] = pop_Q[ind] + 1;
			}		
		}
			
		for(int ind=0;ind<NPop;ind++){
			if(ind!=pSourceQ2){
				ind_prob[ind] = ind_prob[ind]/sumQ;				 
			}	
		}
		
		rnd = Util_Arrays::gen_rand_double(0,1);
		offset=0;
		
		for(int ind=0;ind<NPop;ind++){
			if(ind!=pSourceQ2){
				offset += ind_prob[ind];
				if (rnd < offset) {
					pDestQ1 = ind;
			    	break;
				}
			}		
		}

		copy_rgroup_pop(indoffs, offs_part[indoffs], offs_Q, offs_RGroup_pos[indoffs], offs_RGroup_neg[indoffs], pDestQ1, pop_part[pDestQ1], pop_Q, pop_RGroup_pos[pDestQ1], pop_RGroup_neg[pDestQ1]);		
		//source:pSourceQ2, dest:pDestQ1
		//update only if new dest is better than source
		vi=Util_Arrays::gen_rand_int(0,NSampling-1);
		//lnew=pop_part[pSourceQ2][map_sampl_vert[vi]];	

		//TODO ESWA-> FIND ld as the cluster from pDestQ1 whose inner product with ls is the minimum!
		ls=pop_part[pSourceQ2][map_sampl_vert[vi]];	
		ld=select_crossover_dest_cluster(pop_RGroup_pos[pSourceQ2][ls],pop_RGroup_neg[pSourceQ2][ls],pop_RGroup_pos[pDestQ1],pop_RGroup_neg[pDestQ1]);
		
		for(int i=0;i<NSampling;i++){		
			if(ls == pop_part[pSourceQ2][map_sampl_vert[i]]){
				if(ld!=pop_part[pDestQ1][map_sampl_vert[i]]){		
					lold=pop_part[pDestQ1][map_sampl_vert[i]];
					//update RGroup and update pop_Q					
					offs_Q[indoffs] += move_vert(i, offs_RGroup_pos[indoffs][lold], offs_RGroup_neg[indoffs][lold], offs_part[indoffs], ld, offs_RGroup_pos[indoffs][ld], offs_RGroup_neg[indoffs][ld]);
				}
			}
		}

		#if COUT_POP == true
			cout << "offspring" << endl;
			print_ind(indoffs);	
		#endif

		#if SAVE_DETAILS==true
			print_details_crossover(pSourceQ2,pDestQ1,ls,ld,indoffs);
			char extra[20];
			sprintf(extra,"crossover_%d",g);
			if(g==0){
				print_details_pop_all_group_vectors(extra);
			}
		#endif
	}
}

void Spectral_Clustering::print_details_crossover(int pSourceQ2,int pDestQ1,int ls,int ld,int indoffs){
	FILE* file=fopen(details_file_partitions,"a");
	fprintf(file,"crossover; %d;\ttls; %d;\tld; %d\n",indoffs,ls,ld);
	
	fprintf(file,"crossover; %d;\tsource   ;\t%d;\tQ=%lf;\t",indoffs,pSourceQ2,pop_Q[pSourceQ2]);		
	print_details_part_row(file,pop_part[pSourceQ2],NSampling);		

	fprintf(file,"crossover; %d;\tdest     ;\t%d;\tQ=%lf;\t",indoffs,pDestQ1,pop_Q[pDestQ1]);		
	print_details_part_row(file,pop_part[pDestQ1],NSampling);	

	fprintf(file,"crossover; %d;\toffspring;\t%d;\tQ=%lf;\t",indoffs,indoffs,offs_Q[indoffs]);		
	print_details_part_row(file,offs_part[indoffs],NSampling);	

	fprintf(file,"\n");
	fclose(file);
	file=NULL;
}

void Spectral_Clustering::mutation(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec, int g){
	//sort from 1 to n/2 vertices to change their clusters and sort their clusters
    int num_mutation;
	int vi, lnew, lold;
	int ind;
	int* shuffle = Util_Arrays::alloc_1d_int(NSampling);
	
	num_mutation = Util_Arrays::gen_rand_int(1, floor(NSampling/2));
	ind = Util_Arrays::gen_rand_int(0, NPop-1);
	Util_Arrays::shuffle_array_fisher_yates(shuffle,NSampling);	

	#if COUT_POP == true
		cout << "*** Mutation.. in ind=" << ind << ";**" <<  endl;
		cout << "offs_Qbefore=" << offs_Q[ind] << endl;
	#endif 
	
	#if SAVE_DETAILS==true
		print_details_mutation("before",ind,lold,lnew);
	#endif


	//Sort a individu?
	for(int i=0;i<num_mutation;i++){
		vi = shuffle[i];
		//sort a cluster
		lnew = Util_Arrays::gen_rand_int(0, KNumber-1);	
		lold = offs_part[ind][map_sampl_vert[vi]];

		if(lnew!=lold){							
			//update RGroup and update pop_Q			
			offs_Q[ind] += move_vert(vi, offs_RGroup_pos[ind][lold], offs_RGroup_neg[ind][lold], offs_part[ind], lnew, offs_RGroup_pos[ind][lnew], offs_RGroup_neg[ind][lnew]);
		}
	}

	shuffle = Util_Arrays::destroy_1d(shuffle,NSampling);
		
	#if SAVE_DETAILS==true
		print_details_mutation("after",ind,lold,lnew);
		char extra[20];
		sprintf(extra,"mutation_%d",g);
		if(g==0){
			print_details_all_group_vectors(extra,ind);
		}

		FILE* file=fopen(details_file_partitions,"a");
		fprintf(file,"\n");
		fclose(file);
		file=NULL;
	#endif

	#if COUT_POP == true
		cout << "..P After Mutation" << endl;
		print_ind(ind);
	#endif
}


void Spectral_Clustering::print_details_mutation(const char* extra,int ind,int lold,int lnew){
	FILE* file=fopen(details_file_partitions,"a");	
	fprintf(file,"mutation; %d;\t%s;\tQ=%lf;\t",ind,extra,offs_Q[ind]);		
	print_details_part_row(file,offs_part[ind],NSampling);		

	fclose(file);
	file=NULL;
}



/************************************
 Read functions 
 ************************************/

int Spectral_Clustering::readpajek_adj(char* graph_file){
	int i,j,k,jantes,ret=0;
	int nonZero=0;
	char str[20];
	double aux;
	FILE *f;
	FILE *f2;
	
	f=fopen(graph_file,"r");
 	ret=fscanf(f,"%s %d",str,&NVert);

 	if(ret<1){
 		cout << "ERROR" << endl;
 		exit(1);
 	}

	adj = Util_Arrays::alloc_2d_double(NVert,NVert);
	deg = Util_Arrays::alloc_1d_double(NVert);

	while(strcmp(str,"*Edges")!=0&&strcmp(str,"*Arcslist")!=0&&strcmp(str,"*Arcs")!=0&&strcmp(str,"*Edgeslist")!=0 && strcmp(str,"*Edges0")!=0){
	 	ret=fscanf(f,"%s",str);
	}
	if(strcmp(str,"*Edges")==0){
		directed=false;
		#if POND == false
		 	while(fscanf(f,"%d %d",&i,&j)==2){		 			 
	 			adj[i-1][j-1]=1;
	 			deg[i-1]+=1;		 		
	  			adj[j-1][i-1]=1;	  	
	  			deg[j-1]+=1;
	  
		  		MEdges+=1;		  		
		  		nonZero++;		  		
		 	}
	 	#else
	 	 	while(fscanf(f,"%d %d",&i,&j,&aux)==3){	 		
		 		adj[i-1][j-1]=aux;
		 		deg[i-1]+=aux;
				adj[j-1][i-1]=aux;
		  		deg[j-1]+=aux;		  		
		  		MEdges+=aux;
		  		nonZero++;		  		
	 		}
	 	#endif
	}else if(strcmp(str,"*Edges0")==0){
		directed=false;
		f2=fopen("grafos_nao_direcionados/email-Eu-core_paj.paj","w");		
 		fprintf(f2,"*Vertices 1005\n*Edges\n"); 		
		#if POND == false
		 	while(fscanf(f,"%d %d",&i,&j)==2){
		 		fprintf(f2,"%d %d\n", i+1, j+1);	 	

	 			adj[i][j]=1;
	 			deg[i]+=1;		 		
	  			adj[j][i]=1;	  	
	  			deg[j]+=1;
		  
		  		MEdges+=1;
		  		
		  		nonZero++;		  
		 	}
	 	#else
	 	 	while(fscanf(f,"%d %d",&i,&j,&aux)==3){		 		 	
		 		adj[i-1][j-1]=aux;
		 		deg[i-1]+=aux;
				adj[j-1][i-1]=aux;
		  		deg[j-1]+=aux;		  		
		 		
		  		MEdges+=aux;
		  		nonZero++;
	 		}
	 	#endif

		fclose(f2);
	}else if(strcmp(str,"*Arcs")==0){
		directed=true;
		#if POND == false
		 	while(fscanf(f,"%d %d",&i,&j)==2){		 				 		
		 		adj[i-1][j-1]=1;
		 		deg[i-1]+=1;		 				  					  
		  		MEdges+=1;		  		
		  		nonZero++;		  		
		 	}
	 	#else
	 	 	while(fscanf(f,"%d %d",&i,&j,&aux)==3){		 		
		 		adj[i-1][j-1]=aux;
		 		deg[i-1]+=aux;						 		
		  		MEdges+=aux;
		  		nonZero++;		  		
	 		}
	 	#endif	
	}

	fclose(f);
	f=NULL;	
	f2=NULL;
	nnz_A = nonZero;
 	return nonZero; 
}

int Spectral_Clustering::readpajek_geral_csc(char* graph_file){
	int i,j,k,ind=0,ind_edge=0;
	int nonZero=0;
	char str[20];	

	nnz_Bcsc=(NVert*NVert); //full

	Bcsc = Util_Arrays::alloc_1d_double(nnz_Bcsc);
	irow= Util_Arrays::alloc_1d_int(nnz_Bcsc); // Row index of all nonzero elements of A.
  	pcol=Util_Arrays::alloc_1d_int(NVert+1);  // Pointer to the beginning of each column (in irow and A).


  	edge_csc = Util_Arrays::alloc_1d_double(MEdges*2);
	edge_csc_irow= Util_Arrays::alloc_1d_int(MEdges*2); // Row index of all nonzero elements of A.
  	edge_csc_pcol=Util_Arrays::alloc_1d_int(NVert+1);  // Pointer to the beginning of each column (in irow and A).

  	ind=0;
  	ind_edge=0;
	//Bij = Aij - (di dj)/2*m
	for(j=0;j<NVert;j++){
	  	pcol[j]=ind;
	  	edge_csc_pcol[j]=ind_edge;
		
		for(i=0;i<NVert;i++){
			Bcsc[ind]= (double)adj[i][j]- double(deg[i]*deg[j])/(2.0*MEdges) ;
			irow[ind]=i;

  			ind++;	
  	
  			if(adj[i][j]>= 0){
  				edge_csc[ind_edge]=adj[i][j];
  				edge_csc_irow[ind_edge]=i;
  				ind_edge++;
			}
		}
		
	}
	pcol[NVert]=ind;
	
	pcol[(int)(MEdges)]=ind_edge;
}

int Spectral_Clustering::readpajek_geral_sampling_csc(char* graph_file, int* map_sampl_vert){
	int i,j,k,ind=0,ind_edge=0;
	int sampl_i,sampl_j;
	int nonZero=0;
	char str[20];	

  	ind=0;
  	// ind_edge=0;  	
  	MSampling=0;
	//Bij = Aij - (di dj)/2*m
	for(int j=0;j<NSampling;j++){
	  	pcol[j]=ind;
	  	sampl_j=map_sampl_vert[j];	  	
		for(int i=0;i<NSampling;i++){
			sampl_i=map_sampl_vert[i];
			Bcsc[ind]= mo_weight_Qin * (double)adj[sampl_i][sampl_j]- mo_weight_Qnull * double(deg[sampl_i]*deg[sampl_j])/(2.0*MEdges) ;				
			irow[ind]=i;			
	  		ind++;		  		
		}
		
	}

	pcol[NSampling]=ind;
}

double Spectral_Clustering::readpajek_geral(FILE *f, double **B,int *N){
	int i,j,k,jantes,ret;
	char str[20];

	while(strcmp(str,"*Edges")!=0&&strcmp(str,"*Arcslist")!=0&&strcmp(str,"*Arcs")!=0&&strcmp(str,"*Edgeslist")!=0){
	 	ret=fscanf(f,"%s",str);
	 }
	if(strcmp(str,"*Edges")==0||strcmp(str,"*Arcs")==0){
	 	while(!feof(f)){
	  		ret=fscanf(f,"%d %d",&i,&j);
	  		B[i-1][j-1]=1;
	  		B[j-1][i-1]=1;
	 	}
	}

 return 1;
}

int Spectral_Clustering::get_NVert(){
	return NVert;
}

/************************************ 
 Spectral Clustering 
//  ************************************/

int Spectral_Clustering::calc_eigen_values_vectors(double* eig_val_R, double* eig_val_I, double* eig_vec, int N_size,int nnz, double* matrix_csc, int* matrix_irow, int* matrix_pcol){
	int nconv=0;

	cout << "calc eigen" << endl;
  
 	ARluNonSymMatrix<double, double> matrix(N_size, nnz, matrix_csc, matrix_irow, matrix_pcol);

  	// Defining the eigenvalue problem.
  	ARluNonSymStdEig<double> dprob(PEigen, matrix, "LM", 0, 0.0, 0, NULL, true);  	

  	// Finding eigenvalues.
 	nconv=dprob.EigenValVectors(eig_vec, eig_val_R, eig_val_I);

	//last eigenvalue must be 1 associated with the uniform vector 1=(1,1,...,1)	
	adjust_ortogonal_eig_vec(&nconv, eig_val_R, eig_val_I, eig_vec);
	calc_PEigen_pos_neg(nconv, eig_val_R, eig_val_I, eig_vec);

  	return nconv;
}

int Spectral_Clustering::sampling_random_vert(double perc_sampling){
	NSampling = perc_sampling*NVert;
	cout << "..D NVert=" << NVert << endl << endl;
	int* shuffle = Util_Arrays::alloc_1d_int(NVert);
	Util_Arrays::shuffle_array_fisher_yates(shuffle, NVert);	

	map_sampl_vert = Util_Arrays::alloc_1d_int(NSampling);
	//get the first NSampling nodes.
	for(int i=0;i<NSampling;i++){
		map_sampl_vert[i] = shuffle[i];
	}

	shuffle = Util_Arrays::destroy_1d(shuffle,NVert);
	//return NSampling;

	#if COUT_DEBUG == true
		cout << "..D NSampling=" << NSampling << endl << endl;;
	#endif 
}

//non directed
int Spectral_Clustering::sampling_modularity(double perc_sampling){
	double mean_B=0, aux_B=0;	
	int ind,i,j,vrepr, count=0, count_rem=0, mean_deg=0;
		
	map_group_vert = Util_Arrays::alloc_1d_int(NVert);
	Util_Arrays::init_seq_1d(map_group_vert,NVert); //singleton
    NSampling=NVert;
    map_sampl_vert = Util_Arrays::alloc_1d_int(NSampling);
    Util_Arrays::init_seq_1d(map_sampl_vert,NSampling); //singleton

	#if COUT_DEBUG == true
		cout << "..D NSampling=" << NSampling << endl << endl;
	#endif 
}

double Spectral_Clustering::recover_sampling_modularity(int *final_part, int* sampl_part, double* storeObj){			
	int v_orig, v_sampl, v_repr;
	double final_Q;

	for(int i=0;i<NSampling;i++){
		v_orig=map_sampl_vert[i];
		final_part[v_orig] = sampl_part[i];
	}
	cout << endl;
	
	for(int i=0;i<NVert;i++){		
		if(map_group_vert[i]!=i){		
			v_repr=map_group_vert[i];
			
			final_part[i]=final_part[v_repr];
		}		
	}
	
	final_Q = calc_QModularity(final_part, storeObj);
	cout << "KNumber=" << KNumber << endl; 
	cout << "FINAL Q=" << final_Q << endl; 
	cout << "FINAL Q_classical=" << calc_QModularity_classical(final_part) << endl; 	
}

int Spectral_Clustering::sampling_original_vert(){
	NSampling = NVert;	
	map_sampl_vert = Util_Arrays::alloc_1d_int(NSampling);
	//get the first NSampling nodes.
	for(int i=0;i<NSampling;i++){
		map_sampl_vert[i] = i;
	}	
	
	#if COUT_DEBUG == true
		cout << "..D NSampling=" << NSampling << endl << endl;;
	#endif 
}

int Spectral_Clustering::define_KNumber_square(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec, int NSampling){
	int K=0;
	double radius=0;	
	double mean_d=0, mean_d_aux=0;
	double largest_complex=-1000;
	double largest_real=-1000;	

	
  	for (int i=0; i<nconv; i++) {
		if (eig_val_I[i]==0 && eig_val_R[i] > largest_real) { //complex eigenvalues
	    	largest_real=eig_val_R[i];
	    }	   
	}
	
	radius=sqrt(largest_real);
	cout << "..D largest_real=" << largest_real << ";\t radius=" << radius << endl;

	cout << endl;
 	for (int i=0; i<nconv; i++) {
	    	if((eig_val_R[i]) >= radius){
	    		K++;
	    	}	       
	}

	if(K<=1){
		K=2;
	}else{
		//Maximum number, because the methods might not use all groups
		K=1.25*(double)K;
	}
  	 cout << "..D K=" << K << endl;
  	 return K;
}

int Spectral_Clustering::define_KNumber_square(int nconv, double* eig_val_R, double* eig_val_I,	 double* eig_vec, int NSampling, double multiply){
	int K=0;
	double radius=0;	
	double mean_d=0, mean_d_aux=0;
	double largest_complex=-1000;
	double largest_real=-1000;	

  	for (int i=0; i<nconv; i++) {
	  	if (eig_val_I[i]==0 && multiply*eig_val_R[i] > largest_real) { //complex eigenvalues
	      largest_real=multiply*eig_val_R[i];
	    }	   
	}
	
	radius=sqrt(largest_real);

	cout << endl;
 	for (int i=0; i<nconv; i++) {
	  	if((multiply*eig_val_R[i]) >= radius){
	    	K++;
	    } 
	}

	if(K<=1){
		K=2;
	}else{
		//Maximum number, because the methods might not use all groups
		K=(double)K*1.25;
	}
  	 cout << "..D K=" << K << endl;
  	 return K;
}

int Spectral_Clustering::alloc_Bcsc(){
	nnz_Bcsc=(NSampling*NSampling);

	Bcsc = Util_Arrays::alloc_1d_double(nnz_Bcsc);
	irow= Util_Arrays::alloc_1d_int(nnz_Bcsc); // Row index of all nonzero elements of A.
  	pcol=Util_Arrays::alloc_1d_int(NSampling+1);  // Pointer to the beginning of each column (in irow and A).
}

void Spectral_Clustering::print_part_file(char* results_folder, int it, int* part, int size){
	FILE *f;
	char results_aux[FILE_SIZE];
	
	sprintf(results_aux, "%s_P%d.com",results_folder,it);
	f=fopen(results_aux,"w");

	for(int i=0;i<size;i++){
		fprintf(f,"%d\n",part[i]);		
	}

	fclose(f);	
	f=NULL;
}

void Spectral_Clustering::print_part_file_ensemble(char* results_folder, int* part, int size){
	FILE *f;
	char results_aux[FILE_SIZE];
	sprintf(results_aux, "%s_ensemble.com",results_folder);
	f=fopen(results_aux,"w");
 	
	for(int i=0;i<size;i++){
		fprintf(f,"%d\n",part[i]);
	}
	fclose(f);
	f=NULL;
}


void Spectral_Clustering::print_part_file_mod(char* results_folder, int* part, int size){
	FILE *f;
	char results_aux[FILE_SIZE];
	sprintf(results_aux, "%s_mod.com",results_folder);
	f=fopen(results_aux,"w");
 	
	for(int i=0;i<size;i++){
		fprintf(f,"%d\n",part[i]);
	}
	fclose(f);
	f=NULL;
}


int Spectral_Clustering::alloc_pareto(){
	pareto_parts = Util_Arrays::alloc_2d_int(num_pareto,NVert);
}

int** Spectral_Clustering::get_pareto_parts(){
	return pareto_parts;
}

int** Spectral_Clustering::start_spectral_clustering_MO(double perc_sampling){	 	
 	char fname2[40],aux1[10],str[20], eig_folder[300];
 	int i, nconv;
 	double multiply=0;
 	double storeObj[num_pareto][2];
 	double storeTime[num_pareto];
	struct timespec startIt, endIt;
   	double tTotalIt=0, initTime=0;   
   	char results_aux[FILE_SIZE];
	clock_gettime(CLOCK_MONOTONIC, &startIt);
	FILE *f;	

	#if COUT_GENERAL == true
		cout << "*Start Spectral Clustering MO..*" << endl;
		/*read unweighted graph*/ 	
 		cout << ".. file: " << graph_file << endl; 
 	#endif

 	cout << "num_pareto = " << num_pareto << endl;

	readpajek_adj(graph_file);	
	sampling_modularity(perc_sampling);	

	alloc_Bcsc();
	alloc_pareto();


	//##################################################
	//For each combination of weights	
	double step = 1/(double)(num_pareto-1);
	char extra[10];

	clock_gettime(CLOCK_MONOTONIC, &endIt);
  	initTime = (endIt.tv_sec - startIt.tv_sec);
  	initTime += (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0;

	for(int it=0; it<num_pareto; it++){					
   		clock_gettime(CLOCK_MONOTONIC, &startIt);

		mo_weight_Qin=it*(double)step;
		mo_weight_Qnull = 1-mo_weight_Qin;
		sprintf(extra,"P%d",it);
		set_details_file_extra(extra);

		cout << endl << "################################################" << endl;		
		cout << endl << "### Max " << mo_weight_Qin << "*Qin -" << mo_weight_Qnull << "*Qnull" << endl;					
		readpajek_geral_sampling_csc(graph_file, map_sampl_vert);
			
		#if COUT_GENERAL == true
			cout << ".. N=" << NVert << ", M=" << MEdges << ", nnz_Bcsc=" << nnz_Bcsc <<endl << ".. , PEIG=" << PEigen << endl << endl;
		#endif

		PEigen=floor(PercEigen*(double)(NSampling));

	  	double* eig_val_R = Util_Arrays::alloc_1d_double(PEigen+2);   // Real part of the eigenvalues.
	  	double* eig_val_I = Util_Arrays::alloc_1d_double(PEigen+2); // not used //Util_Arrays::alloc_1d_double(KNumber+1);   // Imaginary part of the eigenvalues.
	  	double* eig_vec = Util_Arrays::alloc_1d_double((PEigen+2)*2*NSampling); //times 2 if has imaginary eigenvalues // new double[272];// Util_Arrays::alloc_1d_double(272); 

	  	nconv = calc_eigen_values_vectors(eig_val_R, eig_val_I, eig_vec, NSampling, NSampling*NSampling, Bcsc, irow, pcol); 		///TODO ESWA SAVE EIGENVECTORS AND EIGENVALUES
		
		multiply=1;
		
		KNumber = define_KNumber_square(nconv, eig_val_R, eig_val_I, eig_vec, NSampling,multiply);			

		cout << "KNumber (square) = " << KNumber << endl;	
		//vector partitioning problem

		///GA!
		#if COUT_GENERAL == true	
			cout << "** Calculate vertex vectors ri... **" << endl;
		#endif
		
		//ri is the same for all solutions and all vertices
		calc_ri_vert_pos(nconv, eig_val_R, eig_val_I, eig_vec);
		calc_ri_vert_neg(nconv, eig_val_R, eig_val_I, eig_vec);
		
		#if COUT_GENERAL == true		
			cout << "** Start GA... **" << endl;
		#endif
		
		alloc_pop_arrays();	
	 	genetic_algorithm(nconv, eig_val_R, eig_val_I, eig_vec);
		
		int *sampl_part = best_from_population_Q();
		
		#if COUT_GENERAL == true		
			cout << "** Recover... **" << endl;
		#endif
				
		recover_sampling_modularity(pareto_parts[it], sampl_part, storeObj[it] );					
 	 	
		//clear iteration			
		eig_val_R = Util_Arrays::destroy_1d(eig_val_R,PEigen+2);
	  	eig_val_I = Util_Arrays::destroy_1d(eig_val_I,PEigen+2);
	  	eig_vec = Util_Arrays::destroy_1d(eig_vec,(PEigen+2)*2*NSampling);
	  	l_eig_val_pos.clear();
	  	l_eig_val_neg.clear();	

		clock_gettime(CLOCK_MONOTONIC, &endIt);
	  	tTotalIt = (endIt.tv_sec - startIt.tv_sec) + (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0 + initTime;
		storeTime[it] = tTotalIt;


		print_part_file(results_folder, it, pareto_parts[it],NVert);		

		sprintf(results_aux, "%s_P%d.txt",time_folder,it);
		f=fopen(results_aux,"w"); 	
		fprintf(f,"%lf\n",tTotalIt);		
		fclose(f);
		f=NULL;

		cout << "print part file ok.." << endl;
		destroy_pop_arrays();
	}
		
	f=fopen(pareto_file,"w"); 	
	cout << "Pareto file = " << pareto_file << endl;
	for(int it=0; it<num_pareto; it++){
      fprintf(f,"Pareto%d \t & \t%lf\t & \t%lf\t & \t%lf\t \\\\ \\hline\n", it, storeTime[it], storeObj[it][0], storeObj[it][1]);      	
	}
	fclose(f);
	f=NULL;
	return pareto_parts;
}

int Spectral_Clustering::start_spectral_clustering_single(double perc_sampling){	 	
 	char fname2[40],aux1[10],str[20], eig_folder[300];
 	int i, nconv;
 	double storeObj[num_pareto][2];
 	double storeTime[num_pareto];
	struct timespec startIt, endIt;
   	double tTotalIt=0, initTime=0;
   	char results_aux[FILE_SIZE];
	clock_gettime(CLOCK_MONOTONIC, &startIt);
	FILE *f;	

	#if COUT_GENERAL == true
		cout << "*Start Spectral Clustering Single..*" << endl;
		/*read unweighted graph*/ 	
 		cout << ".. file: " << graph_file << endl; 
 	#endif

	readpajek_adj(graph_file);	
	sampling_modularity(perc_sampling);	
	
	alloc_Bcsc();
	alloc_pareto();
	
	//##################################################
	//For each combination of weights	
	double step = 1/(double)(num_pareto-1);

	clock_gettime(CLOCK_MONOTONIC, &endIt);
  	initTime = (endIt.tv_sec - startIt.tv_sec);
  	initTime += (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0;
  	
	
   	clock_gettime(CLOCK_MONOTONIC, &startIt);

	mo_weight_Qin=1;
	mo_weight_Qnull=1;
	set_details_file_extra("mod");

	cout << endl << "################################################" << endl;		
	cout << endl << "### Max " << mo_weight_Qin << "*Qin -" << mo_weight_Qnull << "*Qnull" << endl;					
	readpajek_geral_sampling_csc(graph_file, map_sampl_vert);
		
	#if COUT_GENERAL == true
		cout << ".. N=" << NVert << ", M=" << MEdges << ", nnz_Bcsc=" << nnz_Bcsc <<endl << ".. , PEIG=" << PercEigen << endl << endl;
	#endif

	PEigen=floor(PercEigen*(double)(NSampling));	
	cout << "PEigen=" << PEigen<< endl;	

  	double* eig_val_R = Util_Arrays::alloc_1d_double(PEigen+2);   // Real part of the eigenvalues.
  	double* eig_val_I = Util_Arrays::alloc_1d_double(PEigen+2); // Imaginary part of the eigenvalues.
  	double* eig_vec = Util_Arrays::alloc_1d_double((PEigen+2)*2*NSampling); //times 2 if has imaginary eigenvalues

  	nconv = calc_eigen_values_vectors(eig_val_R, eig_val_I, eig_vec, NSampling, NSampling*NSampling, Bcsc, irow, pcol);  	
  	
  	#ifdef SAVE_DETAILS
  	Util_Eigen::print_eigen_values_file(nconv, eig_val_R, eig_val_I, eig_vec, NSampling, details_file_eigenvalues, details_file_eigenvectors);	  	
  	#endif

	KNumber = define_KNumber_square(nconv, eig_val_R, eig_val_I, eig_vec, NSampling);			
	cout << "KNumber (square) = " << KNumber << endl;	
	//vector partitioning problem

	///GA!
	#if COUT_GENERAL == true	
		cout << "** Calculate vertex vectors ri... **" << endl;
	#endif
	
	//ri is the same for all solutions and all vertices
	calc_ri_vert_pos(nconv, eig_val_R, eig_val_I, eig_vec);
	calc_ri_vert_neg(nconv, eig_val_R, eig_val_I, eig_vec);
	
	#if COUT_GENERAL == true		
		cout << "** Start GA... **" << endl;cout << "0" << endl;
	#endif

	alloc_pop_arrays();			 
 	genetic_algorithm(nconv, eig_val_R, eig_val_I, eig_vec);

	//interpolation -> move to function!!		
	int *sampl_part = best_from_population_Q();
	recover_sampling_modularity(pareto_parts[0], sampl_part, storeObj[0] );		


	//clear iteration			
	eig_val_R = Util_Arrays::destroy_1d(eig_val_R,PEigen+2);
  	eig_val_I = Util_Arrays::destroy_1d(eig_val_I,PEigen+2);
  	eig_vec = Util_Arrays::destroy_1d(eig_vec,(PEigen+2)*2*NSampling);
  	l_eig_val_pos.clear();
  	l_eig_val_neg.clear();

	clock_gettime(CLOCK_MONOTONIC, &endIt);
  	tTotalIt = (endIt.tv_sec - startIt.tv_sec) + (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0 + initTime;
	storeTime[0] = tTotalIt;

	print_part_file_mod(results_folder, pareto_parts[0],NVert);		

	destroy_pop_arrays();

	sprintf(results_aux, "%s_mod.txt",time_folder);
	f=fopen(results_aux,"w"); 	
	fprintf(f,"%lf\n",tTotalIt);		
	fclose(f);
	f=NULL;

	return 0;
}


int* Spectral_Clustering::start_spectral_clustering_ensemble_mod(int k_consensus, double** matrix_m, int mult_matrix_m){

 	char fname2[40],aux1[10],str[20], eig_folder[300];
 	int i, nconv;
 	double storeObj[num_pareto][2];
 	double storeTime[num_pareto];
	struct timespec startIt, endIt;
   	double tTotalIt=0, initTime=0;
   	char results_aux[FILE_SIZE];
	clock_gettime(CLOCK_MONOTONIC, &startIt);
	FILE *f;	
	double perc_sampling = 50;

	#if COUT_GENERAL == true

	cout << endl << "################################################" << endl;		
	cout << endl << "################################################" << endl;		
	cout << endl << "################################################" << endl;		

		cout << "*Start Spectral Clustering.. Ensemble modularity*" << endl;
		/*read unweighted graph*/ 	
 		cout << ".. file: " << graph_file << endl; 
 	#endif

	readpajek_adj(graph_file);	
	sampling_modularity(perc_sampling);	

	alloc_Bcsc();
	alloc_pareto();

	//##################################################
	//For each combination of weights	
	double step = 1/(double)(num_pareto-1);

	clock_gettime(CLOCK_MONOTONIC, &endIt);
  	initTime = (endIt.tv_sec - startIt.tv_sec);
  	initTime += (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0;
  	
	
   	clock_gettime(CLOCK_MONOTONIC, &startIt);

	mo_weight_Qin=1;
	mo_weight_Qnull=1;
	set_details_file_extra("ens");
	cout << endl << "################################################" << endl;		
	cout << endl << "### Max " << mo_weight_Qin << "*Qin -" << mo_weight_Qnull << "*Qnull" << endl;					
	
	Util_Arrays::add_matrix(adj, matrix_m, NVert, NVert, mult_matrix_m);
	Util_Arrays::add_matrix_lin_to_array(deg, matrix_m, NVert, NVert, mult_matrix_m);
	MEdges=Util_Arrays::add_matrix_lin_to_num(MEdges, matrix_m, NVert, NVert, mult_matrix_m);


	readpajek_geral_sampling_csc(graph_file, map_sampl_vert);
		
	#if COUT_GENERAL == true
		cout << ".. N=" << NVert << ", M=" << MEdges << ", nnz_Bcsc=" << nnz_Bcsc <<endl << ".. , PEIG=" << PEigen << endl << endl;
	#endif

	PEigen=floor(PercEigen*(double)(NSampling));

  	double* eig_val_R = Util_Arrays::alloc_1d_double(PEigen+2);   // Real part of the eigenvalues.
  	double* eig_val_I = Util_Arrays::alloc_1d_double(PEigen+2);// Imaginary part of the eigenvalues.
  	double* eig_vec = Util_Arrays::alloc_1d_double((PEigen+2)*2*NSampling); //times 2 if has imaginary eigenvalues // new double[272];// Util_Arrays::alloc_1d_double(272); 

  	nconv = calc_eigen_values_vectors(eig_val_R, eig_val_I, eig_vec, NSampling, NSampling*NSampling, Bcsc, irow, pcol);  	


	KNumber = define_KNumber_square(nconv, eig_val_R, eig_val_I, eig_vec, NSampling);			
	cout << "KNumber (square) = " << KNumber << endl;	
	//vector partitioning problem

	///GA!
	#if COUT_GENERAL == true	
		cout << "** Calculate vertex vectors ri... **" << endl;
	#endif
	
	//ri is the same for all solutions and all vertices
	calc_ri_vert_pos(nconv, eig_val_R, eig_val_I, eig_vec);//sqrt(eig_val_R[l])*eig_vec[(l*NVert)+i];			
	calc_ri_vert_neg(nconv, eig_val_R, eig_val_I, eig_vec);//sqrt(eig_val_R[l])*eig_vec[(l*NVert)+i];	
	
	#if COUT_GENERAL == true		
		cout << "** Start GA... **" << endl;cout << "0" << endl;
	#endif

	alloc_pop_arrays();			 
 	genetic_algorithm(nconv, eig_val_R, eig_val_I, eig_vec);

	int *sampl_part = best_from_population_Q();
	recover_sampling_modularity(pareto_parts[0], sampl_part, storeObj[0] );		

	//clear iteration			
	eig_val_R = Util_Arrays::destroy_1d(eig_val_R,PEigen+2);
  	eig_val_I = Util_Arrays::destroy_1d(eig_val_I,PEigen+2);
  	eig_vec = Util_Arrays::destroy_1d(eig_vec,(PEigen+2)*2*NSampling);
  	l_eig_val_pos.clear();
  	l_eig_val_neg.clear();

	clock_gettime(CLOCK_MONOTONIC, &endIt);
  	tTotalIt = (endIt.tv_sec - startIt.tv_sec) + (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0 + initTime;
	storeTime[0] = tTotalIt;

	print_part_file_ensemble(results_folder, pareto_parts[0],NVert);		

	destroy_pop_arrays();

	sprintf(results_aux, "%s_ensemble.txt",time_folder);	
	f=fopen(results_aux,"w"); 	
	fprintf(f,"%lf\n",tTotalIt);		
	fclose(f);
	f=NULL;

	return pareto_parts[0];
}




int** Spectral_Clustering::start_spectral_clustering_ensemble_MO(int k_consensus, double** matrix_m, int mult_matrix_m){

 	char fname2[40],aux1[10],str[20], eig_folder[300];
 	int i, nconv;
 	double storeObj[num_pareto][2];
 	double storeTime[num_pareto];
	struct timespec startIt, endIt;
   	double tTotalIt=0, initTime=0;
   	char results_aux[FILE_SIZE];
	clock_gettime(CLOCK_MONOTONIC, &startIt);
	FILE *f;	
	double perc_sampling = 50;

	#if COUT_GENERAL == true
		cout << endl << "################################################" << endl;		
		cout << endl << "################################################" << endl;		
		cout << endl << "################################################" << endl;		

		cout << "*Start Spectral Clustering.. Ensemble ensemble MO*" << endl;
		/*read unweighted graph*/ 	
 		cout << ".. file: " << graph_file << endl; 
 	#endif

	readpajek_adj(graph_file);	
	sampling_modularity(perc_sampling);	

	alloc_Bcsc();
	alloc_pareto();

//##################################################
	//For each combination of weights	
	double step = 1/(double)(num_pareto-1);
	char extra[10];

	clock_gettime(CLOCK_MONOTONIC, &endIt);
  	initTime = (endIt.tv_sec - startIt.tv_sec);
  	initTime += (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0;

	for(int it=0; it<num_pareto; it++){					
   		clock_gettime(CLOCK_MONOTONIC, &startIt);

		mo_weight_Qin=it*(double)step;
		mo_weight_Qnull = 1-mo_weight_Qin;
		sprintf(extra,"P%d",it);

		set_details_file_extra("ens");
		cout << endl << "################################################" << endl;		
		cout << endl << "### Max " << mo_weight_Qin << "*Qin -" << mo_weight_Qnull << "*Qnull" << endl;					
		
		Util_Arrays::add_matrix(adj, matrix_m, NVert, NVert, mult_matrix_m);
		Util_Arrays::add_matrix_lin_to_array(deg, matrix_m, NVert, NVert, mult_matrix_m);
		MEdges=Util_Arrays::add_matrix_lin_to_num(MEdges, matrix_m, NVert, NVert, mult_matrix_m);


		readpajek_geral_sampling_csc(graph_file, map_sampl_vert);
			
		#if COUT_GENERAL == true
			cout << ".. N=" << NVert << ", M=" << MEdges << ", nnz_Bcsc=" << nnz_Bcsc <<endl << ".. , PEIG=" << PEigen << endl << endl;
		#endif
		PEigen=floor(PercEigen*(double)(NSampling));

	  	double* eig_val_R = Util_Arrays::alloc_1d_double(PEigen+2);   // Real part of the eigenvalues.
	  	double* eig_val_I = Util_Arrays::alloc_1d_double(PEigen+2); // not used //Util_Arrays::alloc_1d_double(KNumber+1);   // Imaginary part of the eigenvalues.
	  	double* eig_vec = Util_Arrays::alloc_1d_double((PEigen+2)*2*NSampling); //times 2 if has imaginary eigenvalues // new double[272];// Util_Arrays::alloc_1d_double(272); 

	  	nconv = calc_eigen_values_vectors(eig_val_R, eig_val_I, eig_vec, NSampling, NSampling*NSampling, Bcsc, irow, pcol);  	
		
		KNumber = define_KNumber_square(nconv, eig_val_R, eig_val_I, eig_vec, NSampling);			
		cout << "KNumber (square) = " << KNumber << endl;	
		//vector partitioning problem

		///GA!
		#if COUT_GENERAL == true	
			cout << "** Calculate vertex vectors ri... **" << endl;
		#endif
		
		//ri is the same for all solutions and all vertices
		calc_ri_vert_pos(nconv, eig_val_R, eig_val_I, eig_vec);//sqrt(eig_val_R[l])*eig_vec[(l*NVert)+i];			
		calc_ri_vert_neg(nconv, eig_val_R, eig_val_I, eig_vec);//sqrt(eig_val_R[l])*eig_vec[(l*NVert)+i];	
		
		#if COUT_GENERAL == true		
			cout << "** Start GA... **" << endl;cout << "0" << endl;
		#endif

		alloc_pop_arrays();			 
	 	genetic_algorithm(nconv, eig_val_R, eig_val_I, eig_vec);
	
		int *sampl_part = best_from_population_Q();
		recover_sampling_modularity(pareto_parts[it], sampl_part, storeObj[it] );		

		//clear iteration			
		eig_val_R = Util_Arrays::destroy_1d(eig_val_R,PEigen+2);
	  	eig_val_I = Util_Arrays::destroy_1d(eig_val_I,PEigen+2);
	  	eig_vec = Util_Arrays::destroy_1d(eig_vec,(PEigen+2)*2*NSampling);
	  	l_eig_val_pos.clear();
	  	l_eig_val_neg.clear();

		clock_gettime(CLOCK_MONOTONIC, &endIt);
	  	tTotalIt = (endIt.tv_sec - startIt.tv_sec) + (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0 + initTime;
		storeTime[0] = tTotalIt;

		print_part_file(results_folder, it, pareto_parts[it],NVert);		

		destroy_pop_arrays();

		sprintf(results_aux, "%s_ensemble_P%d.txt",time_folder,it);
		f=fopen(results_aux,"w"); 	
		fprintf(f,"%lf\n",tTotalIt);		
		fclose(f);
		f=NULL;

	}

	f=fopen(pareto_file,"w"); 	
	cout << "Pareto file = " << pareto_file << endl;
	for(int it=0; it<num_pareto; it++){
      fprintf(f,"Pareto%d \t & \t%lf\t & \t%lf\t & \t%lf\t \\\\ \\hline\n", it, storeTime[it], storeObj[it][0], storeObj[it][1]);      	
	}
	fclose(f);
	f=NULL;
	
	return pareto_parts;
}


int Spectral_Clustering::start_spectral_clustering_matrix(int k_consensus, int nnz, double** matrix_m, double* matrix, int* matrix_irow, int* matrix_pcol){	 	
 	char fname2[40],aux1[10],str[20], eig_folder[300];
 	int i, nconv;
 	struct timespec startIt, endIt;
   	double tTotalIt=0, initTime=0;
   	char results_aux[FILE_SIZE];
	clock_gettime(CLOCK_MONOTONIC, &startIt);
	FILE *f;	

	NSampling=NVert;

	#if COUT_GENERAL == true
		cout << "*Start Spectral Clustering..*" << endl;
		/*read unweighted graph*/ 	
 		cout << ".. Matrix, file: " << graph_file << endl; 
 	#endif

	//##################################################
	//For each combination of weights	
	double step = 1/(double)(num_pareto-1);
	mo_weight_Qin=1;
	mo_weight_Qnull=1;
	double sum=0;
	int ind=0;

	clock_gettime(CLOCK_MONOTONIC, &endIt);
  	initTime = (endIt.tv_sec - startIt.tv_sec);
  	initTime += (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0;

	
   	clock_gettime(CLOCK_MONOTONIC, &startIt);
		
	#if COUT_GENERAL == true
		cout << ".. N=" << NVert << ", M=" << MEdges << ", nnz_Bcsc=" << nnz_Bcsc <<endl << ".. , PEIG=" << PEigen << endl << endl;
	#endif
	PEigen=NSampling-2;

  	double* eig_val_R = Util_Arrays::alloc_1d_double(PEigen+2);   // Real part of the eigenvalues.
  	double* eig_val_I = Util_Arrays::alloc_1d_double(PEigen+2); // not used //Util_Arrays::alloc_1d_double(KNumber+1);   // Imaginary part of the eigenvalues.
  	double* eig_vec = Util_Arrays::alloc_1d_double((PEigen+2)*2*NVert); //times 2 if has imaginary eigenvalues // new double[272];// Util_Arrays::alloc_1d_double(272); 

  	nconv = calc_eigen_values_vectors(eig_val_R, eig_val_I, eig_vec, NVert, nnz, matrix, matrix_irow, matrix_pcol); 

	KNumber = k_consensus;
	cout << "KNumber (square) = " << KNumber << endl;	
	//vector partitioning problem

	#if COUT_GENERAL == true	
		cout << "** Calculate vertex vectors ri... **" << endl;
	#endif
	
	//ri is the same for all solutions and all vertices
	calc_ri_vert_pos(nconv, eig_val_R, eig_val_I, eig_vec);//sqrt(eig_val_R[l])*eig_vec[(l*NVert)+i];			
	calc_ri_vert_neg(nconv, eig_val_R, eig_val_I, eig_vec);//sqrt(eig_val_R[l])*eig_vec[(l*NVert)+i];	
	
	#if COUT_GENERAL == true		
		cout << "** Start GA... **" << endl;
	#endif
	alloc_pop_arrays();
	
  	construct_initial_population(nconv, eig_val_R, eig_val_I,eig_vec);
	local_search(nconv, eig_val_R, eig_val_I,eig_vec,0);	
	local_search(nconv, eig_val_R, eig_val_I,eig_vec,0);	
	local_search(nconv, eig_val_R, eig_val_I,eig_vec,0);	
	local_search(nconv, eig_val_R, eig_val_I,eig_vec,0);	
		int *sampl_part = best_from_population_Q();
	
	//clear iteration			
	eig_val_R = Util_Arrays::destroy_1d(eig_val_R,PEigen+2);
  	eig_val_I = Util_Arrays::destroy_1d(eig_val_I,PEigen+2);
  	eig_vec = Util_Arrays::destroy_1d(eig_vec,(PEigen+2)*2*NVert);
  	l_eig_val_pos.clear();
  	l_eig_val_neg.clear();

	clock_gettime(CLOCK_MONOTONIC, &endIt);
  	tTotalIt = (endIt.tv_sec - startIt.tv_sec) + (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0 + initTime;
	
	cout << "KNumber=" << KNumber << endl; 
	cout << "FINAL Q=" << calc_QModularity(sampl_part) << endl; 
	
	mo_weight_Qin=1;
	mo_weight_Qnull=1;
	cout << "FINAL Q=" << calc_QModularity(sampl_part) << endl; 

	print_part_file_ensemble(results_folder, sampl_part,NVert);		
	
	sprintf(results_aux, "%s_ensemble.txt",time_folder);
	#if COUT_GENERAL == true		
		cout << "** Save Ensemble in " << results_aux << "... **" << endl;
	#endif

	destroy_pop_arrays();

	f=fopen(results_aux,"w"); 	
	fprintf(f,"%lf\n",tTotalIt);		
	fclose(f);
	f=NULL;

	return 0;
}