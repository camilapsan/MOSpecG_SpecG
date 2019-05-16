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

#include "ensemble.h"

int Ensemble::ensemble_clustering(char* results_folder, Spectral_Clustering* sc, int** pareto_parts, int num_pareto, int NVert, char* details_file_base, double p_thr=DEF_THRESHOLD){
	int nnz=0;
	
	int k_consensus=1;
	//read partitions
	//sum the number of times vertices are in the same group
	double** consensus = consensus_parts(pareto_parts, num_pareto, NVert);
	k_consensus = consensus_number(pareto_parts, num_pareto, NVert);
	cout << endl << "Consensus matrix" << endl;
	// Util_Arrays::print_2d(consensus,NVert,NVert);

	cout << endl << "Consensus number=" << k_consensus << endl;

	//threshold + recover conexity
	cout << "Threshold" << endl;
	nnz=threshold(consensus, NVert, p_thr);

	//cout << "****************" << endl;
	//Util_Arrays::print_2d(consensus,NVert,NVert);
	//cout << "****************" << endl;

	int* ens_part=(*sc).start_spectral_clustering_ensemble_mod(k_consensus, consensus, DEF_MULT_CONSENSUS);

	#if SAVE_DETAILS == TRUE 
		print_details_ensemble_file(details_file_base,pareto_parts,ens_part,num_pareto,NVert);
	#endif
	consensus=Util_Arrays::destroy_2d(consensus,NVert,NVert);
	// part_consensus=Util_Arrays::destroy_1d(part_consensus,NVert);
}


int Ensemble::ensemble_clustering_2rounds(char* results_folder, Spectral_Clustering* sc, int** pareto_parts, int num_pareto, int NVert, char* details_file_base){
	int nnz=0;
	
	int k_consensus=1;
	//read partitions
	//sum the number of times vertices are in the same group
	double** consensus = consensus_parts(pareto_parts, num_pareto, NVert);
	k_consensus = consensus_number(pareto_parts, num_pareto, NVert);
	cout << endl << "Consensus matrix" << endl;
	// Util_Arrays::print_2d(consensus,NVert,NVert);

	cout << endl << "Consensus number=" << k_consensus << endl;

	//threshold + recover conexity
	cout << "Threshold" << endl;
	nnz=threshold(consensus, NVert, DEF_THRESHOLD);

	int** pareto_parts_round2=(*sc).start_spectral_clustering_ensemble_MO(k_consensus, consensus, DEF_MULT_CONSENSUS);

	//Ensemble clustering - round 2!
	cout << "Ensemble clustering round 2" << endl;

	ensemble_clustering(results_folder, sc, pareto_parts_round2,num_pareto, NVert, details_file_base);
	consensus=Util_Arrays::destroy_2d(consensus,NVert,NVert);
}

void Ensemble::print_details_ensemble_file(char* details_file_base, int** pareto_parts, int* ens_part, int num_pareto, int NVert){		
	FILE* file;
	char file_name[FILE_SIZE];	
	sprintf(file_name,"%s_ensemble.txt",details_file_base);
	file=fopen(file_name,"w");	

	for(int it=0;it<num_pareto;it++){
		fprintf(file,"part; %d;\t",it);
		for(int i=0;i<NVert;i++){		
			fprintf(file,"%d ",pareto_parts[it][i]);
		}
		fprintf(file,"\n");			
	}

	fprintf(file,"ensemble; ;\t");
	for(int i=0;i<NVert;i++){		
		fprintf(file,"%d ",ens_part[i]);
	}
	
	fprintf(file,"\n");
	fclose(file);
	file=NULL;
}


int Ensemble::threshold(double** consensus, int NVert, double p_thr){
	int i,j,j_larg, nnz=0;
	double sum=0,val_larg=-1;
	double avg_threshold=0;
	
	for(i=0; i<NVert; i++){
		sum=0;

		for(j=0; j<NVert; j++){			
			if(consensus[i][j] > val_larg){
				val_larg=consensus[i][j];
				j_larg=j;
			}

			if(consensus[i][j]<p_thr){
				consensus[i][j] =0;
			}else{
				nnz+=1;
				sum+=consensus[i][j];
			}

			avg_threshold+=consensus[i][j];
		}

		//check if i is disconnected and link to the vert whose consensus value if the highest
		if(sum==0){
			nnz+=1;
			consensus[i][j_larg] = val_larg;
		}
	}
	avg_threshold=avg_threshold/(NVert*NVert);
	cout << "AVG ThRESHOLD = " << avg_threshold << endl;
	
	return nnz; 
}

double** Ensemble::consensus_parts(int** pareto_parts, int num_pareto, int NVert){
	int ind,i,j;
	int k=0;
	int ind_ini, ind_fim, sel_pareto;
	double** consensus = Util_Arrays::alloc_2d_double(NVert,NVert); 
	
	cout << "num_pareto=" << num_pareto << endl;

	cout << 1-DEF_SEL_PARETO << "," << double(num_pareto/2)<< endl; 
	ind_ini = double( double(1-DEF_SEL_PARETO)*double(num_pareto/2.0) ); 
	ind_fim = num_pareto-1-ind_ini;

	cout << "ind_ini=" << ind_ini << ", ind_fim=" << ind_fim << endl;
	sel_pareto=ind_fim-ind_ini+1;

	for(i=0; i<NVert; i++){
		for(j=0; j<NVert; j++){
			consensus[i][j]=0;

			for(ind=ind_ini;ind<=ind_fim;ind++){

				if(pareto_parts[ind][i] == pareto_parts[ind][j]){
					consensus[i][j] += 1;
				}
			}
			consensus[i][j] = consensus[i][j]/sel_pareto;
		}
	}

	return consensus; 
}

int Ensemble::consensus_number(int** pareto_parts, int num_pareto, int NVert){
	int ind,i,j;
	int k_num=0, k_sum=0;
	
	for(ind=0;ind<num_pareto;ind++){	
		k_num=0;
		for(i=0; i<NVert; i++){				
			if(pareto_parts[ind][i] > k_num){
				k_num+=1;
			}			
		}
		k_num+=1;
		k_sum+=k_num;
	}

	k_sum=k_sum/num_pareto;

	return k_sum; 
}
int main(int argc, char *argv[]){
	int MO=0;
	int max_it=10;	
	double p_leading=1;
	char* graph_file;
	char* results_folder;	
	char* time_folder;	
	char* pareto_file;	
	char* details_file_base;
	struct timespec startIt, endIt;
   	clock_gettime(CLOCK_MONOTONIC, &startIt);
   	double tTotalIt=0;
   	int num_pareto=1;
   	int NVert=0;
		
	int n_gen=0;
	int n_pop=0;
	double p_offs=0;
	double p_thr=0;
	
	char results_folder_param[FILE_SIZE];
	char time_folder_param[FILE_SIZE];
	char pareto_file_param[FILE_SIZE];

	if(argc>=5){
		MO = atoi(argv[1]);
		max_it = atoi(argv[2])*10; //multiply by 10 because local search algorithm /10.		
		p_leading = atof(argv[3]);
		graph_file = argv[4];
		results_folder = argv[5];
		time_folder = argv[6];
		pareto_file = argv[7];

		#if COUT_GENERAL == true
		cout << "MO: " << MO << endl;
		cout << "graph_file: " << graph_file << endl;
		cout << "results:" << results_folder << endl;		
		#endif		

		int** pareto_parts = NULL;

		if(MO==10){
			cout << "Parameters - SINGLE" << endl;			
			
			n_gen = atoi(argv[10]);
			n_pop = atoi(argv[11]);
			p_offs = atof(argv[12]);

			printf("PARAMS ADJUST:\n,n_gen=%d\n,n_pop=%d\n,p_offs=%.1lf\n,p_leading=%.1lf\n,max_it/10=%d\n\n",n_gen,n_pop,p_offs,p_leading,max_it/10);							
			sprintf(results_folder_param, "%s_par_%d_%d_%.1lf_%.1lf_%d", results_folder,n_gen,n_pop,p_offs,p_leading,max_it/10);
			sprintf(time_folder_param, "%s_par_%d_%d_%.1lf_%.1lf_%d", time_folder, n_gen,n_pop,p_offs,p_leading,max_it/10);
			sprintf(pareto_file_param, "%s_par_%d_%d_%.1lf_%.1lf_%d", pareto_file, n_gen,n_pop,p_offs,p_leading,max_it/10);
			printf("RESULTS FOLDER: %s_par_%d_%d_%.1lf_%.1lf_%d\n\n", results_folder,n_gen,n_pop,p_offs,p_leading,max_it/10);	

			Spectral_Clustering sc(max_it, p_leading, -1, -1, graph_file, results_folder_param, time_folder_param, pareto_file_param, num_pareto, n_gen, n_pop, p_offs);								
			sc.start_spectral_clustering_single(-1);
		}
		else if(MO==11){
			cout << "Parameters test - MO" << endl;			
			num_pareto = atoi(argv[8]);
			p_thr = atof(argv[9]);
			n_gen = atoi(argv[10]);
			n_pop = atoi(argv[11]);
			p_offs = atof(argv[12]);

			printf("PARAMS ADJUST:\nnum_pareto=%d\n,p_thr=%.1lf\n,n_gen=%d\n,n_pop=%d\n,p_offs=%.1lf\n,p_leading=%.1lf\n,max_it/10=%d\n\n",num_pareto,p_thr,n_gen,n_pop,p_offs,p_leading,max_it/10);				
			printf("RESULTS FOLDER: %s_par_%d_%.1lf_%d_%d_%.1lf_%.1lf_%d\n\n", results_folder,num_pareto,p_thr,n_gen,n_pop,p_offs,p_leading,max_it/10);	
			sprintf(results_folder_param, "%s_%d_%.1lf_%d_%d_%.1lf_%.1lf_%d", results_folder,	num_pareto,p_thr, n_gen,n_pop,p_offs,p_leading,max_it/10);
			sprintf(time_folder_param, "%s_%d_%.1lf_%d_%d_%.1lf_%.1lf_%d", time_folder,num_pareto,p_thr, n_gen,n_pop,p_offs,p_leading,max_it/10);
			sprintf(pareto_file_param, "%s_%d_%.1lf_%d_%d_%.1lf_%.1lf_%d", pareto_file,num_pareto,p_thr, n_gen,n_pop,p_offs,p_leading,max_it/10);


			Spectral_Clustering sc(max_it, p_leading, -1, -1, graph_file, results_folder_param, time_folder_param, pareto_file_param, num_pareto, n_gen, n_pop, p_offs);	
			pareto_parts=sc.start_spectral_clustering_MO(-1);
			NVert = sc.get_NVert();

			//Ensemble clustering
			cout << "Ensemble clustering" << endl;
			Ensemble ensemble;
			ensemble.ensemble_clustering(results_folder_param, &sc, pareto_parts,num_pareto, NVert, details_file_base, p_thr);

		}
	
		clock_gettime(CLOCK_MONOTONIC, &endIt);
  		tTotalIt = (endIt.tv_sec - startIt.tv_sec);
  		tTotalIt += (endIt.tv_nsec - startIt.tv_nsec) / 1000000000.0;

  		cout << "* Total Time = " << tTotalIt << " *" << endl;

  		FILE *f;	
  		char results_aux[FILE_SIZE];
		sprintf(results_aux, "%s.txt",time_folder);	
		f=fopen(results_aux,"w"); 	
		fprintf(f,"%lf\n",tTotalIt);		
		fclose(f);
		f=NULL;

	}
	else{

	}
}
