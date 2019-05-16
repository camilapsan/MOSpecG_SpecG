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

#include "flow_matrix.h"

Flow_Matrix::Flow_Matrix(){
	matrix_csc = NULL;
	irow = NULL;
	pcol = NULL;
}

Flow_Matrix::Flow_Matrix(int param_size){
	size = param_size;		
	nnz_max = size*size;
}

Flow_Matrix::~Flow_Matrix(){	
	if(matrix_csc!=NULL){matrix_csc=Util_Arrays::destroy_1d(matrix_csc,nnz);}
	if(irow!=NULL){irow=Util_Arrays::destroy_1d(irow,nnz);}	
	if(pcol!=NULL){pcol=Util_Arrays::destroy_1d(pcol,size+1);}
}

int Flow_Matrix::destroy_data(){	
	matrix_csc=Util_Arrays::destroy_1d(matrix_csc,nnz);
	irow=Util_Arrays::destroy_1d(irow,nnz);
	pcol=Util_Arrays::destroy_1d(pcol,size+1);

	return 0;
}


int Flow_Matrix::initialize_params(int param_size){
	size = param_size;		
	nnz_max = size*size;
	cout << "..D Flow Matrix: nnz_max=" << nnz_max << endl;

	is_newman=FALSE;
}


int Flow_Matrix::initialize_params(int param_size, bool param_newman){
	size = param_size;		
	nnz_max = size*size;
	cout << "..D Flow Matrix: nnz_max=" << nnz_max << endl;

	is_newman=param_newman;
}

//O(m*m)
int Flow_Matrix::construct_calc_KNumber(int nnz_edge_csc, double* edge_csc, int* edge_csc_irow, int* edge_csc_pcol, double* deg, int NVert, char* eig_folder){	
	cout << "D.. nnz_max=" << nnz_max << endl;
	construct_flow_matrix(nnz_edge_csc,edge_csc, edge_csc_irow, edge_csc_pcol, deg);
	cout << "D.. nnz=" << nnz << endl;
	
	if(is_newman==FALSE){
		return calc_K_number(deg, NVert, eig_folder);
	}else{
		return calc_K_number_newman(deg, NVert, eig_folder);
	}
}

int Flow_Matrix::calc_K_number_newman(double* deg, int NVert, char* eig_folder){

	int K=0;
	double radius=0;
	int nconv=0;
	int PEigen=NVert/10 	;///NVert;
	double mean_d=0, mean_d_aux=0;
	double largest_complex=-1000;//TODO
	double largest_real=-1000;//TODO	

	
	double* eig_val_R = Util_Arrays::alloc_1d_double(PEigen+2);   // Real part of the eigenvalues.	
	double* eig_val_I = Util_Arrays::alloc_1d_double(PEigen+2); // not used //Util_Arrays::alloc_1d_double(KNumber+1);   // Imaginary part of the eigenvalues.  	
  	double* eig_vec = Util_Arrays::alloc_1d_double((PEigen+2)*2*size); //times 2 if has imaginary eigenvalues // new double[272];// Util_Arrays::alloc_1d_double(272);   		 

  	nconv = Util_Eigen::calc_eigen_values(eig_val_R, eig_val_I, eig_vec, size, nnz, matrix_csc, irow, pcol, PEigen, "LR");  	
 	
  	// Util_Eigen::print_eigen_values(nconv, eig_val_R, eig_val_I, eig_vec);
  	FILE *f = fopen(eig_folder,"w");

  	cout << "(only val) print file: " << eig_folder << endl;
  	Util_Eigen::print_eigen_values_file(nconv, eig_val_R, eig_val_I, eig_vec, f);
  	cout << "print file ok" << endl;
  	fclose(f);
  	f=NULL;

  	//largest eigenvalue
	
  	//if eigenvalue is complex
  	for (int i=0; i<nconv; i++) {
	    //cout << "\t[" << (i) << "]: " << eig_val_R[i] << ", ";
	    //if (eig_val_I[i]!=0.0 && eig_val_R[i] > largest_real) { //complex eigenvalues
  		if (eig_val_I[i]==0 && eig_val_R[i] > largest_real) { //complex eigenvalues
	    	//complex
	      largest_real=eig_val_R[i];
	      //largest_complex=eig_val_I[i];
	    }	   
	}
	
	radius=sqrt(largest_real);
	cout << "..D largest_real=" << largest_real << ";\t radius=" << radius << endl;

	//radius=sqrt( mean( d/(d-1) )/ mean(d) )
	mean_d=0;
	mean_d_aux=0;
	for(int i=0; i<NVert; i++){
		if(deg[i]>1){			
			//mean_d_aux+=(deg[i]/(deg[i]-1));
			//TEST 
			mean_d_aux+=(largest_real/(deg[i]-1));
		}else{			
			//mean_d_aux+=deg[i]/deg[i];
			mean_d_aux+=largest_real/deg[i];

		}

		mean_d+=deg[i];
	}
	mean_d_aux=mean_d_aux/NVert;
	mean_d=mean_d/NVert;
	cout << "..D mean d=" << mean_d <<endl;
	cout << "..D mean d aux=" << mean_d_aux << "square = " << sqrt(mean_d_aux/mean_d) <<endl;
	radius=sqrt(mean_d_aux/mean_d);

	//TEST 
	
	// cout << "..D c=" << largest_real << "+ i=" << largest_complex << " g" << "  o autovalor complexo que apresenta o maior valor na parte real, c.  "<< endl;
	// radius = sqrt(largest_real);
	
	cout << endl;
 	for (int i=0; i<nconv; i++) {
	    //cout << "\t[" << (i) << "]: " << eig_val_R[i] << ", ";
	    //if (eig_val_I[i]==0.0) {  //real eigenvalues
	    	if(fabs(eig_val_R[i]) >= radius){
	    		K++;
	    	}	      
	    //}	   
	}
	// cout << "radius = " << radius << endl;
  	 cout << "..D K=" << K << endl;

  	eig_val_R = Util_Arrays::destroy_1d(eig_val_R,PEigen+2);
  	eig_val_I = Util_Arrays::destroy_1d(eig_val_I,PEigen+2);
  	eig_vec = Util_Arrays::destroy_1d(eig_vec,(PEigen+2)*2*size);

  	return K;
}

int Flow_Matrix::calc_K_number(double* deg, int NVert, char* eig_folder){
	 cout << "calc k number" << endl;	 
	int K=0;
	double radius=0;
	int nconv=0;
	//int PEigen=size/5;///NVert;
	int PEigen=NVert;///NVert;
	double mean_d=0, mean_d_aux=0;
	double largest_complex=-1000;//TODO
	double largest_real=-1000;//TODO	


	//O número de clusters pode ser estimado pelo número de autovalores reais da matriz que são maiores que c.
	//calc largest eigen
	double* eig_val_R = Util_Arrays::alloc_1d_double(PEigen+2);   // Real part of the eigenvalues.	
	double* eig_val_I = Util_Arrays::alloc_1d_double(PEigen+2); // not used //Util_Arrays::alloc_1d_double(KNumber+1);   // Imaginary part of the eigenvalues.  	
  	double* eig_vec = Util_Arrays::alloc_1d_double((PEigen+2)*2*size); //times 2 if has imaginary eigenvalues // new double[272];// Util_Arrays::alloc_1d_double(272);   		 
  	nconv = Util_Eigen::calc_eigen_values_vectors(eig_val_R, eig_val_I, eig_vec, size, nnz, matrix_csc, irow, pcol, PEigen, "LR");  	
 	
  	// Util_Eigen::print_eigen_values(nconv, eig_val_R, eig_val_I, eig_vec);
  	FILE *f = fopen(eig_folder,"w");

  	cout << "print file: " << eig_folder << endl;
  	Util_Eigen::print_eigen_values_file(nconv, eig_val_R, eig_val_I, eig_vec, f);
  	cout << "print file ok" << endl;
  	fclose(f);
  	f=NULL;

  	//largest eigenvalue
	
  	//if eigenvalue is complex
  	for (int i=0; i<nconv; i++) {
	    //cout << "\t[" << (i) << "]: " << eig_val_R[i] << ", ";
	    //if (eig_val_I[i]!=0.0 && eig_val_R[i] > largest_real) { //complex eigenvalues
  		if (eig_val_I[i]==0 && eig_val_R[i] > largest_real) { //complex eigenvalues
	    	//complex
	      largest_real=eig_val_R[i];
	      //largest_complex=eig_val_I[i];
	    }	   
	}
	
	radius=sqrt(largest_real);
	cout << "..D largest_real=" << largest_real << ";\t radius=" << radius << endl;

	cout << endl;
 	for (int i=0; i<nconv; i++) {
	    //cout << "\t[" << (i) << "]: " << eig_val_R[i] << ", ";
	    //if (eig_val_I[i]==0.0) {  //real eigenvalues
	    	if(fabs(eig_val_R[i]) >= radius){
	    		K++;
	    	}	      
	    //}	   
	}
	// cout << "radius = " << radius << endl;
  	// cout << "..D K=" << K s<< endl;

  	eig_val_R = Util_Arrays::destroy_1d(eig_val_R,PEigen+2);
  	eig_val_I = Util_Arrays::destroy_1d(eig_val_I,PEigen+2);
  	eig_vec = Util_Arrays::destroy_1d(eig_vec,(PEigen+2)*2*size);

  	return K;
}

//O(m*m)
void Flow_Matrix::construct_flow_matrix(int nnz_edge_csc, double* edge_csc, int* edge_csc_irow, int* edge_csc_pcol, double* deg){	
	int b=0, ind1=0,ind2=0,ind_f=0;
	int i=0, j=0, k=0, l=0;  
	int count_nnz=0;
	assert(matrix_csc!=NULL);
	 cout << "..T construct.. nnz=" << nnz_edge_csc << endl;
	 cout << "..T construct.. size=" << size << endl;
	//for each col (edge) ind1


	l=0;
	for(ind2=0;ind2<nnz_edge_csc;ind2++){		
		//ind2: (k,l)

		k = edge_csc_irow[ind2];	
		
		//assert(l+1 <= NVert);			
		if(edge_csc_pcol[l] < edge_csc_pcol[l+1]){
			//only if edge exists.... edge_csc do not evaluate if edge exists!!! TODO: new structure
			//column of edge ind 2 starts at position ind_f			
			j=0;
			//for each line (edge) ind1
			for(ind1=0;ind1<nnz_edge_csc;ind1++){					
				//ind1:(i,j)		
				i = edge_csc_irow[ind1];	
				//assert(j+1 <= NVert);	
				
				//the corresponding row of position ind_f is (edge) ind1
				if(edge_csc_pcol[j] < edge_csc_pcol[j+1]){

					//(ind1,ind2) -> (i,j)(k,l)
					if(i!=l && j==k){
						count_nnz++;
					}
					
					if((ind1+1) >= edge_csc_pcol[j+1]){
						j++;
					}
				}
				else{
					j++;
				}
			}			
			if((ind2+1) >= edge_csc_pcol[l+1]){
				l++;
			}
		}
		else{
			l++;
		}
	}

	assert(nnz==count_nnz);
	nnz=count_nnz;
	
	matrix_csc = Util_Arrays::alloc_1d_double(nnz);
	irow= Util_Arrays::alloc_1d_int(nnz); // Row index of all nonzero elements of A.
  	pcol=Util_Arrays::alloc_1d_int(size+1);  // Pointer to the beginning of each column (in irow and A).

  	///TODO: function to construct is_newman separately

	l=0;
	for(ind2=0;ind2<nnz_edge_csc;ind2++){		
		//ind2: (k,l)
		k = edge_csc_irow[ind2];	
		
		//assert(l+1 <= NVert);	
		pcol[ind2]=ind_f; 
		if(edge_csc_pcol[l] < edge_csc_pcol[l+1]){
			//only if edge exists.... edge_csc do not evaluate if edge exists!!! TODO: new structure

			//column of edge ind 2 starts at position ind_f			
			j=0;
			//for each line (edge) ind1
			for(ind1=0;ind1<nnz_edge_csc;ind1++){					
				//ind1:(i,j)		
				i = edge_csc_irow[ind1];	

				//assert(j+1 <= NVert);	
				//the corresponding row of position ind_f is (edge) ind1
				if(edge_csc_pcol[j] < edge_csc_pcol[j+1]){

					//(ind1,ind2) -> (i,j)(k,l)
					if(i!=l && j==k){
						if(is_newman==FALSE){
							matrix_csc[ind_f]=(1.0);
						}else{
							if(deg[i]>1){		
								// cout << "Ind=" << ind1 << "(i=" << i << ",j=" << j << ")" << " --> ";
								// cout << "Ind=" << ind2 << "(k=" << k << ",l=" << l << ")" << " --> ";			
								// cout << "= 1" << endl;
								matrix_csc[ind_f]=(1.0/(deg[i]-1));
							}
							else{
								matrix_csc[ind_f]=(1.0/(deg[i]));
							}
						}
						irow[ind_f]=ind1;
						ind_f++;
					}
					// else{
					// 	cout << "Ind=" << ind1 << "(i=" << i << ",j=" << j << ")" << " --> ";
					// 	cout << "Ind=" << ind2 << "(k=" << k << ",l=" << l << ")" << " --> ";			
					// 	cout << "= 0" << endl;
					// }
					// Do not store 0 positions!! else{//assert
					// 	matrix_csc[ind_f]=0;
					// }
					
					if((ind1+1) >= edge_csc_pcol[j+1]){
						j++;
					}
				}
				else{
					j++;
				}
			}
			
			if((ind2+1) >= edge_csc_pcol[l+1]){
				l++;
			}
		}
		else{
			l++;
		}
	}

	pcol[size]=ind_f;
	nnz=ind_f;
	// cout << "pcol" << endl;
	// Util_Arrays::print_1d(pcol, size+1);

	// cout << "irow" << endl;
	// Util_Arrays::print_1d(irow, size);

	// cout << "data" << endl;

	// cout << "..size=" << size << endl<<endl;
	// cout << "..T DEBUG TRY ind_f=" << ind_f << endl<<endl;
}