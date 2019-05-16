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

#include "util_arrays.h"

/************************************ 
	General alloc and destroy 
************************************/	

double*** Util_Arrays::alloc_3d_double(int n1,int n2, int n3){
	double*** mat=NULL;
	if ((mat = (double ***) calloc(sizeof(double **),n1)) == NULL) {
		printf("Erro na alocacao de memoria n1.\n");
		exit(1);
	}
	for(int i=0; i<n1; i++) {
		if ((mat[i] = (double **) calloc(sizeof(double*),n2)) == NULL) {
			printf("Erro na alocacao de memoria n2.\n");
			exit(1);
		}
		for(int j=0; j<n2; j++) {
			if ((mat[i][j] = (double *) calloc(sizeof(double),n3)) == NULL) {
				printf("Erro na alocacao de memoria n2.\n");
				exit(1);
			}	
		}
	}	
	return mat;
}

int** Util_Arrays::alloc_2d_int(int n1,int n2){
	int** mat=NULL;
	if ((mat = (int **) calloc(sizeof(int *),n1)) == NULL) {
		printf("Erro na alocacao de memoria n1.\n");
		exit(1);
	}
	for(int i=0; i<n1; i++) {
		if ((mat[i] = (int *) calloc(sizeof(int),n2)) == NULL) {
			printf("Erro na alocacao de memoria n2.\n");
			exit(1);
		}		
	}	
	return mat;
}

double** Util_Arrays::alloc_2d_double(int n1,int n2){
	double** mat=NULL;
	if ((mat = (double **) calloc(sizeof(double *),n1)) == NULL) {
		printf("Erro na alocacao de memoria n1.\n");
		exit(1);
	}
	for(int i=0; i<n1; i++) {
		if ((mat[i] = (double *) calloc(sizeof(double),n2)) == NULL) {
			printf("Erro na alocacao de memoria n2.\n");
			exit(1);
		}		
	}	
	return mat;
}

double* Util_Arrays::alloc_1d_double(int size){
	double* ret1d=NULL;
	if ((ret1d = (double *) calloc(sizeof(double ),size)) == NULL) {
		printf("Erro na alocacao de memoria.\n");
		exit(1);
	}
	return ret1d;
}

int* Util_Arrays::alloc_1d_int(int size){
	int* ret1d=NULL;
	if ((ret1d = (int *) calloc(sizeof(int ),size)) == NULL) {
		printf("Erro na alocacao de memoria.\n");
		exit(1);
	}
	return ret1d;
}

double*** Util_Arrays::destroy_3d(double***mat, int n1,int n2, int n3){	
	for(int i=0; i<n1; i++) {
		for(int j=0; j<n2; j++) {
			free(mat[i][j]);			
			mat[i][j]=NULL;
		}
		free(mat[i]);
		mat[i]=NULL;
	}	
	free(mat);
	mat=NULL;
	return mat;
}

double** Util_Arrays::destroy_2d(double**mat, int n1,int n2){	
	for(int i=0; i<n1; i++) {
		free(mat[i]);			
		mat[i]=NULL;
	}	
	free(mat);
	mat=NULL;
	return mat;
}

int** Util_Arrays::destroy_2d(int**mat, int n1,int n2){	
	for(int i=0; i<n1; i++) {
		free(mat[i]);			
		mat[i]=NULL;
	}	
	free(mat);
	mat=NULL;
	return mat;
}

double* Util_Arrays::destroy_1d(double*mat, int n1){	

	free(mat);
	mat=NULL;
	return mat;
}

int* Util_Arrays::destroy_1d(int*mat, int n1){	

	free(mat);
	mat=NULL;
	return mat;
}

/************************************
 General init 
 ************************************/
void Util_Arrays::init_seq_1d(int* mat, int n1){
	//initialize.
	for(int i=0;i<n1;i++){
		mat[i]=i;
	}	
}

/************************************
 General print 
 ************************************/

void Util_Arrays::print_1d(double* mat, int n1){
	for(int i=0;i<n1;i++){
		cout << "[" << i << "]" << mat[i] << ", ";
	}
	cout << endl;
}


void Util_Arrays::print_1d_file(char* file_name,double* mat, int n1){
	FILE* file = fopen(file_name,"w");

	for(int i=0;i<n1;i++){
		fprintf(file,"%lf", mat[i]);
	}
	fclose(file);
	file=NULL;
}

void Util_Arrays::print_1d(int* mat, int n1){
	for(int i=0;i<n1;i++){
		cout << "[" << i << "]" << mat[i] << ", ";		
	}
	cout << endl;
}


void Util_Arrays::print_2d(double** mat, int n1, int n2){
	for(int i=0;i<n1;i++){
		cout << "i=" << i << endl;
		for(int j=0;j<n2;j++){
			//cout << "[" << i << "][" << j << "]="  << mat[i][j] << ", ";
			cout << "[" << j << "]="  << mat[i][j] << ", ";
		}
		cout << endl << endl;
	}	
}

void Util_Arrays::print_2d_file(char* file_name,double** mat, int n1, int n2){
	FILE* file = fopen(file_name,"w");

	for(int i=0;i<n1;i++){
		fprintf(file,"i=%d; ",i);

		for(int j=0;j<n2;j++){
			fprintf(file,"%lf; ",mat[i][j]);			
		}
		fprintf(file,"\n");					
	}
	fclose(file);
	file=NULL;
}
/************************************
 General random generator 
 ************************************/

double Util_Arrays::gen_rand_double(double min, double max) {
    random_device rd;
     uniform_real_distribution<double> distribution(min, max);
     return distribution(rd);
}

int Util_Arrays::gen_rand_int(int min, int max) {	
     random_device rd;
     uniform_int_distribution<int> distribution(min, max);
     return distribution(rd);
}

void Util_Arrays::shuffle_array_fisher_yates(int* mat, int size) {    
	int r;
	int aux;
	for(int i=0; i<size; i++) {  
		mat[i]=i;
	}
 	for(int i=size-1; i>0; i--) { 
        r = Util_Arrays::gen_rand_int(0, (i)); 
        aux=mat[i];
        mat[i]=mat[r];
        mat[r]=aux;
    }    
}

/************************************
 Vetorial functions
 ************************************/
/*
	Check difference between source and dest
*/	
bool Util_Arrays::compare_equal_array(double* dest, double* source, int size){
	double sum_diff=0;
	for(int l=0;l<size;l++){
		sum_diff = sum_diff + (dest[l] - source[l]);
	}
	if(abs(sum_diff) < SMALLZERO){
		return true;
	}
	return false;
}

/*
	Add ri contribution to RGroup of cluster s
*/	
void Util_Arrays::add_array(double* dest, double* source, int size){
	for(int l=0;l<size;l++){
		dest[l] += source[l];
	}
}

/*
	Sub ri contribution to RGroup of cluster s
*/	
void Util_Arrays::sub_array(double* dest, double* source, int size){
	for(int l=0;l<size;l++){
		dest[l] -= source[l];
	}
}

double Util_Arrays::norm2_1d(double* array, int size){
	double norm=0;
	for(int l=0;l<size;l++){
		norm+=array[l]*array[l];
	}
	return norm;
}

double Util_Arrays::inner_product_1d(double*array1, double* array2, int size){
	double prod=0;
	for(int i=0;i<size;i++){
		prod+=array1[i]*array2[i];
	}
	return prod;
}

int Util_Arrays::convert_matrix_to_csc(double** source, int n1, int n2, double* dest_csc, int* irow, int* pcol){

	int	ind=0;  
	for(int j=0;j<n2;j++){
	  	pcol[j]=ind;
	  	
		for(int i=0;i<n1;i++){
			if(source[i][j]>0){
				dest_csc[ind]= source[i][j];
				irow[ind]=i;
	  			ind++;				  			
	  		}	  			
		}	
	}
	
	pcol[n1]=ind;

	return ind;
}

int Util_Arrays::add_matrix(double **dest_matrix, double **add_matrix,int n1, int n2, int mult_add_matrix){
	for(int i=0;i<n1;i++){
		for(int j=0;j<n2;j++){
			dest_matrix[i][j] += (double(mult_add_matrix)*add_matrix[i][j]);
		}		
	}	
	return 0;
}

int Util_Arrays::add_matrix_lin_to_array(double *dest_array, double **add_matrix,int n1, int n2, int mult_add_matrix){
	double sum=0;
	for(int i=0;i<n1;i++){
		sum=0;
		for(int j=0;j<n2;j++){
			sum += (double(mult_add_matrix)*add_matrix[i][j]);
		}	
		dest_array[i]+= sum;
	}	
	return 0;
}

double Util_Arrays::add_matrix_lin_to_num(double num, double **add_matrix,int n1, int n2, int mult_add_matrix){
	double sum=0;

	for(int i=0;i<n1;i++){	
		for(int j=0;j<n2;j++){
			sum += (double(mult_add_matrix)*add_matrix[i][j]);
		}	
	}
	num+= sum;
		
	return num;
}
