/**************************************************************************
*  
*
*  Copyright (c) 2018 Camila P.S. Tautenhain, Mariá C.V. Nascimento
*
*  This file is part of MOSpecG/SpecG software.
*
*    MOSpecG/SpecG software is free software: you can redistribute it and/or modify
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

#include "util_eigen.h"

void Util_Eigen::print_eigen_values(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec){
	int indexeig_vec=0;
  	cout << "..P Eigenvalues:" << endl;
  	for (int i=0; i<nconv; i++) {
	    cout << "\t[" << (i) << "]: " << eig_val_R[i] << ", ";
	    if (eig_val_I[i]>=0.0) {
	      cout << " + " << eig_val_I[i] << " I" << endl;
	    }
	    else {
	      cout << " - " << fabs(eig_val_I[i]) << " I" << endl;
	    }
		
  	}
  	cout << endl;
}


void Util_Eigen::print_eigen_values_file(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec, FILE* file){
	int indexeig_vec=0;
  	cout << "..P Eigenvalues:" << endl;
  	cout << file << endl;
  	for (int i=0; i<nconv; i++) {	    	      
	      fprintf(file,"%.25lf;%.25lf\n",eig_val_R[i],eig_val_I[i]);

  	}
  	//cout << endl;
}

void Util_Eigen::print_eigen_values_file(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec, int NSampling, char* filename_eigenvalues, char* filename_eigenvectors){
	int indexeig_vec=0;
  	cout << "..P Eigenvalues:" << endl;

  	FILE* file_eigenvalues;
	file_eigenvalues=fopen(filename_eigenvalues,"w");

	cout << file_eigenvalues << endl;
	
  	for (int i=0; i<nconv; i++) {
	      fprintf(file_eigenvalues,"%.25lf;%.25lf\n",eig_val_R[i],eig_val_I[i]);
  	}
  	cout << "try..." << endl;
  	fclose(file_eigenvalues);
  	file_eigenvalues=NULL;
  	
  	
  	cout << "ok..." << endl;
}

int Util_Eigen::calc_eigen_values_vectors(double* eig_val_R, double* eig_val_I, double* eig_vec, int size, int nnz, double* matrix_csc, int* irow, int* pcol, int PEigen, const char* which){
	int nconv=0;
  
  ARluNonSymMatrix<double, double> matrix(size, nnz, matrix_csc, irow, pcol);
  ARluNonSymStdEig<double> dprob(PEigen, matrix, which, 0, 0.0, 0, NULL, true);
 nconv=dprob.EigenValVectors(eig_vec, eig_val_R, eig_val_I);

 	assert(nconv==PEigen);
  	return nconv;
}


int Util_Eigen::calc_eigen_values(double* eig_val_R, double* eig_val_I, double* eig_vec, int size, int nnz, double* matrix_csc, int* irow, int* pcol, int PEigen, const char* which){
	int nconv=0;
  
  ARluNonSymMatrix<double, double> matrix(size, nnz, matrix_csc, irow, pcol);

  // Defining the eigenvalue problem.

  ARluNonSymStdEig<double> dprob(PEigen, matrix, which, 0, 0.0, 0, NULL, true);

  // Finding eigenvalues.

 	nconv=dprob.Eigenvalues(eig_val_R, eig_val_I);
 	assert(nconv==PEigen);

  	return nconv;
}