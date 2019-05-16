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

#ifndef Util_Eigen_H
#define Util_Eigen_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string.h>
#include <list>

#include "../arpackpp-master/include/ardsmat.h"
#include "../arpackpp-master/include/ardnsmat.h"
#include "../arpackpp-master/include/arlnsmat.h"
#include "../arpackpp-master/include/arlsnsym.h"

#include "../defines.h"
#include "util_arrays.h"

using namespace std;

class Util_Eigen{

private:
	
public: 
	static int calc_eigen_values_vectors(double* eig_val_R, double* eig_val_I,	 double* eig_vec,int size, int nnz, double* matrix_csc, int*irow, int*pcol, int PEigen, const char* which);
	static int calc_eigen_values(double* eig_val_R, double* eig_val_I,	 double* eig_vec,int size, int nnz, double* matrix_csc, int*irow, int*pcol, int PEigen, const char* which);
	static void print_eigen_values(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec);	
	static void print_eigen_values_file(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec, int NSampling, char* filename_eigenvalues, char* filename_eigenvectors);	
	static void print_eigen_values_file(int nconv, double* eig_val_R, double* eig_val_I, double* eig_vec, FILE* file);	
};

#endif
