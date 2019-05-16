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

#ifndef Flow_Matrix_H
#define Flow_Matrix_H

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
#include "util_eigen.h"

using namespace std;

class Flow_Matrix{

private:
	//Non backtracking matrix
	int size;
	double* matrix_csc;	
	int *irow;
	int *pcol;
	int nnz_max;	
	int nnz;
	bool is_newman;

public: 
	Flow_Matrix();
	Flow_Matrix(int param_size);
	~Flow_Matrix();

	int destroy_data();

	int initialize_params(int param_size, bool param_newman);
	int initialize_params(int param_size);
	int construct_calc_KNumber(int nnz_edge_csc, double* edge_csc, int* edge_csc_irow, int* edge_csc_pcol, double* deg, int NVert, char* eig_folder);
	int calc_K_number(double* deg, int NVert, char* eig_folder);
	int calc_K_number_newman(double* deg, int NVert, char* eig_folder);
	void construct_flow_matrix(int nnz_edge_csc, double* edge_csc, int* edge_csc_irow, int* edge_csc_pcol, double* deg);

};

#endif
