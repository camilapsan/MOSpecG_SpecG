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

#ifndef Util_Arrays_H
#define Util_Arrays_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string.h>
#include <list>
#include <random>

#include "../defines.h"


using namespace std;

class Util_Arrays{

private:


public:	
	/** General alloc and destroy **/
	static double* alloc_1d_double(int size);
	static int* alloc_1d_int(int size);
	static double ** alloc_2d_double(int n1,int n2);
	static int** alloc_2d_int(int n1,int n2);
	static double*** alloc_3d_double(int n1,int n2,int n3);
	
	static double* destroy_1d(double*mat, int n1);
	static int* destroy_1d(int*mat, int n1);
	static double** destroy_2d(double**mat, int n1,int n2);
	static int** destroy_2d(int**mat, int n1,int n2);	
	static double*** destroy_3d(double***mat, int n1,int n2,int n3);
	
	/** General init **/
	static void init_seq_1d(int* mat, int n1);

	/** General print **/
	static void print_1d(double* mat, int n1);
	static void print_1d(int* mat, int n1);
	static void print_2d(double** mat, int n1, int n2);
	
	static void print_1d_file(char* file_name,double* mat, int n1);
	static void print_2d_file(char* file_name,double** mat, int n1, int n2);

	/** General random generator **/
	static int gen_rand_int(int min, int max);
	static double gen_rand_double(double min, double max);
	static void shuffle_array_fisher_yates(int* mat, int size);

	/** Vetorial functions **/
	static void add_array(double* dest, double* source, int size);
	static void sub_array(double* dest, double* source, int size);
	static double norm2_1d(double* array, int size);	
	static double inner_product_1d(double*array1, double* array2, int size);
	static bool compare_equal_array(double* dest, double* source, int size);

	static int convert_matrix_to_csc(double** source, int n1, int n2, double* dest_csc, int* irow, int* pcol);
	static int add_matrix(double **dest_matrix,  double **add_matrix, int n1, int n2, int mult_add_matrix);
	static int add_matrix_lin_to_array(double *dest_array, double **add_matrix,int n1, int n2, int mult_add_matrix);
	static double add_matrix_lin_to_num(double num, double **add_matrix,int n1, int n2, int mult_add_matrix);
};

#endif
