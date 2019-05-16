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
#ifndef Ensemble_H
#define Ensemble_H

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string.h>
#include <list>

#include "../defines.h"
#include "../spectral_clustering/spectral_clustering_aggr.h"

#define DEF_THRESHOLD 0.5
#define DEF_SEL_PARETO 0.8
#define DEF_MULT_CONSENSUS 2

using namespace std;

class Ensemble{
	private:
		int threshold(double** consensus, int NVert, double p_thr);
		double** consensus_parts(int** pareto_parts, int num_pareto, int NVert);
		int consensus_number(int** pareto_parts, int num_pareto, int NVert);		
		void print_details_ensemble_file(char* details_file_base, int** pareto_parts, int* ens_part, int num_pareto, int NVert);

	public:
		int ensemble_clustering(char* results_folder, Spectral_Clustering* sc, int** pareto_parts, int num_pareto, int NVert, char* details_file_base, double p_thr);
		int ensemble_clustering_2rounds(char* results_folder, Spectral_Clustering* sc, int** pareto_parts, int num_pareto, int NVert, char* details_file_base);
	};

#endif
