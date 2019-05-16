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

---------------------------------------------------------

Guide for MOSpecG and SpecG codes

---------------------------------------------------------

1 Install Dependencies: 

	ARPACK++
		source: https://github.com/m-reuter/arpackpp
		Unpack arpackpp-master.zip and follow the instructions in README.md.
		Obs: make sure the dependencies LAPACK, BLAS, ARPACK and SuperLU are installed.
			 The script installers might require aditional packages such as curl and cmake.

		You can also try to run:
			sudo apt-get install libarpack++2-dev 

	C++ 2011

---------------------------------------------------------
2. Compile MOSpecG and SpecG:
	make

---------------------------------------------------------
3. To run MOSpecG and SpecG
	./spectral_clustering.out <MO> <IT> <p> <instance> <res_base> <time_base> <pareto_file> NF tau NG NP NO
	
	Obs: < > is used only to indicate the parameters.

	A bash example is provided in file: <run_example.sh>. To run:
		bash run_example.sh

	---------------
	3.1. OUTPUTS

		Partitions: files beginning in <res_base>
		Execution times: files beginning in <time_base> 

	When MO=10, the output partition of MOSpecG-mod is:
		<res_base>"_par_"NG"_"NP"_"NO"_"p"_"IT"_mod.com"
		<time_base>...

	When MO=11, the output of MOSpecG and SpecG are:
		Solutions in the pareto sets:
			<res_base>"_"NF"_"thr"_"NG"_"NP"_"NO"_"p"_"IT"_P0.com": partition gamma_1=0, gamma_2=1-gamma_1
			<res_base>"_"NF"_"thr"_"NG"_"NP"_"NO"_"p"_"IT"_P1.com": partition gamma_1=1*inc, gamma_2=1-gamma_1,  inc=1/(<NF>-1)
			<res_base>"_"NF"_"thr"_"NG"_"NP"_"NO"_"p"_"IT"_P2.com": partition gamma_1=2*inc, gamma_2=1-gamma_1
			...
			<res_base>"_"NF"_"thr"_"NG"_"NP"_"NO"_"p"_"IT"_P"<NF-1>".com": partition gamma_1=1, gamma_2=1-gamma_1
			
			<time_base>...

		Consensus partition: 
			<res_base>"_"NF"_"thr"_"NG"_"NP"_"NO"_"p"_"IT"_ensemble.com"			
			
			<time_base>...

		Values of the objective functions: 
			File <pareto_file>

	---------------
	3.2. INPUTS

	The parameters for the algorithm are: 

		1. MO
			MO=10 to run MOSpecG for maximizing the classical modularity (MOSpecG-mod).
		    MO=11 to run MOSpecG to find Pareto sets and obtain a final consensus partition with SpecG algorithm.

		2. IT: maximum number of iterations for the local search multiplied by 10*****.

		3. p: consider only the largest p*n eigenvalues in absolute value, where n is the number of vertices. 

		4. <instance>: full path to the input graph in pajek format.

		5. <res_base>: base name of the file to output the partitions in between "".
			The algorithm will automatically generate output files.

		6. <time_base>: path and base name of the file to output the execution time in seconds in between "".
			The algorithm will automatically generate files containing the execution times in between "".

		7. <pareto_file>: full path to file to save the objective values

		8. NF: number of solutions in the Pareto set --  When M0=10 this input can be set to any value.

		9. tau: threshold parameter of SpecG - When M0=10 this input can be set to any value.

		10. NG: number of generations of the memetic algorithm.

		11. NP: size of the population of the memetic algorithm.

		12. NO: the NO% fittest individuals from the offspring replace the NO% least fit individuals from a current population in the memetic algorithm.

---------------------------------------------------------
