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

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <cmath>
#include <iomanip>
	
#define NDEBUG 

#include <assert.h>
#include <string.h>


//TODO: USE RANDOM LIBRARY
//#include <random>
//#include <time.h>
//#include <thread>

#define SMALLZERO 0.0000000005
#define LARGEZERO 0.000005
#define COUT_DEBUG true
#define COUT_POP false
#define COUT_EIGEN false
#define COUT_GENERAL true
#define POND false
#define FILE_SIZE 400
#define SAVE_DETAILS false

#define SAVE_CPLEX_MODEL true
