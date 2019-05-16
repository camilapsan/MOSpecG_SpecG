# ARPACK++ v1.2 2/18/2000
# c++ interface to ARPACK code.
# examples/areig/sym directory makefile.

# including other makefiles.

include arpackpp-master/Makefile.inc

# defining areig directory.

AREIG_DIR = $(ARPACKPP_DIR)/examples/areig

# defining cscmat directory.

CSCMAT_DIR = $(ARPACKPP_DIR)/examples/matrices/sym

# compiling and linking all examples.

all: spectral_clustering

spectral_clustering: ensemble.o spectral_clustering3.o util_arrays.o util_eigen.o flow_matrix.o defines.h
	g++ -O3 -o spectral_clustering.out ensemble.o spectral_clustering3.o util_arrays.o util_eigen.o flow_matrix.o -lpthread -lm  -I$(AREIG_DIR) -I$(CSCMAT_DIR) $(SUPERLU_LIB) $(ALL_LIBS) -std=c++11 -fopenmp -Wformat=0

ensemble.o: ensemble/ensemble.cpp ensemble/ensemble.h defines.h
	g++ -O3 -o ensemble.o -c ensemble/ensemble.cpp -std=c++11 -fopenmp -Wformat=0 -lpthread -lm -I$(AREIG_DIR) -I$(CSCMAT_DIR) $(SUPERLU_LIB) $(ALL_LIBS) -std=c++11 -fopenmp -Wformat=0

util_arrays.o: spectral_clustering/util_arrays.cpp spectral_clustering/util_arrays.h defines.h
	g++ -O3 -o util_arrays.o -c spectral_clustering/util_arrays.cpp -std=c++11 -fopenmp -Wformat=0 -lpthread -lm -I$(AREIG_DIR) -I$(CSCMAT_DIR) $(SUPERLU_LIB) $(ALL_LIBS) -std=c++11 -fopenmp -Wformat=0

util_eigen.o: spectral_clustering/util_eigen.cpp spectral_clustering/util_eigen.h defines.h spectral_clustering/areig.h
	g++ -O3 -o util_eigen.o -c spectral_clustering/util_eigen.cpp -std=c++11 -fopenmp -Wformat=0 -lpthread -lm -I$(AREIG_DIR) -I$(CSCMAT_DIR) $(SUPERLU_LIB) $(ALL_LIBS) -std=c++11 -fopenmp -Wformat=0

flow_matrix.o: spectral_clustering/flow_matrix.cpp spectral_clustering/flow_matrix.h defines.h
	g++ -O3 -o flow_matrix.o -c spectral_clustering/flow_matrix.cpp -std=c++11 -fopenmp -Wformat=0 -lpthread -lm -I$(AREIG_DIR) -I$(CSCMAT_DIR) $(SUPERLU_LIB) $(ALL_LIBS) -std=c++11 -fopenmp -Wformat=0

spectral_clustering3.o: spectral_clustering/spectral_clustering_aggr.cpp spectral_clustering/spectral_clustering_aggr.h spectral_clustering/areig.h defines.h 
	g++ -O3 -o spectral_clustering3.o -c spectral_clustering/spectral_clustering_aggr.cpp -std=c++11 -fopenmp -Wformat=0 -lpthread -lm -I$(AREIG_DIR) -I$(CSCMAT_DIR) $(SUPERLU_LIB) $(ALL_LIBS) -std=c++11 -fopenmp -Wformat=0

# defining cleaning rule.

.PHONY:	clean
clean:
	rm -f *~ *.o *g.out spectral_clustering/*.o

# defining pattern rules.
%.o:	%.cc
	$(CPP) -I$(AREIG_DIR) -I$(CSCMAT_DIR) -I$(SUPERLU_DIR) $<

