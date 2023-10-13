# by default select gcc
CXX=g++
CXXFLAGS= -O3 -g -DNDEBUG -march=native -mtune=native -funroll-loops
CXXFLAGS_DEBUG= -O0 -g 

#make sure the path to the boost library is correct
CXXINCLUDES= -I /home/xniu/boost_1_81_0/
# for error checking, use: -fsanitize=address

#specify input variables

#file containing the COO storage of the tensor.  Cannot be left blank
TENSOR= tensors/nips.tns

#If the indices start with 1, use OFFSET=1
#If the indices start with 0, use OFFSET=0
OFFSET= 1

#separator between each column entry of tensor file.  Cannot be left blank
SEP= ' '

#tensor dimensions, e.g. for a 100x100x100 tensor: TENSOR_DIM = 100,100,100.  Cannot be left blank
TENSOR_DIM= 2482,2862,14036,17

#block size for blocked Bit-If. Cannot be left blank
#use computation about what blocksize is best given dimension
BLOCK_SIZE= 1000

#specify on which mode the TVM should be computed
#e.g. for TVM along the first mode: MODE = 1
#e.g. for TVM along the second mode: MODE = 2
MODE= 2

#specity traversal curve. Default is hierachical traversal.
#available traversals: hierachical (TRAVERSAL = hierarchical), Z-curve (TRAVERSAL = zcurve), Hilbert curve (TRAVERSAL = hilbert)
TRAVERSAL= zcurve

#file containing the values of the vector for TVM. Must be a single column file. 
#please make sure to provide the correct vector file for the TVM, 
#i.e. the number of entries in the vector file should match the tensor dimension
#if left blank, no TVM will be performed and only the Bit-IF datastructure is constructed
VECTOR= vectors/vec_2862.tns #check whether vec size == dimension

.SUFFIXES: .cpp

.PHONY: all run clean

all: bitif bitif_block coo

bitif:
	${CXX} ${CXXINCLUDES} main_bitif.cpp src/*.cpp ${CXXFLAGS} -o bitif

bitif_block: 
	${CXX} ${CXXINCLUDES} main_bitif_blocked.cpp src/*.cpp ${CXXFLAGS} -o bitif_block
 
coo:
	${CXX} ${CXXINCLUDES} main_coo.cpp src/*.cpp ${CXXFLAGS} -o coo

run_blocked: bitif_block
	@./bitif_block $(TENSOR) $(OFFSET) $(SEP) $(TENSOR_DIM) $(BLOCK_SIZE) $(MODE) $(TRAVERSAL) $(VECTOR)

run: bitif
	@./bitif $(TENSOR) $(OFFSET) $(SEP) $(TENSOR_DIM) $(MODE) $(TRAVERSAL) $(VECTOR)

run_coo: coo
	@./coo $(TENSOR) $(OFFSET) $(SEP) $(TENSOR_DIM) $(MODE) $(TRAVERSAL) $(VECTOR)

clean:
	rm -f bitif bitif_block coo
