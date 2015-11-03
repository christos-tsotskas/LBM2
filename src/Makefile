CUDA_SRC= bounceback.cu main.cu streaming.cu check_density.cu convective_bc.cu  global_defines.cu redistribute.cu LBM_cuda.cu relaxation.cu

CUDA_OBJ= $(CUDA_SRC:.cu=.o)

CXX_SRC=v_field.cpp read_reactor.cpp write_full_tecplot.cpp redistribute2.cpp write_tecplot.cpp write_vorticity.cpp redistribute_initialise.cpp init_density.cpp

CXX_OBJ=$(CXX_SRC:.cpp=.o)


EXECUTABLE = cudaLBM

COMPILER= nvcc
CUDA_FLAGS= --compile -O3 -gencode arch=compute_20 code=sm_20 -x cu 
CXX_FLAGS=-O3 --compile  -x c++ 

LINKER=nvcc

EXECUTION = mpirun
EXECUTION_FLAGS=-n 4 

all: $(CUDA_SRC) $(CXX_SRC) $(EXECUTABLE)

$(EXECUTABLE): $(CUDA_OBJ) $(CXX_OBJ)
	$(LINKER) -link $(CUDA_OBJ) $(CXX_OBJ) -o $@

.cu.o:
	$(COMPILER)  $(CUDA_FLAGS) $< -o $@

.cpp.o:
	$(COMPILER)  $(CXX_FLAGS) $< -o $@
	

run: $(EXECUTABLE)
	./$(EXECUTABLE)

clean: $(EXECUTABLE)
	rm $(EXECUTABLE) $(CUDA_OBJ) $(CXX_OBJ)
