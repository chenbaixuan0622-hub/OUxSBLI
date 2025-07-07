# Makefile
CUDA_PATH       ?= /usr/local/cuda
PYTHON_VERSION  := 3.13
PYTHON_INC      := $(shell python$(PYTHON_VERSION)-config --includes)
CUDA_INC        := -I$(CUDA_PATH)/include
PYBIND11_INC    := $(shell python3 -m pybind11 --includes)

NVCC            := $(CUDA_PATH)/bin/nvcc
CXX             := g++
CXXFLAGS        := -O3 -fPIC -std=c++14

all: cufd.so

flux.o: flux.cu
	$(NVCC) -Xcompiler -fPIC -O3 -c $< -o $@

rk.o: rk.cu
	$(NVCC) -Xcompiler -fPIC -O3 -c $< -o $@

wrapper.o: wrapper.cpp
	$(CXX) $(CXXFLAGS) $(PYTHON_INC) $(CUDA_INC) $(PYBIND11_INC) -c $< -o $@

cufd.so: rk.o flux.o wrapper.o
	$(CXX) -shared -o $@ $^ -L$(CUDA_PATH)/lib64 -lcudart

clean:
	rm -f *.o *.so

