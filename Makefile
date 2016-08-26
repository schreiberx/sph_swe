

CFLAGS=-I./local_software/local/include -DSWEET_THREADING=1 -DNUMA_BLOCK_ALLOCATOR_TYPE=2 -fopenmp -pedantic -std=c++14  -Isrc/include

LDFLAGS=-fopenmp -L./local_software/local/lib -lshtns_omp -lrt -lm -lnuma -llapack -lrefblas

FFTW_LINK=-lfftw3_omp -lfftw3


all:	release

debug:
	mkdir -p build
	$(CXX) -g -O0 $(CFLAGS) -c src/main.cpp -o build/main.o
	$(CXX) -o build/sh_example build/main.o $(LDFLAGS) $(FFTW_LINK) -lgfortran

release:
	mkdir -p build
	$(CXX) -DNDEBUG=1 -g -O2 $(CFLAGS) -c src/main.cpp -o build/main.o
	$(CXX) -DNDEBUG=1 -o build/sh_example build/main.o $(LDFLAGS) $(FFTW_LINK) -lgfortran

intel:
	# Link also with MKL instead of FFTW3
	mkdir -p build
	icpc -DNDEBUG=1 -g -O2 $(CFLAGS) -c src/main.cpp -o build/main.o
	icpc -DNDEBUG=1 -o build/sh_example build/main.o  -mkl $(LDFLAGS)

knl:
	# Link also with MKL instead of FFTW3
	mkdir -p build
	icpc -DNDEBUG=1 -mic -g -O2 $(CFLAGS) -c src/main.cpp -o build/main.o -DNUMA_BLOCK_ALLOCATOR_TYPE=4
	icpc -DNDEBUG=1 -mic -o build/sh_example build/main.o $(LDFLAGS)

clean:
	rm -f build/main.o
	rm -f build/sh_example
	rm -f O_*
