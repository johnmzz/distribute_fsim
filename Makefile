# makefile syntax for a rule:
# target: prerequisites
# <TAB> recipe

CODE = ./code/
BOOST_ROOT = /import/vldb01/2/scratch/mazhuo/boost/boost_1_77_0

# flags:
# -c 			- compile source code to object file
# -Ofast 		- highest level of optimization, longer compilation time
# -fopenmp		- 多核编程框架, link OpenMP runtime library
# -I/path		- path to header files
# -o 			- output executable file
# std=c++14		- use C++14 standard

all: calsim
calsim: main_MPI.o Graph.o Initializer.o InitializerMPI.o MatchMethod.o Denominator.o FracSimulation.o DependencyGraph.o
	mpicxx -Ofast  -fopenmp  main_MPI.o Graph.o Initializer.o InitializerMPI.o MatchMethod.o Denominator.o FracSimulation.o DependencyGraph.o -o calsim -std=c++14
	rm *.o
main_MPI.o:
	mpicxx -c  -Ofast  -fopenmp  $(CODE)main_MPI.cpp -I$(BOOST_ROOT) -std=c++14
FracSimulation.o:
	mpicxx -c  -Ofast  -fopenmp  $(CODE)FracSimulation.cpp -I$(BOOST_ROOT) -std=c++14
Denominator.o:
	g++ -c  -Ofast  -fopenmp  $(CODE)Denominator.cpp -I$(BOOST_ROOT) -std=c++14
MatchMethod.o:
	g++ -c  -Ofast  -fopenmp  $(CODE)MatchMethod.cpp -I$(BOOST_ROOT) -std=c++14
Initializer.o:
	g++ -c  -Ofast  -fopenmp  $(CODE)Initializer.cpp -I$(BOOST_ROOT) -std=c++14
InitializerMPI.o:
	g++ -c  -Ofast  -fopenmp  $(CODE)InitializerMPI.cpp -I$(BOOST_ROOT) -std=c++14
Graph.o:
	g++ -c  -Ofast  -fopenmp  $(CODE)Graph.cpp -I$(BOOST_ROOT) -std=c++14
DependencyGraph.o:
	mpicxx -c  -Ofast  -fopenmp  $(CODE)DependencyGraph.cpp -I$(BOOST_ROOT) -std=c++14
