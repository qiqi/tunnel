euler2d:	euler2d.o
	g++ -std=c++11 -o euler2d euler2d.o -L/usr/lib/openmpi/lib -lmpi -lmpi_cxx

euler2d.o:	euler2d.cpp PngWriter.hpp
	g++ -std=c++11 -I/usr/include/openmpi -c euler2d.cpp
