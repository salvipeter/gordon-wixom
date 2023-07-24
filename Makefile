all: gordon-wixom-test

LIBGEOM=../libgeom

# Debug
# CXXFLAGS=-std=c++20 -Wall -pedantic -O0 -g -fsanitize=address -I$(LIBGEOM)
# LIBS=-L$(LIBGEOM)/release -lgeom -lasan

# Release
CXXFLAGS=-std=c++20 -Wall -pedantic -O3 -I$(LIBGEOM)
LIBS=-L$(LIBGEOM)/release -lgeom

gordon-wixom-test: gordon-wixom.o gordon-wixom-test.o
	$(CXX) -o $@ $^ $(LIBS)
