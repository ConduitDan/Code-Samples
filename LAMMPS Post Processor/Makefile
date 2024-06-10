DEPS = $(wildcard ../*.h)
src = $(wildcard *.cpp) $(wildcard ../*.cpp)
OBJ = $(src:.cpp=.o)
CXXFLAGS = -I. -I./.. -fopenmp -std=c++17
# CXX = nvcc
LDFLAGS=-L/share/software/utilities/BOOST/1.7.2_gcc/lib/ -lboost_filesystem
make: CXXFLAGS += -O3
make: $(OBJ)
	$(CXX) -o post_processor $^ $(LDFLAGS) $(CXXFLAGS)
	make clean
	
debug: CXXFLAGS += -O0 -ggdb
debug: $(OBJ)
	$(CXX) -o post_processor $^ $(LDFLAGS) $(CXXFLAGS)
	
clean:
	rm -rf *.o
