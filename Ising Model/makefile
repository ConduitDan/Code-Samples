TEST = $(wildcard test/*.cpp)
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)

CXXFLAGS = -I. -I./src
LDFLAGS = 

make: $(OBJ)
	$(CXX) -o Ising $^ $(CXXFLAGS) $(LDFLAGS)
	make clean


debug: CXXFLAGS += -O0 -ggdb
debug: make
debug: gdb Ising

check: LDFLAGS += -lboost_unit_test_framework
check: $(TEST)
		$(CXX) -o test/Test $^ $(CXXFLAGS) $(LDFLAGS)
		./test/Test -p
		
clean:
	rm -rf src/*.o