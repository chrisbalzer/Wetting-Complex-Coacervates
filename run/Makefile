# Select compiler and libraries to link
CXX=g++ -fopenmp -std=c++17
DIR_SRC=../src
DIR_INC=-I ../src/external/
CFLAGS=-Wall -g -O3 $(DIR_INC)
DIR_LIB=
LDFLAGS=

# Location of source files
CPP_FILES=$(wildcard ${DIR_SRC}/*.cpp) $(wildcard ${DIR_SRC}/*/*.cpp) $(wildcard ${DIR_SRC}/*/*/*.cpp) $(wildcard ${DIR_SRC}/*/*/*/*.cpp)
CPP_FILES:=$(filter-out ${DIR_SRC}/external/%, $(CPP_FILES))
SRC=$(CPP_FILES)
OBJ=$(SRC:.cpp=.o)

# Name of executable
TARGET=WettingCC

$(TARGET): $(OBJ)
	$(CXX) ${DIR_LIB} -o $@ $(OBJ) ${LDFLAGS}

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

default: $(TARGET)
	-rm $(OBJ)

.PHONY : clean
clean:
	rm $(TARGET) $(OBJ)
