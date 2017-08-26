CPP   = g++
FLAGS = -std=c++14 -stdlib=libc++ -O3 -Wall -Wextra -ffast-math -fvectorize -funroll-loops -fstrict-vtable-pointers -fstrict-enums -Rpass=loop-vectorize -I ~/cpp_libs

src  := $(shell find .. -maxdepth 2 -iname '*.cpp')
obj  := $(src:.cpp=.o)

all: diplom

diplom: $(obj)
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm -rf *.o diplom
	rm -rf ../*.o
