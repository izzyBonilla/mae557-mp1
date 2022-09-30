# C++ Compiler Choice
CXX = g++

# Will look into adding parallelism later
OPTFLAGS = -O3
DEBUGFLAGS = -g

submit: burger.cpp burger.hpp
	$(CXX) $(OPTFLAGS) -o burger.out burger.cpp

debug: burger.cpp burger.hpp
	$(CXX) $(DEBUGFLAGS) -o  burger.out burger.cpp

.PHONY: clean
clean:
	rm *.out