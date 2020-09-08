exec: main.o solve.o
	g++ main.cpp include/solve.cpp -o exec
main.o: main.cpp
	g++ -c main.cpp
solve.o: include/solve.cpp
	g++ -c include/solve.cpp