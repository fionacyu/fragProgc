main: mgraph.o input.o load_data.o gcharacter.o main.o
	g++-11 -std=c++17 -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -o main main.o mgraph.o input.o load_data.o gcharacter.o

mgraph.o: mgraph.cpp
	g++-11 -std=c++17 -c mgraph.cpp

input.o: input.cpp 
	g++-11 -std=c++17 -c input.cpp

load_data.o: load_data.cpp
	g++-11 -std=c++17 -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -c load_data.cpp

gcharacter.o: gcharacter.cpp
	g++-11 -std=c++17 -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -c gcharacter.cpp

main.o: main.cpp
	g++-11 -std=c++17 -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 -c main.cpp

clean:
	rm *.o