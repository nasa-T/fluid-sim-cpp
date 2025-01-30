fluid: fluid.cpp fluid.h
	g++ -std=c++17 -o fluid fluid.cpp -lSDL2

display: simDisplayer.cpp
	g++ -std=c++17 -o display simDisplayer.cpp -lSDL2