# /n/fink1/sportillo/software needs to be replaced by the install location of DNest3 and RJObject
default:
	g++ -Wall -Wextra -pedantic -m64 -O3 -flto -march=native -funroll-loops -I /n/fink1/sportillo/software/include/dnest3 -I/n/fink1/sportillo/software/include/rjobject -c *.cpp
	g++ -Wall -Wextra -pedantic -m64 -O3 -flto -march=native -funroll-loops *.o  -o main -rdynamic -L/n/fink1/sportillo/software/lib -lrjobject -ldnest3 -lgsl -lgslcblas -lboost_thread-mt -lboost_system-mt
	rm -f *.o

profiling:
	g++ -Wall -Wextra -pedantic -m64 -O3 -flto -march=native -funroll-loops -I /n/fink1/sportillo/software/include/dnest3 -I/n/fink1/sportillo/software/include/rjobject -c *.cpp -pg -g
	g++ -Wall -Wextra -pedantic -m64 -O3 -flto -march=native -funroll-loops *.o  -o mainprof -rdynamic -L/n/fink1/sportillo/software/lib -lrjobject -ldnest3 -lgsl -lgslcblas -lboost_thread-mt -lboost_system-mt -pg -g
	rm -f *.o
