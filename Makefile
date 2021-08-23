CC=g++ -std=c++11 -pthread
OPT=-O3
CPPFLAGS=-I$(BOOSTINCLUDEDIR)
LDFLAGS=-L$(BOOSTLIBDIR) -lboost_filesystem -lboost_system -lboost_program_options
DBG=-g

.PHONY: default all clean debug

default: all

all: NB.run

test: Class.cpp Class.hpp NB.cpp NB.hpp test.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp Read.cpp Read.hpp
	$(CC) $(OPT)  Read.cpp Class.hpp Genome.cpp NB.cpp Diskutil.cpp test.cpp -o test.run -lboost_filesystem -lboost_system -lboost_program_options

proteus: Class.cpp Class.hpp NB.cpp NB.hpp main.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp
	$(CC) $(OPT) $(CPPFLAGS) Genome.cpp Class.hpp NB.cpp Diskutil.cpp main.cpp -o NB.run $(LDFLAGS)
	chmod +x NB.run

debug: Class.cpp Class.hpp NB.cpp NB.hpp main.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp Read.cpp Read.hpp
	$(CC) $(OPT) $(DBG) -ggdb Read.cpp Genome.cpp Class.hpp NB.cpp Diskutil.cpp main.cpp -o NB.run -lboost_filesystem -lboost_system -lboost_program_options
	chmod +x NB.run

NB.run: Class.cpp Class.hpp NB.cpp NB.hpp main.cpp Genome.hpp Genome.cpp Diskutil.hpp Diskutil.cpp Read.cpp Read.hpp
	$(CC) $(OPT) Read.cpp Genome.cpp Class.hpp NB.cpp Diskutil.cpp main.cpp -o NB.run -lboost_filesystem -lboost_system -lboost_program_options
	chmod +x NB.run
clean:
	@rm NB.run
