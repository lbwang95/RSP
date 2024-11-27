CC = gcc
CXX = g++
LIBS = -lrt
CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 -mcmodel=large
TARGET = ./all

all: index index_PR TBS ERSP SDRSP SMOGA Maintenance
index: index.cpp cdfppf.h
	$(CXX) $(CFLAGS) -o index index.cpp cdfppf.h

index_PR: index_PR.cpp cdfppf.h
	$(CXX) $(CFLAGS) -o index_PR index_PR.cpp cdfppf.h
	
TBS: TBS.cpp cdfppf.h
	$(CXX) $(CFLAGS) -o TBS TBS.cpp cdfppf.h

ERSP: ERSP-Astar.cpp cdfppf.h
	$(CXX) $(CFLAGS) -o ERSP ERSP-Astar.cpp cdfppf.h

SDRSP: SDRSP-HAstar.cpp cdfppf.h
	$(CXX) $(CFLAGS) -o SDRSP SDRSP-HAstar.cpp cdfppf.h

SMOGA: SMOGA.cpp cdfppf.h
	$(CXX) $(CFLAGS) -o SMOGA SMOGA.cpp cdfppf.h
	
Maintenance: index_maintenance.cpp cdfppf.h
	$(CXX) $(CFLAGS) -o Maintenance index_maintenance.cpp cdfppf.h

.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm index 2> /dev/null || true
