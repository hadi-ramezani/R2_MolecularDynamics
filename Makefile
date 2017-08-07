CC = icpc
CFLAGS  = -std=c++11 -Wall -O3 #-ipo #-g #-xavx #-g  
#LDFLAGS = -opt-report=4 -opt-report-phase=loop,vec -opt-report-file=stderr
HEAD = $(wildcard ./src/*.h)
SOURCES = $(wildcard ./src/*.cpp) 
EXECUTABLE = ./bin/MD
MKL = /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl

OBJECTS = $(SOURCES:.cpp=.o)
$(info Object Files = $(OBJECTS))
$(info header Files = $(HEAD))

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ -Wl,--start-group $(MKL)/lib/intel64/libmkl_intel_lp64.a $(MKL)/lib/intel64/libmkl_intel_thread.a $(MKL)/lib/intel64/libmkl_core.a -Wl,--end-group -L/$(MKL)/../compiler/lib/intel64 -liomp5 -lm -ldl -lpthread

%.o: %.cpp
	$(CC) $(CFLAGS) $(DEBUG) -c -o $@ $< 
clean:
	rm ./build/*.o $(EXECUTABLE)

trj:
	g++ $(DEBUG) -c -o ./build/Trajectory.o ./src/Trajectory.cpp
rand:
	icpc -std=c++11 -w -I$(MKL)/include -c -o Thermostat.o Thermostat.h -Wl,--start-group $(MKL)/lib/intel64/libmkl_intel_lp64.a $(MKL)/lib/intel64/libmkl_intel_thread.a $(MKL)/lib/intel64/libmkl_core.a -Wl,--end-group -L/$(MKL)/../compiler/lib/intel64 -liomp5 -lm -ldl -lpthread 
