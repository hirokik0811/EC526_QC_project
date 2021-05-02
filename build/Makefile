EIGEN = ../include/eigen/
VPATH = ../src/

CCX = g++
CXXFLAGS = -g -Wall -std=c++11 -I $(EIGEN)
#============================================================

TARGET1 = QFT2
C_OBJS1 = QFT2.o

TARGET2 = QFT_Benchmark
C_OBJS2 = QFT_Benchmark.o

TARGET3 = Ising_TE
C_OBJS3 = Ising_TE.o

C_SOURCES =  QFT2.cpp QFT_Benchmark.cpp Ising_TE.cpp
MY_INCLUDES = QCircuit.h matMul.h tensorProd.h

ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)

DATA_FILES = time_evolution_TFIM.txt Benchmark.txt
#============================================================
all: $(TARGET1) $(TARGET2) $(TARGET3)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<

$(TARGET1) :   $(C_OBJS1)
	$(CCX) $(CXXFLAGS)  $^ -o $@

$(TARGET2) :   $(C_OBJS2)
	$(CCX) $(CXXFLAGS)  $^ -o $@

$(TARGET3) :   $(C_OBJS3)
	$(CCX) $(CXXFLAGS)  $^ -o $@


# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================


clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3) $(C_OBJS1) $(C_OBJS2) $(C_OBJS3) $(DATA_FILES)*~