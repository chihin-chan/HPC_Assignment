CC=mpicxx
CXXFLAGS=-Wall -O0 -pedantic
LDLIBS=-llapack -lblas
OBJS = LidDrivenCavitySolver
TARGET=LidDrivenCavity
HDRS=LidDrivenCavity.h

default: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).o $(OBJS)

.PHONY: clean doc

doc:
	doxygen Doxyfile

clean:
	rm -rf $(TARGET) *.o html latex
