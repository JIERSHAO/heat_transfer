## Option
CXX = g++
CXXFLAGS = -Wall -O2

TARGET = grid
SRCS = gridmain.cpp transform.cpp
OBJS = gridmain.o
RES = grid.dat

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^
clean:
	rm -f $(OBJS) $(TARGET) $(RES)

.PHONY: all clean