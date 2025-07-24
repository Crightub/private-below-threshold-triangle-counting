CXX = g++-15
BOOST_ROOT = /usr/local/opt/boost
BOOST_INCLUDE = -I$(BOOST_ROOT)/include
INCLUDE = -Iinclude

CXXFLAGS = -std=c++17 -O2 -fopenmp $(BOOST_INCLUDE) $(INCLUDE)
LDFLAGS = -fopenmp

TARGET = graph_generator

SRCS = main.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(TEST_TARGET)
