# Compiler and flags
CXX ?= g++
CXXFLAGS = -O3 -std=c++17 -pthread -DNDEBUG -flto -ffast-math
LDFLAGS =

# On Windows, we might need to link against ws2_32 for some compatibility libraries
ifeq ($(OS),Windows_NT)
    CXXFLAGS += -msse4.2 -march=native
    LDFLAGS += -lws2_32
    TARGET = dedup.exe
else
    ifneq ($(shell uname -s),Darwin)
        CXXFLAGS += -msse4.2 -march=native
    endif
    TARGET = dedup
endif

# Source and target
SRC = main.cc

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Smoke tests
test: $(TARGET)
	./tests/smoke.sh

# Clean up
clean:
	rm -f $(TARGET)

.PHONY: all test clean
