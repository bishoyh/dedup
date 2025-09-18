# Compiler and flags
CXX ?= g++
CXXFLAGS = -O3 -std=c++17 -pthread -DNDEBUG -flto -ffast-math -msse4.2
LDFLAGS =

# On Windows, we might need to link against ws2_32 for some compatibility libraries
ifeq ($(OS),Windows_NT)
    LDFLAGS += -lws2_32
    TARGET = dedup.exe
else
    TARGET = dedup
endif

# Source and target
SRC = main.cc

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Clean up
clean:
	rm -f $(TARGET)

.PHONY: all clean
