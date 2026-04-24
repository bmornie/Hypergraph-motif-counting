CC      = g++
CFLAGS  = -std=c++23 -O3 -Wall -Wextra -flto=auto -march=native -DNDEBUG

TARGET1 = motif_counting
TARGET2 = compute_abundance

SRC1 := motif_counting.cpp hypergraph_motif.cpp
SRC2 := compute_abundance.cpp random_models.cpp hypergraph_motif.cpp

OBJ1 := $(SRC1:.cpp=.o)
OBJ2 := $(SRC2:.cpp=.o)

DEP := $(OBJ1:.o=.d) $(OBJ2:.o=.d)

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJ1)
	$(CC) $(CFLAGS) -o $@ $^

$(TARGET2): $(OBJ2)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

-include $(DEP)

clean:
	rm -f $(OBJ1) $(OBJ2) $(DEP) $(TARGET1) $(TARGET2)

run: $(TARGET1)
	./$(TARGET1)

run2: $(TARGET2)
	./$(TARGET2)

.PHONY: all clean run run2