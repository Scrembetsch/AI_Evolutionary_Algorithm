SRC_DIR := src
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,%.o,$(SRC_FILES))
LDFLAGS := -lm
CXXFLAGS := -Wall -Wextra -Werror -pedantic -O3

.DEFAULT_GOAL := MagicSquare

test: MagicSquare
	./MagicSquare --socket

MagicSquare: $(OBJ_FILES)
	g++ $(LDFLAGS) -o $@ $^

%.o: $(SRC_DIR)/%.cpp
	g++ $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm *.o MagicSquare