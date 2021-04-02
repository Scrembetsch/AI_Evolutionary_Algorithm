SRC_DIR := src
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,%.o,$(SRC_FILES))
LDFLAGS := -lm
CXXFLAGS := -Wall -Wextra -Werror -pedantic -O3

.DEFAULT_GOAL := MagicSquare

test: MagicSquare
	./MagicSquare --mode 0 --size 3
	./MagicSquare --mode 0 --size 4
	./MagicSquare --mode 0 --size 5
	./MagicSquare --mode 1 --size 9 --m 3 --n 3
	./MagicSquare --mode 2 --size 8 --power 2

MagicSquare: $(OBJ_FILES)
	g++ $(LDFLAGS) -o $@ $^

%.o: $(SRC_DIR)/%.cpp
	g++ $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	-rm *.o MagicSquare