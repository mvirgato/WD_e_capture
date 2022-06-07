#FILES = README.md COPYING.md pcubature.c hcubature.c cubature.h clencurt.h vwrapper.h converged.h test.c clencurt_gen.c NEWS.md

# CFLAGS = -pg -O3 -fno-inline-small-functions -Wall -ansi -pedantic
# CFLAGS = -g -Wall -ansi -pedantic
# CFLAGS = -O3 -Wall #-ansi -pedantic

# CC=gcc-9

# CFLAGS =-O3 -fopenmp

# LDFLAGS =-lgsl -lgslcblas -lm

# SOURCES=multiscatter.c finite_T.c crate_apprx_funcs.c cap_interp.c cap_funcs.c elec_wd_main.c
# OBJECTS=$(SOURCES:.c=.o)
# EXEC=run_capelec


# all: $(EXEC)

# $(EXEC): $(OBJECTS)
# 	$(CC) $(CFLAGS)  -o $@ $^  $(LDFLAGS)
# 	rm *.o

# .o:
# 	$(CC) $(CFLAGS) $@ -c $<

# clean:
# 	rm *.o $(EXEC)


CC = gcc-9
SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin

EXE := $(BIN_DIR)/run_capelec
SRC := $(wildcard $(SRC_DIR)/*.c)
OBJ := $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

INC = -Iinclude
CXXFLAGS :=  -MMD -MP
CFLAGS   := -O3 -fopenmp
LDFLAGS  := -Llib
LDLIBS   := -lgsl -lgslcblas -lm

.PHONY: all clean

all: $(EXE)

$(EXE):	$(OBJ) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(OBJ_DIR)/%.o:	$(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INC) $(CXXFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)

-include $(OBJ:.o=.d)
