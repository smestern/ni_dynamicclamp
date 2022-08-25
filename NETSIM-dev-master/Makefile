# ================================================================================
# ================================================================================
# ================================================================================
#
#   Makefile: NETSIM
#
#   Copyright (C) 2018-2020 Lyle Muller
#	http://mullerlab.ca
#
#	COMMAND LINE OPTIONS:
#
# 	EXTERNAL_INPUT=yes		compile in support for external Poisson noise
# 	RELEASE_PROBABILITY=yes		for probabilistic release
#
# ================================================================================
# ================================================================================
# ================================================================================


SRCDIR = src
BUILDDIR = ./obj

#
# compiler settings
#

CC = gcc

#
# complation type
#

CFLAGS = -std=gnu11
LDFLAGS = 

# debug
ifeq ($(COMPILE_TYPE),debug)
    $(info Using performance flags)
    CFLAGS += -ggdb -Wall
endif

# profiler
ifeq ($(COMPILE_TYPE),profiler)
    $(info Using profiler flags)
    CFLAGS += -g
endif

# 
ifeq ($(COMPILE_TYPE),profilegen)
    $(info Using profile generate flags)
    CFLAGS += -O3 -Wall -march=native -flto -fprofile-generate -fprofile-arcs -ftest-coverage
    LDFLAGS += -fprofile-arcs
endif

ifeq ($(COMPILE_TYPE),profileuse)
    $(info Using profile use flags)
    CFLAGS += -O3 -Wall -march=native -flto -fprofile-use
endif

# performance
ifeq ($(COMPILE_TYPE),performance)
    $(info Using performance flags)
    CFLAGS += -O3 -Wall -march=native -flto
endif

#
# differential compile options
#

ifeq ($(EXTERNAL_INPUT),yes)
    $(info Using external input)
    FLAGS += -DEXTERNAL_INPUT
endif

ifeq ($(RELEASE_PROBABILITY),yes)
    $(info Using releae probability)
    FLAGS += -DRELEASE_PROBABILITY
endif

#
# library settings
#

ifeq ($(shell uname),Linux)
	LIBS += -lm  
endif

LIBS += libnidaqmx.so
#PACK IN THE NI INTERFACE?




#
# files
#

SRC = $(wildcard $(SRCDIR)/*.c)
INC = $(wildcard $(SRCDIR)/*.h)





OBJS = $(patsubst %.c, %.o, $(filter %.c, $(subst $(SRCDIR), $(BUILDDIR), $(SRC))))

TARGET = netsim

#
# rules
#

.PHONY: build
build: $(BUILDDIR) $(OBJS) 
	@$(CC) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)
	@echo [LD] Linked $^ into $(TARGET)


$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	@$(CC) $(CFLAGS) -c $^ -o $@
	@echo [CC] Compiled $^ into $@

.PHONY: clean
clean:
	@rm -f $(OBJS) $(TARGET)
	@echo Cleaned $(OBJS) and $(TARGET)

.PHONY: rebuild
rebuild: clean build

$(BUILDDIR):
	mkdir -p $@
