# TALYS Makefile
#
# Examples:
#
#   make
#   make FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"
#   make FC=ifx FFLAGS="-O3"
#   make clean
#   make check

# Compiler-specific flags can be supplied by the user.
ifeq ($(origin FC), default)
  FC = gfortran
endif

ifeq ($(origin FFLAGS), undefined)
  ifneq ($(findstring gfortran,$(FC)),)
    FFLAGS = -w -O3 -ffp-contract=off
  else
    FFLAGS =
  endif
endif

# Installation directories.

TALYS_DIR := $(abspath ..)
BIN_DIR  ?= $(TALYS_DIR)/bin
TARGET   := $(BIN_DIR)/talys

# TALYS source files.
#
# Free-form source files must be named so that alphabetical order places
# module definitions, such as A0_talys_mod.f90, before their users.

FIXED_SOURCES := $(sort $(wildcard *.f))
FREE_SOURCES  := $(sort $(wildcard *.f90))
SOURCES       := $(FIXED_SOURCES) $(FREE_SOURCES)

.PHONY: all install clean check

all install: $(TARGET)

$(TARGET): $(SOURCES)
	@echo "Compiler: $(FC)"
	@echo "Flags:    $(FFLAGS)"
	@mkdir -p "$(BIN_DIR)"
	$(FC) $(FFLAGS) -o "$(TARGET)" $(SOURCES)
	@rm -f *.o *.mod
	@echo "Installed: $(TARGET)"

clean:
	rm -f *.o *.mod talys "$(TARGET)"

check: all
	cd "$(TALYS_DIR)/samples" && \
	  TALYS_DIR="$(TALYS_DIR)" ./verify
