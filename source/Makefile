# TALYS Makefile
#
# Examples:
#
#   make
#   make FC=ifx
#   make FFLAGS="-O3 -march=native"
#   make FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"
#   make clean
#   make check

# GNU make has a built-in default FC=f77. Replace that built-in default
# with gfortran, while respecting settings from the environment or the
# command line.

ifeq ($(origin FC), default)
  FC = gfortran
endif

# Default compilation flags. Can be overridden on the command line or
# through the environment.

FFLAGS ?= -w -O3 -ffp-contract=off

# Installation directories.

TALYS_DIR := $(abspath ..)
BIN_DIR    ?= $(TALYS_DIR)/bin
TARGET     := $(BIN_DIR)/talys

# Source files.
#
# Keep fixed-form files first, matching the existing build procedure.
# Free-form files are sorted alphabetically; module files such as
# A0_talys_mod.f90 therefore appear before routines that use them.

FIXED_SOURCES := $(sort $(wildcard *.f))
FREE_SOURCES  := $(sort $(wildcard *.f90))
SOURCES       := $(FIXED_SOURCES) $(FREE_SOURCES)

.PHONY: all install clean check

all: $(TARGET)

$(TARGET): $(SOURCES)
	@echo "Compiler: $(FC)"
	@echo "Flags:    $(FFLAGS)"
	@mkdir -p "$(BIN_DIR)"
	$(FC) $(FFLAGS) -o "$(TARGET)" $(SOURCES)
	@rm -f *.o *.mod
	@echo "Installed: $(TARGET)"

install: all

clean:
	rm -f *.o *.mod talys "$(TARGET)"

check: all
	cd "$(TALYS_DIR)/samples" && ./verify
