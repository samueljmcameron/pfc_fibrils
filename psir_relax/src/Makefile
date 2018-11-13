# ------------------------------------------------
# Generic Makefile
#
# Author: yanick.rochon@gmail.com
# Date  : 2011-08-10
#
# Changelog :
#   2010-11-05 - first version
#   2011-08-10 - added structure : sources, objects, binaries
#                thanks to http://stackoverflow.com/users/128940/beta
#   2017-04-24 - changed order of linker params
# ------------------------------------------------

# project name (generate executable with this name)
TARGET   = psivsr

CC       = gcc
# compiling flags here
CFLAGS   = -std=c99 -I.

LINKER   = gcc
# linking flags here
LFLAGS   = -Wall -I. -lm -lgsl -lgslcblas

SRCDIR   = ../../src
OBJDIR   = ../obj
BINDIR   = ../bin

SOURCES  := $(wildcard $(SRCDIR)/*.c)
MAINSRC  := $(wildcard *.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
MAINOBJECT := $(MAINSRC:%.c=$(OBJDIR)/%.o)
rm = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS) $(MAINOBJECT)
	@$(LINKER) $(OBJECTS) $(MAINOBJECT) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(MAINOBJECT) : $(MAINSRC) 
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"


.PHONY: clean
clean:
	@$(rm) $(OBJECTS) $(MAINOBJECT)
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"

print-%	: ; @echo $* = $($*)
