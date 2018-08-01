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
TARGET   = EvsL

CC       = gcc
# compiling flags here
CFLAGS   = -std=c99 -I.

LINKER   = gcc
# linking flags here
LFLAGS   = -Wall -I. -lm

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

#$(OBJDIR)/bksub.o : ../../src/bksub.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/difeq.o : ../../src/difeq.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/energy.o : ../../src/energy.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/energylandscape.o : ../../src/energylandscape.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/nrutil.o : ../../src/nrutil.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/pinvs.o : ../../src/pinvs.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/polint.o : ../../src/polint.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/qromb.o : ../../src/qromb.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/red.o : ../../src/red.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/solvde.o : ../../src/solvde.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

#$(OBJDIR)/trapzd.o : ../../src/trapzd.c 
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

# this is where I'm leaving off, to get things working fine
# uncomment everything above, but comment out the next three
# lines below.

$(OBJECTS): $(SOURCES)
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