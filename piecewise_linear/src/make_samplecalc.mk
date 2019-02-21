CC = gcc
CFLAGS = -I. -O2
LINKER = gcc
LFLAGS = -Wall -I. -O2 -lm -lgsl -lgslcblas
TARGET = samplecalc

BINDIR = ../bin
GSL_SRCDIR = edited_gsl_src
ENERGY_SRCDIR = funcs
OBJDIR = ../obj

LOCAL_SRC = samplecalc.c

ENERGY_SRC := energyderivs.c energyfunc.c f1_functions.c f2_functions.c \
              g1_functions.c g2_functions.c u_functions.c v_functions.c \
              nrutil.c

EDITED_GSL_SRC := $(wildcard $(GSL_SRCDIR)/*.c)

INCLUDES := headerfile.h
EDITEDGSL_INC := $(wildcard $(GSL_SRCDIR)/*.h)


EDITED_GSL_OBJS := $(EDITED_GSL_SRC:$(GSL_SRCDIR)/%.c=$(OBJDIR)/%.o)
ENERGY_OBJS := $(ENERGY_SRC:%.c=$(OBJDIR)/%.o)
LOCAL_OBJS := $(LOCAL_SRC:%.c=$(OBJDIR)/%.o)
OBJECTS := $(ENERGY_OBJS) $(EDITED_GSL_OBJS) $(LOCAL_OBJS)
rm = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $(OBJECTS) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(LOCAL_OBJS) : $(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(ENERGY_OBJS) : $(OBJDIR)/%.o: $(ENERGY_SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(EDITED_GSL_OBJS) : $(OBJDIR)/%.o: $(GSL_SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

PHONY: clean
clean:
	$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"
