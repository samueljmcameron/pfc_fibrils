CC = gcc
CFLAGS = -I. -O2
LINKER = gcc
LFLAGS = -Wall -I. -O2 -lm -lgsl -lgslcblas
TARGET = driver

BINDIR = ../bin
GSL_SRCDIR = edited_gsl_src
SHARED_SRCDIR = ../shared_src
OBJDIR = ../obj

LOCAL_SRC = driver.c

SHARED_SRC := pinvs.c red.c shooting.c bksub.c energy.c nrutil.c polint.c trapzd.c \
              scaling.c solvde.c difeq.c finite_differences.c qromb.c shared.c

EDITED_GSL_SRC := $(wildcard $(GSL_SRCDIR)/*.c)

INCLUDES := headerfile.h
EDITEDGSL_INC := $(wildcard $(GSL_SRCDIR)/*.h)


EDITED_GSL_OBJS := $(EDITED_GSL_SRC:$(GSL_SRCDIR)/%.c=$(OBJDIR)/%.o)
SHARED_OBJS := $(SHARED_SRC:%.c=$(OBJDIR)/%.o)
LOCAL_OBJS := $(LOCAL_SRC:%.c=$(OBJDIR)/%.o)
OBJECTS := $(SHARED_OBJS) $(EDITED_GSL_OBJS) $(LOCAL_OBJS)
rm = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $(OBJECTS) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(LOCAL_OBJS) : $(LOCAL_SRC)
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(SHARED_OBJS) : $(OBJDIR)/%.o: $(SHARED_SRCDIR)/%.c
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
