# -*- mode: Makefile -*-
include ../Makefile.inc

LIBNAME = analysistools

DEPDIR = .dep
OBJDIR = .obj

LOCAL_INC += -I./include
LOCAL_INC += -I$(BASE)/rootvariable/include
LOCAL_LIBS += -L$(LIBDIR) -lrootvariable

LOCAL_INC += -I$(PYTHONINC)

CCSRCS = $(wildcard src/*.cc)
PYSRCS = $(wildcard py/*.pyx)
PYSRCDIR = py
LDPYLIBS = $(ROOT_LIBS) $(LOCAL_LIBS) $(PYTHONLIB)
VPATH += 

include ../Makefile.auto

