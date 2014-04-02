# -*- mode: Makefile -*-
include ../Makefile.config

LIBNAME = analysistools

BASEDIR = ..
DEPDIR = .dep
OBJDIR = .obj

LOCAL_INC += -I./include
LOCAL_INC += -I$(BASEDIR)/rootvariable/include
LOCAL_LIBS += -L$(LIBDIR) -lrootvariable

LOCAL_INC += -I$(PYTHONINC)

CCSRCS = $(wildcard src/*.cc)
PYSRCS = $(wildcard py/*.pyx)
PYSRCDIR = py
LDPYLIBS = $(ROOT_LIBS) $(LOCAL_LIBS) $(PYTHONLIB)
VPATH += 

include ../Makefile.auto

