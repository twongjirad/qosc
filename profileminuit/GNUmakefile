# -*- mode: Makefile -*-
include ../Makefile.config

LIBNAME = profileminuit

VPATH = src/

BASEDIR = ..
OBJDIR = .obj
DEPDIR = .dep

LOCAL_INC += -I./include
LOCAL_INC += -I$(BASEDIR)/analysistools/include -I$(BASEDIR)/analysisfw/include -I$(BASEDIR)/rootvariable/include
LOCAL_INC += -I$(PYTHONINC)

CCSRCS = $(wildcard src/*.cc)

OPT = -g -pg
#OPT = -O2

EXE_LIBS += -L$(LIBDIR) -lanalysistools -lanalysisfw -lrootvariable

EXESRCS = $(wildcard tests/*.cc)
EXEDIR = ./bin
VPATH += exesrc

#CC_FLAGS += -fno-omit-frame-pointer

include ../Makefile.auto

