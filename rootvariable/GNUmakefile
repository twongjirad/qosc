# -*- mode:Makefile
include ../Makefile.config

DEPDIR = .dep
OBJDIR = .obj

BASEDIR = ..
LIBNAME = rootvariable
LOCAL_INC += -I./include
CCSRCS = $(wildcard src/*.cc)
EXE_LIBS += $(LIB)
EXESRCS += $(wildcard tests/*.cc)
VPATH += exesrc

include ../Makefile.auto