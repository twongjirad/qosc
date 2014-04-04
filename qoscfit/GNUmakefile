# -*- mode: Makefile -*-
include ../Makefile.config

LIBNAME = qoscfit
EXEDIR = ./bin
OBJDIR = .obj
DEPDIR = .dep

ifeq ($(BASEDIR),)
	BASEDIR = ../
endif

# ----------------------------------------------------
# SETUP DEPENDENCIES

LOCAL_INC += -I./include

# qosc Libraries
LOCAL_INC += -I$(BASEDIR)/rootvariable/include/ -I$(BASEDIR)/analysistools/include
LOCAL_INC += -I$(BASEDIR)/profileminuit/include
LOCAL_INC += -I$(BASEDIR)/analysisfw/include
LOCAL_LIBS += -L$(LIBDIR) -lrootvariable -lanalysistools -lanalysisfw -lprofileminuit

# T2K Reweight Headers: systematic uncertainty package for reweighting MC predictions
# Can I remove this from the main package? This is application specific.
ifneq ($(T2KRW_VER),0)
T2KREWEIGHT_INC = -I$(T2KREWEIGHT)/src -I$(T2KREWEIGHT)/NIWGReWeight -I$(T2KREWEIGHT)/JReWeight  -I$(T2KREWEIGHT)/SKDetReWeight
T2KREWEIGHT_INC += -I$(NEUT_ROOT)/src/reweight -I$(NEUT_ROOT)/include -I$(NEUT_ROOT)/inc -I$(NEUT_ROOT)/src/neutclass
LOCAL_INC += $(T2KREWEIGHT_INC)
LOCAL_INC += -I$(BASEDIR)/syserrs2012/include
LOCAL_LIBS += -lodsyserrors2012
endif

# Google Performance Tools
ifeq ($(USE_GPERFTOOLS),1)
LOCAL_INC += -I$(GPERFTOOLS)/include
endif

# Prob3++: Oscillation Probability Calculator
LOCAL_INC += -I$(PROBTPP_DIR)
LOCAL_LIBS += $(PROBTPP_LIB)

# T2K Reweight Libraries
ifneq ($(T2KRW_VER),0)
LOCAL_LIBS += -L$(T2KREWEIGHT)/lib -lT2KReWeight -L$(T2KREWEIGHT)/JReWeight -lJReWeight -L$(T2KREWEIGHT)/NIWGReWeight -lNIWGReWeight -L$(T2KREWEIGHT)/SKDetReWeight -lSKDetReWeight
LOCAL_LIBS += -L$(NEUT_ROOT)/lib -ltauola -lneutcore -lnuccorrspl -lnuceff -lpartnuck -lskmcsvc -lneutclass 
LOCAL_LIBS += -L$(NEUT_ROOT)/src/reweight_head -lNReWeight ${CERNLIBS} -gfortran
LOCAL_LIBS += `cernlib -pdflib804`
endif

# CERNLIB
LOCAL_LIBS += $(CERNLIBS)

# PYTHON BINDINGS
# Python Bindings via Cython
ifeq ($(BUILD_CYTHON_LIBRARIES),1)
LOCAL_LIBS += $(PYTHONLIB)
LOCAL_INC += -I$(PYTHONINC)
PYSRCS = $(wildcard py/*.pyx)
PYSRCDIR = py
LDPYLIBS = $(ROOT_LIBS) $(LOCAL_LIBS)
#LDPYLIBS += $(SKOFL_ROOT)/lib/libzbs.a $(CERNLIBS) -lgfortran -lstdc++
endif

EXE_LIBS += ${LOCAL_LIBS} ${ROOT_LIBS} ${CERNLIBS} 

CCSRCS = $(wildcard src/*.cc)
EXESRCS = $(wildcard exesrc/*.cc)
EXESRCS += $(wildcard tests/*.cc)
VPATH += exesrc

include ../Makefile.auto

