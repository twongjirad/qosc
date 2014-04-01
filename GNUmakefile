# -*- mode: Makefile -*-
# Don't hard code into this file

include Makefile.config

ifeq (${BASEDIR},)
	BASEDIR=${PWD}
endif

LIBS = rootvariable

.PHONY: basedir libs all

all: basedir libraries

lib%.so:
	@mkdir -p $(BASEDIR)/lib
	@echo 'Building library '$@,' from '$(subst .so,,$(subst lib,,$@))
	@cd $(subst .so,,$(subst lib,,$@)); gmake
#	@cd ${BASEDIR}

basedir:
	@echo 'Base directory set to ${BASEDIR}'

libraries: $(addprefix lib, $(addsuffix .so, $(LIBS)))
	@echo 'Made packages: '$^

