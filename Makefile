SHELL=/bin/sh

ifndef TARGET
  TARGET := $(error TARGET undefined)UNDEFINED
endif

all:
	cd src; $(MAKE) all

install:
	cd src; $(MAKE) install

clean:
	cd src; $(MAKE) clean
