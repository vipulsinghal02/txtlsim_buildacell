# Makefile for MATLAB TXTL modeling library
# RMM, 29 Sep 2012

VERSION = 0.2a
FILES = ChangeLog README *.csv txtl_init.m
SUBDIRS = components core doc examples 

build:
	mkdir -p txtl-$(VERSION)
	cp -prf $(FILES) $(SUBDIRS) txtl-$(VERSION)
	rm -rf `find txtl-$(VERSION) -name .svn`
	tar zcf txtl-$(VERSION).tgz txtl-$(VERSION)

