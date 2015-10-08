
# ------------------------------------------------------------------------------
# Build-time configuration
#
# Note 1: don't change the settings below. Instead, write you customized
# settings to Makefile.conf. This will make upgrading to a new version of
# the MeatAxe easier.
#
# Note 2: even if you use the default settings you must create an empty
# Makefile.conf to build the MeatAxe.
# ------------------------------------------------------------------------------

# Directory where binaries and run-time files will be installed.
MTXROOT = ${CURDIR}

# Directory where executables are installed.
# This must be an absolute path (starting with '/').
MTXBIN=${MTXROOT}/bin

# MeatAxe library directory (obsolete, not used anymore)
MTXLIB=${MTXROOT}/lib

# C compiler and common compiler/linker flags
# The default is to use GNU C
CC=gcc	
CFLAGS1=-std=c99 -g -O3 -Wall
LDFLAGS1=-g -Wall
# For debugging:
#CFLAGS1=-std=c99 -g -Wall -Werror -DDEBUG -DPARANOID
#CFLAGS1=-std=c99 -g -Wall -Werror -DASM_MMX -DDEBUG

# Select which kernel you want to use.
# Standard kernel, up to GF(256)
ZZZ=0
# Big kernel, up to GF(2^16)  -- NOTE: THIS IS NO LONGER AVAILABLE
#ZZZ=1

# configuration overrides
include Makefile.conf

# ------------------------------------------------------------------------------
# Other settings
# ------------------------------------------------------------------------------

V=0
SILENT0=@
SILENT1=
SILENT=${SILENT${V}}

MTXVERSION = 2.5.0-SNAPSHOT

CFLAGS=$(CFLAGS1) -I"${MTXROOT}/include" -Itmp -DZZZ=${ZZZ}
LDFLAGS=$(LDFLAGS1)

PROGRAMS = \
  cfcomp checksum chop decomp genmod mkcycl mkdotl mkgraph mkhom mkhom_old\
  mkinc mksub mktree orbrep precond pseudochop pwkond rad soc symnew tcond tuc \
  zad zbl zcf zcl zcp zct zcv zef zev zfr ziv zkd zmo zmu zmw znu zor zpo zpr \
  zpt zqt zro zsc zsi zsp zsy ztc zte ztm ztr zts zuk zvp

# ------------------------------------------------------------------------------
# Main targets
# ------------------------------------------------------------------------------

all build: $(PROGRAMS:%=${MTXBIN}/%) ${MTXROOT}/lib/libmtx.a ${MTXROOT}/include/meataxe.h

clean:
	rm -rf bin tmp *.zzz check.ma1 check.pe1 check.po1
	-cd tests; make clean

Makefile.conf:
	@echo "------"
	@echo "Create Makefile.conf (which may be empty) and try again"
	@echo "Read the Makefile for more information"
	@echo "------"
	@false

# ------------------------------------------------------------------------------
# Compile C sources
# ------------------------------------------------------------------------------

tmp/%.o: tmp/mk.dir src/%.c ${MTXROOT}/include/meataxe.h
	${SILENT}echo "CC $*.c -> $@"
	${SILENT}${CC} $(CFLAGS) -c src/$*.c -o $@

tmp/mk.dir:
	mkdir -p tmp
	touch $@

${MTXBIN}/mk.dir:
	mkdir -p "${MTXBIN}"
	touch $@

# ------------------------------------------------------------------------------
# Link programs
# ------------------------------------------------------------------------------

${MTXBIN}/%: ${MTXBIN}/mk.dir tmp/%.o ${MTXROOT}/lib/libmtx.a
	${SILENT}mkdir -p "${MTXBIN}"
	${SILENT}echo "LD $@"
	${SILENT}$(CC) $(LFLAGS) -o $@ tmp/$*.o "${MTXROOT}/lib/libmtx.a"


# ------------------------------------------------------------------------------
# Library
# ------------------------------------------------------------------------------

LIB_OBJS=\
	args berlekmp \
	bsand bscore bsdup bsissub bsmatch bsminus \
	bsop bsor bsprint bsread bswrite \
	cfinfo \
	charpol chbasis \
	error \
	ffio \
	fpcore fpdup fpmul fpmul2 fpprint \
	gcd genseed\
	grmaprow grmatcore grtable \
	homcomp \
	imatcore imatread imatwrite\
	init intio issub \
	isisom kernel-$(ZZZ) \
	ldiag \
	maddmul mat2vec matadd matclean matcmp \
	maketab-$(ZZZ) \
	matcopy matcore matcut \
	matdup matech matid matins matinv matmul \
	matnull matorder \
	matpivot\
	matprint matpwr matread mattr mattrace matwrite \
	message \
	mfcore mfread mfreadlong mfwrite mfwritelong \
	minpol mkendo\
	mmulscal mraddgen mrcore mrread mrtranspose mrwrite \
	msclean mscore \
	mtensor mtxobj os \
	permcmp permcore permdup perminv permmul permorder\
	permprint permpwr permread permwrite poladd\
	polcmp polcore polderive poldiv poldup\
	polgcd polmul polprint polread polwrite \
	quotient random rdcfgen \
	saction setcore setinsert settest \
	spinup spinup2 \
	split stabpwr stfcore \
	stfread stfwrite \
	string \
	sumint \
	temap \
	tkinfo vec2mat \
	wgen \
	zcleanrow zcmprow zgap zpermrow \
	zzz2 \
	version

${MTXROOT}/lib/libmtx.a: $(LIB_OBJS:%=tmp/%.o)
	${SILENT}echo "AR $@"
	${SILENT}mkdir -p "${MTXROOT}/lib"
	${SILENT}rm -f "$@"
	${SILENT}ar r "$@" $(LIB_OBJS:%=tmp/%.o)

${MTXROOT}/include/meataxe.h: src/meataxe.h.in tmp/genconfig
	${SILENT}echo "CF src/meataxe.h.in -> $@"
	${SILENT}mkdir -p "${MTXROOT}/include"
	${SILENT}tmp/genconfig <src/meataxe.h.in >"$@"

tmp/genconfig: Makefile Makefile.conf src/genconfig.c
	${SILENT}echo "CL src/genconfig.c -> $@"
	${SILENT}mkdir -p tmp
	${SILENT}$(CC) $(CFLAGS) -DZZZ=${ZZZ} -DMTXVERSION="${MTXVERSION} \
	   "$(LFLAGS) -o "$@" src/genconfig.c
   

# ------------------------------------------------------------------------------
# Test suite
# ------------------------------------------------------------------------------

TS_OBJS1=c-args c-bitstring c-charpol\
	c-ffio c-fileio c-ffmat c-ffrow c-fpoly \
	c-grease c-kernel c-matins c-matrix c-matset\
	c-os c-perm c-poly c-pseed c-quot c-random \
	c-sets c-stf c-tensor c-zzz

TS_OBJS=$(TS_OBJS1:%=tmp/%.o) ${MTXROOT}/lib/libmtx.a

tmp/c-%.o: tests/c-%.c
	${SILENT}$(CC) -c $(CFLAGS) -o "$@" "$<"

${MTXBIN}/zzztest: $(TS_OBJS)
	${SILENT}echo "LD $@"
	${SILENT}mkdir -p "${MTXBIN}"
	${SILENT}$(CC) $(LDFLAGS) -o "$@" $(TS_OBJS)

${MTXBIN}/checksum: tmp/checksum.o
	${SILENT}mkdir -p "${MTXBIN}"
	${SILENT}echo "LD $@"
	${SILENT}$(CC) $(CFLAGS) -o "$@" tmp/checksum.o

TESTS=0000\
  0100 0100 0105 0106 0107 0108 0109 0110 0111 0112\
  0113 0114 0115 0116 0117 0118 0200 0201 0210 0211\
  0212 0213 0214 0215

#tests: $(TESTS:%=tmp/t-%.done)
test: tmp/zzztest.done 

tmp/test_table.c: tmp/tex $(TS_OBJS1:%=tests/%.c)
	tmp/tex $(TS_OBJS1:%=tests/%.c) >$@

tmp/tex: tests/tex.c
	@$(CC) -Itests -Isrc $(CFLAGS) $(LFLAGS) $< -o $@

tmp/c-%.o: tests/c-%.c src/meataxe.h tmp/config.h
	${SILENT}mkdir -p tmp
	${SILENT}echo "CC $*.c -> $@"
	${SILENT}$(CC) $(CFLAGS) -Itests -Isrc -c "tests/c-$*.c" -o "$@"

tmp/c-zzz.o: tests/c-zzz.c tmp/test_table.c
	${SILENT}$(CC) $(CFLAGS) -Itests -Isrc -c tests/c-zzz.c -o $@

tmp/zzztest.done: ${MTXBIN}/zzztest
	mkdir -p tmp
	cd tmp && ${MTXBIN}/zzztest
	touch $@

tmp/t-%.done: tmp/mk.dir tests/t-% tmp/t.config ${MTXBIN}/checksum build
	@echo "t-$* `grep '^#:' tests/t-$* | cut -c 3-100`"
	@cd tmp && ../tests/t-$*
	@touch $@

tmp/t.config: tmp/mk.dir tests/config
	cp tests/config $@

.PHONY: tar clean install test
.PRECIOUS: tmp/%.o


# ------------------------------------------------------------------------------
# Package maintenance
# ------------------------------------------------------------------------------

info:
	@echo "TS_OBJS=${TS_OBJS}"


# ------------------------------------------------------------------------------
# Documentation (requires Doxygen)
# ------------------------------------------------------------------------------

.PHONY: doc

docDir = doc/${MTXVERSION}
docDocs = src/changelog.doc src/meataxe.doc src/sections.doc
docProducts = ${docDir}/index.html ${docDir}/pages.html ${docDir}/classes.html

doc: ${docProducts}

${docProducts}: \
   etc/Doxyfile etc/layout.xml $(PROGRAMS:%=src/%.c) $(LIB_OBJS:%=src/%.c) \
   ${MTXROOT}/include/meataxe.h  \
   ${docDocs} src/meataxe.doc src/changelog.doc
	mkdir -p ${docDir}
	doxygen etc/Doxyfile


# ------------------------------------------------------------------------------
# Releasing
# ------------------------------------------------------------------------------

EXPORTED_FILES =\
  $(PROGRAMS:%=src/%.c)\
  $(LIB_OBJS:%=src/%.c)\
  src/meataxe.h.in src/genconfig.c\
  Makefile README.md COPYING\
  src/meataxe.doc src/changelog.doc\

  #tests/check.h $(TS_OBJS1:%=tests/%.c) $(TS_OBJS1:%=tests/%.h) \

tar: all doc
	rm -f meataxe-${MTXVERSION} meataxe-${MTXVERSION}.tar meataxe-${MTXVERSION}.tar.gz \
	&& ln -s . meataxe-${MTXVERSION} \
	&& tar cf meataxe-${MTXVERSION}.tar $(EXPORTED_FILES:%=meataxe-${MTXVERSION}/%) \
	&& rm meataxe-${MTXVERSION} \
	&& gzip meataxe-${MTXVERSION}.tar \
	&& echo "Created meataxe-${MTXVERSION}.tar.gz"

