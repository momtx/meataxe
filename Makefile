# ------------------------------------------------------------------------------
# Build-time configuration
#
# Don't change the settings here. Instead, write your customized settings to
# Makefile.conf. This will make upgrading to a new version of the MeatAxe
# easier.
# ------------------------------------------------------------------------------

# Directory where binaries and run-time files will be installed.
MTXROOT = ${CURDIR}

# C compiler and common compiler/linker flags
# The default is to use GNU C
CC=gcc
CFLAGS1=-std=c99 -D_DEFAULT_SOURCE -g -O3 -Wall
LDFLAGS1=-g -Wall
# For debugging:
#CFLAGS1=-std=c99 -D_DEFAULT_SOURCE -g -Wall -Werror -DMTX_DEBUG
#CFLAGS1=-std=c99 -D_DEFAULT_SOURCE -g -Wall -Werror -DASM_MMX -DMTX_DEBUG

# Flags to pass to the ar utility
ARFLAGS1=
# For AIX 64-bit:
#ARFLAGS1=-X 64

# Select which kernel you want to use.
# Standard kernel, up to GF(256)
ZZZ=0
# Big kernel, up to GF(2^16)  -- NOTE: THIS IS NO LONGER AVAILABLE
#ZZZ=1

# Verbose output (echo all commands)
V=0

# configuration overrides
include Makefile.conf

# ------------------------------------------------------------------------------
# Other settings
# ------------------------------------------------------------------------------

SILENT0=@
SILENT1=
SILENT=${SILENT${V}}

MTXBIN = ${MTXROOT}/bin

CFLAGS=$(CFLAGS1) -I"${MTXROOT}/include" -Itmp -DMTX_ZZZ=${ZZZ}
LDFLAGS=$(LDFLAGS1)

PROGRAMS = \
  cfcomp chop decomp genmod mkcycl mkdotl mkgraph mkhom \
  mkinc mksub mktree orbrep precond pseudochop pwkond rad soc symnew tcond tuc \
  zad zbl zcf zcl zcp zct zcv zef zev zfr ziv zkd zmo zmu zmw znu zor zpo zpr \
  zpt zqt zro zsc zsi zsp zsy ztc zte ztm ztr zts zuk zvp

# ------------------------------------------------------------------------------
# Main targets
# ------------------------------------------------------------------------------

all build: $(PROGRAMS:%=${MTXBIN}/%) ${MTXBIN}/zzztest \
   ${MTXROOT}/lib/libmtx.a

clean:
	rm -rf bin tmp *.zzz check.ma1 check.pe1 check.po1

Makefile.conf:
	@echo "*** Creating empty Makefile.conf"
	@echo "*** See Makefile for more information"
	@sleep 5
	touch "$@"

# ------------------------------------------------------------------------------
# Compile C sources
# ------------------------------------------------------------------------------

tmp/%.o: src/%.c src/meataxe.h Makefile Makefile.conf
	@echo "# CC $*.c -> $@"
	${SILENT}mkdir -p tmp
	${SILENT}${CC} $(CFLAGS) -c src/$*.c -o $@

# ------------------------------------------------------------------------------
# Create directories
# ------------------------------------------------------------------------------

tmp/_mkdir:
	mkdir -p tmp
	touch "$@"

# ------------------------------------------------------------------------------
# Link programs
# ------------------------------------------------------------------------------

${MTXBIN}/%: tmp/%.o ${MTXROOT}/lib/libmtx.a
	@echo "# LD $@"
	${SILENT}mkdir -p "${MTXBIN}"
	${SILENT}$(CC) $(LDFLAGS) -o $@ tmp/$*.o "${MTXROOT}/lib/libmtx.a"

# ------------------------------------------------------------------------------
# Library
# ------------------------------------------------------------------------------

LIB_OBJS=\
	args berlekmp \
	bitstring \
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
	mfcore \
	mkendo\
	mmulscal mraddgen mrcore mrread mrtranspose mrwrite \
	msclean mscore \
	mtensor mtxobj os \
	permcmp permcore permdup perminv permmul permorder\
	permprint permpwr permread permwrite poladd\
	polcmp polcore polderive poldiv poldup\
	polgcd polmul polprint polread polwrite \
	quotient random rdcfgen \
	saction \
	spinup spinup2 \
	split stabpwr stfcore \
	stfread stfwrite \
	string \
	sumint \
	temap \
	tkinfo types vec2mat \
	wgen \
	zcleanrow zcmprow zgap zpermrow \
	zzz2 \
	version

${MTXROOT}/lib/libmtx.a: $(LIB_OBJS:%=tmp/%.o)
	@echo "# AR $@"
	${SILENT}mkdir -p "${MTXROOT}/lib"
	${SILENT}rm -f "$@"
	${SILENT}ar ${ARFLAGS1} r "$@" $(LIB_OBJS:%=tmp/%.o)

# ------------------------------------------------------------------------------
# Test suite
# ------------------------------------------------------------------------------

TESTS=\
  0001-zzztest \
  0010 \
  0020 \
  0100 \
  0105 \
  0106 \
  0107 \
  0108 \
  0109 \
  0110 \
  0111 \
  0112 \
  0113 \
  0114 \
  0115 \
  0116 \
  0117 \
  0118 \
  0119 \
  0120 \
  0121 \
  0200 \
  0201 \
  0202 \
  0210 \
  0211 \
  0212 \
  0213 \
  0214 \
  0215



test: $(TESTS:%=tmp/test-%.done)

# meataxe library tests

TS_OBJS1=c-args c-bitstring c-cfinfo c-charpol\
	c-ffio c-fileio c-ffmat c-ffrow c-fpoly \
	c-grease c-kernel c-matins c-matrix c-matset\
	c-os c-perm c-poly c-pseed c-quot c-random \
	c-stf c-tensor c-zzz

TS_OBJS=$(TS_OBJS1:%=tmp/%.o) ${MTXROOT}/lib/libmtx.a

${MTXBIN}/zzztest: $(TS_OBJS) tmp/_mkdir
	@echo "# LD $@"
	${SILENT}mkdir -p "${MTXBIN}"
	${SILENT}$(CC) $(LDFLAGS) -o "$@" $(TS_OBJS)

tmp/test_table.c: tmp/_mkdir tmp/tex $(TS_OBJS1:%=tests/%.c)
	tmp/tex $(TS_OBJS1:%=tests/%.c) >$@

tmp/tex: tests/tex.c tmp/_mkdir
	@$(CC) -Itests -Isrc $(CFLAGS) $(LDFLAGS) $< -o $@

tmp/c-%.o: tests/c-%.c tests/testing.h src/meataxe.h Makefile Makefile.conf
	@echo "# CC c-$*.c -> $@"
	${SILENT}mkdir -p tmp
	${SILENT}$(CC) $(CFLAGS) -Itests -Isrc -c "tests/c-$*.c" -o "$@"

tmp/c-zzz.o: tests/c-zzz.c tests/testing.h src/meataxe.h tmp/test_table.c Makefile Makefile.conf
	@echo "# CC tests/c-zzz.c -> $@"
	${SILENT}$(CC) $(CFLAGS) -Itests -Isrc -c "tests/c-zzz.c" -o "$@"

tmp/test-%.done: tests/common.sh ${MTXBIN}/zzztest tests/%/run \
   $(PROGRAMS:%=${MTXBIN}/%) tmp/_mkdir
	@cd tmp && MTXBIN="${MTXBIN}" MTXLIB="${MTXROOT}/lib" MTX_TEST_ID="$*" \
           MTX_ZZZ="${ZZZ}" \
	   PATH="${MTXBIN}:/usr/bin:/bin" ../tests/$*/run
	@touch "$@"


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

doc ${docProducts}: \
   etc/Doxyfile etc/layout.xml $(PROGRAMS:%=src/%.c) $(LIB_OBJS:%=src/%.c) \
   src/meataxe.h  \
   ${docDocs} src/meataxe.doc src/changelog.doc
	cp etc/Doxyfile tmp/Doxyfile.auto
	echo "PROJECT_NUMBER=${MTXVERSION}" >>tmp/Doxyfile.auto
	echo "OUTPUT_DIRECTORY=${docDir}" >>tmp/Doxyfile.auto
	mkdir -p ${docDir}
	doxygen tmp/Doxyfile.auto >tmp/doxygen.log

doxy:
	cp etc/Doxyfile tmp/Doxyfile.auto
	echo "PROJECT_NUMBER=${MTXVERSION}" >>tmp/Doxyfile.auto
	echo "OUTPUT_DIRECTORY=${docDir}" >>tmp/Doxyfile.auto
	mkdir -p ${docDir}
	doxygen tmp/Doxyfile.auto >tmp/doxygen.log


# ------------------------------------------------------------------------------
# Releasing
# ------------------------------------------------------------------------------

MTXVERSION=$(shell grep "^#define MTX_VERSION" src/meataxe.h | cut -d '"' -f 2)

EXPORTED_FILES =\
  $(PROGRAMS:%=src/%.c)\
  $(LIB_OBJS:%=src/%.c)\
  src/meataxe.h src/genconfig.c\
  Makefile README.md COPYING\
  src/meataxe.doc src/changelog.doc\

tar: all doc
	rm -f meataxe-${MTXVERSION} meataxe-${MTXVERSION}.tar meataxe-${MTXVERSION}.tar.gz \
	&& ln -s . meataxe-${MTXVERSION} \
	&& tar cf meataxe-${MTXVERSION}.tar $(EXPORTED_FILES:%=meataxe-${MTXVERSION}/%) \
	&& rm meataxe-${MTXVERSION} \
	&& gzip meataxe-${MTXVERSION}.tar \
	&& echo "Created meataxe-${MTXVERSION}.tar.gz"
