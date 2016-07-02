# ------------------------------------------------------------------------------
# Build-time configuration
#
# Don't change the settings here. Instead, write your customized settings to
# Makefile.conf. This will make upgrading to a new version of the MeatAxe
# easier.
# ------------------------------------------------------------------------------

# Directory where binaries and run-time files will be installed.
MTXROOT = ${CURDIR}

# Directory where executables are installed.
# This must be an absolute path (starting with '/').
MTXBIN=${MTXROOT}/bin

# C compiler and common compiler/linker flags
# The default is to use GNU C
CC=gcc
CFLAGS1=-std=c99 -g -O3 -Wall
LDFLAGS1=-g -Wall
# For debugging:
#CFLAGS1=-std=c99 -g -Wall -Werror -DDEBUG -DPARANOID
#CFLAGS1=-std=c99 -g -Wall -Werror -DASM_MMX -DDEBUG

# Flags to pass to the ar utility
ARFLAGS1=
# For AIX 64-bit:
#ARFLAGS1=-X 64

# Select which kernel you want to use.
# Standard kernel, up to GF(256)
ZZZ=0
# Big kernel, up to GF(2^16)
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

MTXVERSION = 2.4.x-SNAPSHOT

CFLAGS=$(CFLAGS1) -I"${MTXROOT}/include" -Itmp \
 -DZZZ=${ZZZ} -DMTXLIB="${MTXROOT}/lib" -DMTXBIN="${MTXBIN}"
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

Makefile.conf:
	@echo "*** Creating empty Makefile.conf"
	@echo "*** See Makefile for more information"
	@sleep 5
	touch "$@"

# ------------------------------------------------------------------------------
# Compile C sources
# ------------------------------------------------------------------------------

tmp/%.o: src/%.c ${MTXROOT}/include/meataxe.h
	@echo "# CC $*.c -> $@"
	${SILENT}mkdir -p tmp
	${SILENT}${CC} $(CFLAGS) -c src/$*.c -o $@

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
	maketabF \
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
	@echo "# AR $@"
	${SILENT}mkdir -p "${MTXROOT}/lib"
	${SILENT}rm -f "$@"
	${SILENT}ar ${ARFLAGS1} r "$@" $(LIB_OBJS:%=tmp/%.o)

${MTXROOT}/include/meataxe.h: src/meataxe.h.in tmp/genconfig
	@echo "# GEN src/meataxe.h.in -> $@"
	${SILENT}mkdir -p "${MTXROOT}/include"
	${SILENT}tmp/genconfig <src/meataxe.h.in >"$@"

tmp/genconfig: Makefile Makefile.conf src/genconfig.c
	@echo "# CC src/genconfig.c -> $@"
	${SILENT}mkdir -p tmp
	${SILENT}$(CC) $(CFLAGS) -DZZZ=${ZZZ} -DMTXVERSION="${MTXVERSION} \
	   "$(LDFLAGS) -o "$@" src/genconfig.c
   

# ------------------------------------------------------------------------------
# Test suite
# ------------------------------------------------------------------------------

TESTS=0100 0105 0106 0107 0108 0109 0110 0111 0112 0113 0114 0115 0116 0117 0118 \
   0200 0201 0210 0211 0212 0213 0214 0215

test: tmp/zzztest.done $(TESTS:%=tmp/test-%.done)

# meataxe library tests

TS_OBJS1=c-args c-bitstring c-charpol\
	c-ffio c-fileio c-ffmat c-ffrow c-fpoly \
	c-grease c-kernel c-matins c-matrix c-matset\
	c-os c-perm c-poly c-pseed c-quot c-random \
	c-sets c-stf c-tensor c-zzz

TS_OBJS=$(TS_OBJS1:%=tmp/%.o) ${MTXROOT}/lib/libmtx.a

tmp/c-%.o: tests/c-%.c tests/check.h
	${SILENT}mkdir -p tmp
	${SILENT}$(CC) -c $(CFLAGS) -o "$@" "$<"

zzztest: ${MTXBIN}/zzztest
.PHONY: zzztest

${MTXBIN}/zzztest: $(TS_OBJS)
	@echo "# LD $@"
	${SILENT}mkdir -p "${MTXBIN}"
	${SILENT}$(CC) $(LDFLAGS) -o "$@" $(TS_OBJS)

${MTXBIN}/checksum: tmp/checksum.o
	${SILENT}mkdir -p "${MTXBIN}"
	@echo "# LD $@"
	${SILENT}$(CC) $(CFLAGS) -o "$@" tmp/checksum.o

tmp/test_table.c: tmp/tex $(TS_OBJS1:%=tests/%.c)
	tmp/tex $(TS_OBJS1:%=tests/%.c) >$@

tmp/tex: tests/tex.c
	@$(CC) -Itests -Isrc $(CFLAGS) $(LDFLAGS) $< -o $@

tmp/c-%.o: tests/c-%.c src/meataxe.h tmp/config.h
	${SILENT}mkdir -p tmp
	@echo "# CC $*.c -> $@"
	${SILENT}$(CC) $(CFLAGS) -Itests -Isrc -c "tests/c-$*.c" -o "$@"

tmp/c-zzz.o: tests/c-zzz.c tmp/test_table.c
	${SILENT}$(CC) $(CFLAGS) -Itests -Isrc -c tests/c-zzz.c -o $@

tmp/zzztest.done: ${MTXBIN}/zzztest
	mkdir -p tmp
	cd tmp && ${MTXBIN}/zzztest
	touch $@

# other tests

tmp/test-%.done: tests/common.sh tests/%/run $(PROGRAMS:%=${MTXBIN}/%)
	@echo "Running test $*"
	${SILENT}mkdir -p tmp
	@cd tmp && MTXBIN="${MTXBIN}" MTXLIB="${MTXROOT}/lib" MTX_TEST_ID="$*" \
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

doc: ${docProducts}

${docProducts}: \
   etc/Doxyfile etc/layout.xml $(PROGRAMS:%=src/%.c) $(LIB_OBJS:%=src/%.c) \
   ${MTXROOT}/include/meataxe.h  \
   ${docDocs} src/meataxe.doc src/changelog.doc
	cp etc/Doxyfile tmp/Doxyfile.auto
	echo "PROJECT_NUMBER=${MTXVERSION}" >>tmp/Doxyfile.auto
	echo "OUTPUT_DIRECTORY=${docDir}" >>tmp/Doxyfile.auto
	mkdir -p ${docDir}
	doxygen tmp/Doxyfile.auto >tmp/doxygen.log

# ------------------------------------------------------------------------------
# Releasing
# ------------------------------------------------------------------------------

EXPORTED_FILES =\
  $(PROGRAMS:%=src/%.c)\
  $(LIB_OBJS:%=src/%.c)\
  src/meataxe.h.in src/genconfig.c\
  Makefile README.md COPYING\
  src/meataxe.doc src/changelog.doc\

tar: all doc
	rm -f meataxe-${MTXVERSION} meataxe-${MTXVERSION}.tar meataxe-${MTXVERSION}.tar.gz \
	&& ln -s . meataxe-${MTXVERSION} \
	&& tar cf meataxe-${MTXVERSION}.tar $(EXPORTED_FILES:%=meataxe-${MTXVERSION}/%) \
	&& rm meataxe-${MTXVERSION} \
	&& gzip meataxe-${MTXVERSION}.tar \
	&& echo "Created meataxe-${MTXVERSION}.tar.gz"
