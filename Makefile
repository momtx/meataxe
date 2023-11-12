# ------------------------------------------------------------------------------
# Build-time configuration
#
# Don't change the settings here. Instead, write your customized settings to
# Makefile.conf. This will make upgrading to a new version of the MeatAxe
# easier.
# ------------------------------------------------------------------------------

# Directory where binaries and run-time files will be installed.
MTXROOT = ${CURDIR}

# ------------------------------------------------------------------------------
# Common compiler and linker flags
# ------------------------------------------------------------------------------
# The default is to use GNU C
CC=gcc
CFLAGS1=-std=c99 -D_DEFAULT_SOURCE -g -O3 -Wall
LDFLAGS1=-g -Wall
# For debugging:
#CFLAGS1=-std=c99 -D_DEFAULT_SOURCE -g -Wall -Werror -DMTX_DEBUG
#CFLAGS1=-std=c99 -D_DEFAULT_SOURCE -g -Wall -Werror -DASM_MMX -DMTX_DEBUG

# ------------------------------------------------------------------------------
# Multithreading support
# ------------------------------------------------------------------------------
# Leave empty to disable:
CFLAGS_THREADS=
LDFLAGS_THREADS=
# Enable with four threads by default:
#CFLAGS_THREADS=-pthread -DMTX_DEFAULT_THREADS=4
#LDFLAGS_THREADS=-pthread

# ------------------------------------------------------------------------------
# Flags to pass to the ar utility
# ------------------------------------------------------------------------------
ARFLAGS1=
# For AIX 64-bit:
#ARFLAGS1=-X 64

# ------------------------------------------------------------------------------
# Arithmetic kernel
# ------------------------------------------------------------------------------
# Standard kernel, up to GF(256)
ZZZ=0
# Big kernel, up to GF(2^16)
#ZZZ=1

# Verbose output (echo all commands)
V=0

MTXDOCDIR=tmp/doc

# configuration overrides
include Makefile.conf

# ------------------------------------------------------------------------------
# Other settings
# ------------------------------------------------------------------------------

SILENT0=@
SILENT1=
SILENT=${SILENT${V}}

MTXBIN = ${MTXROOT}/bin

CFLAGS=$(CFLAGS1) $(CFLAGS_THREADS) -I"include" -Itmp -DMTX_ZZZ=${ZZZ}
LDFLAGS=$(LDFLAGS1) $(LDFLAGS_THREADS)

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
	touch "$@"
	@sleep 5

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
	hashlittle2 \
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
	mmulscal \
	mprintf \
	mraddgen mrcore mrread mrtranspose mrwrite \
	msclean mscore \
	mtensor mtxobj os \
	permcmp permcore permdup perminv permmul permorder\
	permprint permpwr permread permwrite \
	pex poladd\
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
  0001_zzztest \
  0010_test_data \
  0020_tables \
  0100_zpr_zcv \
  0105_zcf \
  0106_zor \
  0107_ztr \
  0108_zad_zmu \
  0109_zmo \
  0110_zmo_zkd \
  0111_zsi \
  0112_zuk \
  0113_zsy \
  0114_zcp \
  0115_zte \
  0116_zpr_gap \
  0117_zmw \
  0118_zpt \
  0119_zef_zcl \
  0120_zfr \
  0121_zev \
  0122_ztm \
  0123_zsp \
  0200_lattice_m11 \
  0201_lattice_ac \
  0202 \
  0203_lattice_m \
  0210 \
  0211 \
  0212_soc_a5reg \
  0213_endo_m11 \
  0214 \
  0215_m11_x_m11



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

docDocs = src/changelog.md src/mainpage.md
docProducts = ${MTXDOCDIR}/index.html ${MTXDOCDIR}/pages.html ${MTXDOCDIR}/classes.html

doc ${docProducts}: \
   etc/Doxyfile etc/layout.xml $(PROGRAMS:%=src/%.c) $(LIB_OBJS:%=src/%.c) \
   src/meataxe.h  \
   ${docDocs} 
	mkdir -p tmp
	cp etc/Doxyfile tmp/Doxyfile.auto
	echo "PROJECT_NUMBER=${MTXVERSION}" >>tmp/Doxyfile.auto
	echo "OUTPUT_DIRECTORY=${MTXDOCDIR}" >>tmp/Doxyfile.auto
	mkdir -p ${MTXDOCDIR}
	doxygen tmp/Doxyfile.auto >tmp/doxygen.log

doxy:
	cp etc/Doxyfile tmp/Doxyfile.auto
	echo "PROJECT_NUMBER=${MTXVERSION}" >>tmp/Doxyfile.auto
	echo "OUTPUT_DIRECTORY=${MTXDOCDIR}" >>tmp/Doxyfile.auto
	mkdir -p ${MTXDOCDIR}
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
  src/meataxe.md src/changelog.md

tar: all doc
	rm -f meataxe-${MTXVERSION} meataxe-${MTXVERSION}.tar meataxe-${MTXVERSION}.tar.gz \
	&& ln -s . meataxe-${MTXVERSION} \
	&& tar cf meataxe-${MTXVERSION}.tar $(EXPORTED_FILES:%=meataxe-${MTXVERSION}/%) \
	&& rm meataxe-${MTXVERSION} \
	&& gzip meataxe-${MTXVERSION}.tar \
	&& echo "Created meataxe-${MTXVERSION}.tar.gz"
