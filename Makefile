include Makefile.conf

MTX_VERSION = 2.4

CFLAGS=$(CFLAGS1) -Itmp
LFLAGS=$(LFLAGS1)

PROGRAMS = \
  cfcomp checksum chop decomp genmod mkcycl mkdotl mkgraph mkhom mkhom_old\
  mkinc mksub mktree orbrep precond pseudochop pwkond rad soc symnew tcond tuc \
  zad zbl zcf zcl zcp zct zcv zef zev zfr ziv zkd zmo zmu zmw znu zor zpo zpr \
  zpt zqt zro zsc zsi zsp zsy ztc zte ztm ztr zts zuk zvp

all build: $(PROGRAMS:%=bin/%)

clean:
	rm -rf bin tmp
	-cd tests; make clean


# ------------------------------------------------------------------------------
# Compile C sources
# ------------------------------------------------------------------------------

tmp/%.o: tmp/mk.dir src/%.c src/meataxe.h tmp/config.h
	@echo "Compiling $*.c"
	@$(CC) $(CFLAGS) -c src/$*.c -o $@

tmp/mk.dir:
	mkdir -p tmp
	touch $@

bin/mk.dir:
	mkdir -p bin
	touch $@

# ------------------------------------------------------------------------------
# Link programs
# ------------------------------------------------------------------------------

bin/%: bin/mk.dir tmp/%.o tmp/libmtx.a
	@echo "Linking $@"
	@$(CC) $(LFLAGS) -o $@ tmp/$*.o tmp/libmtx.a


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

tmp/libmtx.a: $(LIB_OBJS:%=tmp/%.o)
	@echo "Creating $@"
	@rm -f $@
	@ar r $@ $(LIB_OBJS:%=tmp/%.o)


# ------------------------------------------------------------------------------
# Test suite
# ------------------------------------------------------------------------------

TS_OBJS1=c-args c-bitstring c-charpol\
	c-ffio c-fileio c-ffmat c-ffrow c-fpoly \
	c-grease c-kernel c-matins c-matrix c-matset\
	c-os c-perm c-poly c-pseed c-quot c-random \
	c-sets c-stf c-tensor

TS_OBJS=tmp/zzztest.o $(TS_OBJS1:%=tmp/%.o) tmp/libmtx.a

bin/zzztest: bin/mk.dir $(TS_OBJS)
	@echo "Linking $@"
	@$(CC) $(CFLAGS) -o $@ $(TS_OBJS)

bin/checksum: bin/mk.dir tmp/checksum.o
	@echo "Linking $@"
	@$(CC) $(CFLAGS) -o $@ tmp/checksum.o

TESTS=0000\
  0100 0100 0105 0106 0107 0108 0109 0110 0111 0112\
  0113 0114 0115 0116 0117 0118 0200 0201 0210 0211\
  0212 0213 0214 0215

check: tmp/zzztest.done $(TESTS:%=tmp/t-%.done)

tmp/zzztest.done: tmp/mk.dir bin/zzztest
	cd tmp && ../bin/zzztest
	touch $@

tmp/t-%.done: tmp/mk.dir tests/t-% tmp/t.config bin/checksum build
	@echo "t-$* `grep '^#:' tests/t-$* | cut -c 3-100`"
	@cd tmp && ../tests/t-$*
	@touch $@

tmp/t.config: tmp/mk.dir tests/config
	cp tests/config $@

# ------------------------------------------------------------------------------
# config.h
# ------------------------------------------------------------------------------

tmp/config.h: tmp/mk.dir Makefile Makefile.conf src/genconfig.c
	@echo "Generating config.h"
	@$(CC) $(CFLAGS) $(LFLAGS) -o tmp/genconfig src/genconfig.c
	@tmp/genconfig >$@


.PHONY: tar clean install check
.PRECIOUS: tmp/%.o


# ------------------------------------------------------------------------------
# Documentation
# ------------------------------------------------------------------------------

.PHONY: doc

docDir = doc/${MTX_VERSION}
docDocs = src/changelog.doc src/meataxe.doc src/sections.doc
docProducts = ${docDir}/index.html ${docDir}/pages.html ${docDir}/classes.html

doc: ${docProducts}

${docProducts}: \
   etc/Doxyfile etc/layout.xml $(PROGRAMS:%=src/%.c) $(LIB_OBJS:%=src/%.c) src/meataxe.h  \
   ${docDocs} src/meataxe.doc src/changelog.doc
	mkdir -p ${docDir}
	doxygen etc/Doxyfile


# ------------------------------------------------------------------------------
# Releasing
# ------------------------------------------------------------------------------

EXPORTED_FILES =\
  $(PROGRAMS:%=src/%.c)\
  $(LIB_OBJS:%=src/%.c)\
  src/meataxe.h src/genconfig.c\
  src/check.h $(TS_OBJS1:%=src/%.c) $(TS_OBJS1:%=src/%.h) src/zzztest.c\
  $(TESTS:%=tests/t-%) tests/config tests/data\
  Makefile Makefile.conf.dist README COPYING\
  src/meataxe.doc src/changelog.doc\

tar: all rebuild-doc
	V=2.4.`cat svnversion` \
	&& rm -f meataxe-$$V meataxe-$$V.tar meataxe-$$V.tar.gz \
	&& ln -s . meataxe-$$V \
	&& tar cf meataxe-$$V.tar $(EXPORTED_FILES) \
	&& rm meataxe-$$V \
	&& gzip meataxe-$$V.tar \
	&& echo "Created meataxe-$$V.tar.gz"



