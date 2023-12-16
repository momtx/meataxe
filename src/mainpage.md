@mainpage Documentation

The MeatAxe is a set of programs for working with matrices over finite fields and,
to some degree, with permutations and polynomials. It was originally written to
calculate modular character tables, but it can be used for other purposes such as
investigating subgroup structure, module structure etc.
Indeed, there is a set of programs to compute the submodule lattice of a given module.
See @ref pg_progs_index for a list of all programs.

To build the MeatAxe, change to the top level directory (containing `Makefile`)
and execute
```
make clean test
```
This command builds all programs in the `bin/` directory and runs the tests.
It should work on Linux, macos and probably other UNIX operating systems
with a C compiler supporting the C99 standard.

Here are some simple examples for using the programs:

* Show the MeatAxe version (works with any program)
  ```
  bin/zpr --version
  ```
* Print a matrix in text format:
  ```
  bin/zpr tests/data/0/C0.1
  ```
* Calculate the order of a matrix:
  ```
  bin/zor tests/data/0/C0.2
  ```
* Calculate the characteristic and minimal polynomial of a matrix:
  ```
  bin/zcp -f tests/data/0/C0.2
  bin/zcp -m -f tests/data/0/C0.2
  ```
* Find the irreducible constituents of a 112 dimensional
  representation of M11 over GF(2):
  ```
  cd tmp
  cp ../tests/data/0/m11.1 ../tests/data/0/m11.2 .
  ../bin/chop m11
  ```

<b>More Information</b>

- @subpage pg_userguide -- Information about installing and using the MeatAxe programs.
- @subpage pg_programming -- Read this if you want to modify or extend the MeatAxe, or
     use the MeatAxe library in your programs.
- @subpage pg_changelog
- @subpage pg_bib


@page pg_building Building and Installing the MeatAxe

# System Requirements

The MeatAxe is developed under Linux and regularly built on Linux and macos.
GNU make and a C compiler supporting the C99 standard is required. 
Building the documentation requires Doxygen.

# Obtaining the Source Code

The MeatAxe source code can be downloaded from https://github.com/momtx/meataxe.
Cloning with Git or unpacking the downloaded archive will create a top-level
directory containing, among other files, @c Makefile and two subdirectories,
@c src and @c tests.

# Building

To build the MeatAxe, change to the top level directory (containing `Makefile`)
and execute
```
make clean test
```
This command builds all programs in the `bin/` directory and runs the tests.


# Compile-Time Configuration (Makefile.conf)

The build process can be customized by creating the file @c Makefile.conf
next to @c Makefile and entering make variable settings in this file.
This is recommended instead of modifying the @c Makefile directly because
@c Makefile.conf will not be overwritten when you download a new version.

The variables that can be overwritten are described at the beginning of @c Makefile.
Here is a (possibly incomplete) list:

* **MTXINSTALLDIR** is the directory where run-time files are installed.
  If you are only working locally in the source directory, make install is not required
  and `${MTXINSTALLDIR}` can be left empty.

  If you want to install the runtime files at a different location, set
  `${MTXINSTALLDIR}` to the name of an exsting directory. Running `make install` will
  copy programs and library files to the `bin`, `lib`, and `include` subdirectory of
  `${MTXINSTALLDIR}`. These subdirectories are created when necessary, but the top-level
  directory must exist.
  
* **CC**, **CFLAGS1**, and **LDFLAGS1** set the compiler command, compile options,
  and link options, respectively.

* **ZZZ** controls which arithmentic module will be used. There are two modules available:
  * @c ZZZ=0 selects the standard arithmetic, supporting fields up to GF(256).
  * @c ZZZ=1 selecte the "large fields" arithmetic module which supports fields up to
    GF(2<sup>16</sup>) but is considerably slower than the standard arithmetic.

* **MTXDOCDIR** sets the directory where the HTML documentation will be created.
  Default is "./doc". If you have define an installation directory (see above) may want
  to build the documentation in a subdirectory of the installation directory. 
  For example: `MTXDCDIR=${MTXINSTALLDIR}/doc`.

* **CFLAGS_THREADS** and **LDFLAGS_THREADS** control support for multithreading.
  Thread support is an experimental feature and disabled by default.
  Turning it on can speed up the pwkond program in some cases.
  To enable threading, the preprocessor symbol @c MTX_DEFAULT_THREADS must be defined
  and should be set to the available number of CPU cores. Depending on the operating
  system you may also need compiler and linker options to enable the POSIX threads
  library. Under Linux, the following settings should work:
  ```
  CFLAGS_THREADS=-pthread -DMTX_DEFAULT_THREADS=4
  LDFLAGS_THREADS=-pthread
  ```

@page pg_userguide User's Guide

- @subpage pg_building
- @subpage pg_using
- @subpage prog_stdopts
- @subpage pg_progs_index
- @subpage pg_file_formats

See also @ref pg_bib.



@page pg_bib Bibliography


@anchor HEBA05
@par [HEBA05]
   D.F. Holt; B. Eick; E.A. O'Brien:
   Handbook of computational group theory, Discrete mathematics and its applications 24 (2005).
   CRC Press, ISBN 978-1-58488-372-2

@anchor HR94
@par [HR94]
    D. F. Holt, S. Rees:
    <em>Testing modules for irreducibility</em>.
    J. Austral. Math. Soc. Ser. A 57 (1994), 1–16.

@anchor LI98
@par [LI98]
    Klaus Lux, Gábor Ivanyos:
    <em>Treating the exceptional cases of the Meataxe</em>.
    February 19, 1998 (unpublished).

@anchor LMR94
@par [LMR94]
    Klaus Lux, Jürgen Müller, and Michael Ringe:
    <em>Peakword condensation and submodule lattices: An application of the Meat-axe</em>.
    J. Symbolic Computation, 17:529-544, 1994.

@anchor Ri98
@par [Ri98]
    Michael Ringe:
    <em>Bemerkungen zur Kondensation von Tensorprodukten irreduzibler Moduln</em>.
    (Nov. 1998, unpublished)

@anchor Pa84
@par [Pa84]
    R. A. Parker:
    <em>The computer calculation of modular characters (the Meat-Axe)</em>.
    In: Computational Group Theory (1984), Academic Press, 267–274.

@anchor Sz98
@par [Sz98]
    Magdolna Szöke:
    <em>Examining Green Correspondents of Weight Modules</em>.
    Aachener Beiträge zur Mathematik, Band 24.
    Wissenschaftsverlag Mainz, Aachen, 1998.

@anchor Wie94 
@par [Wie94]
    Markus Wiegelmann.
    <em>Fixpunktkondensation von Tensorproduktmoduln</em>.
    Diplomarbeit, Lehrstuhl D für Mathematik der RWTH Aachen, 1994.

@anchor BC85
@par [BC85]
    D.J. Benson, J.H. Conway.
    <em>Diagrams for Modular Lattices.</em>
    Journal of Pure and Applied Algebra 37 (1985), 111-116.



@page pg_progs_index Program Index

Thd following list contains a detailed description of each MeatAxe program.
Some command line options, which are common to all programs, are
documented on a separate page (see @ref prog_stdopts).

- @subpage prog_cfcomp
- @subpage prog_chop
- @subpage prog_decomp
- @subpage prog_genmod
- @subpage prog_mkcycl
- @subpage prog_mkdotl
- @subpage prog_mkgraph
- @subpage prog_mkhom
- @subpage prog_mkinc
- @subpage prog_mksub
- @subpage prog_mktree
- @subpage prog_orbrep
- @subpage prog_precond
- @subpage prog_pseudochop
- @subpage prog_pwkond
- @subpage prog_rad
- @subpage prog_soc
- @subpage prog_tcond
- @subpage prog_tuc
- @subpage prog_zad
- @subpage prog_zbl
- @subpage prog_zcf
- @subpage prog_zcl
- @subpage prog_zcp
- @subpage prog_zct
- @subpage prog_zcv
- @subpage prog_zef
- @subpage prog_zev
- @subpage prog_zfr
- @subpage prog_ziv
- @subpage prog_zkd
- @subpage prog_zmo
- @subpage prog_zmu
- @subpage prog_zmw
- @subpage prog_znu
- @subpage prog_zor
- @subpage prog_zpc
- @subpage prog_zpo
- @subpage prog_zpr
- @subpage prog_zpt
- @subpage prog_zqt
- @subpage prog_zro
- @subpage prog_zsc
- @subpage prog_zsi
- @subpage prog_zsp
- @subpage prog_zsy
- @subpage prog_ztc
- @subpage prog_zte
- @subpage prog_ztm
- @subpage prog_ztr
- @subpage prog_zts
- @subpage prog_zuk
- @subpage prog_zvp


# File Conversion Programs
These programs are used to convert between different file formats.

- @subpage prog_zcf
- @subpage prog_zct
- @subpage prog_zcv
- @subpage prog_zpr
- @subpage prog_zpt


# The Lattice Programs {#sec_progs_lattice}

These programs are used to investigate the structure of matrix representations over
finite fields. The first step is always to find the irreducible constituents of the
representation with @ref prog_chop "chop", and to find the corresponding peak words with
@ref prog_pwkond "pwkond".
Here are some examples of what you can do with the programs:
- Calculate the complete submodule lattice of a module (using @ref prog_chop "chop", 
  @ref prog_pwkond "pwkond", @ref prog_mkcycl "mkcycl", @ref prog_mkinc "mkinc",
  @ref prog_mkdotl "mkdotl", @ref prog_mksub "mksub", @ref prog_genmod "genmod").
- Find out if a given irreducible module occurs as a constituent of a module
  (@ref prog_cfcomp "cfcomp").
- Draw the submodule lattice in graphical form (@ref prog_mkgraph "mkgraph").
- Calculate a module's socle and radical series (@ref prog_chop "chop",
  @ref prog_pwkond "pwkond",
  @ref prog_soc "soc", @ref prog_rad "rad").
- Calculate all homomorphisms between two modules, or the endomorphism
  ring of a module (@ref prog_chop "chop", @ref prog_pwkond "pwkond", @ref prog_mkhom "mkhom").
- Decompose a module into direct summands (@ref prog_chop "chop", @ref prog_pwkond "pwkond",
  @ref prog_mkhom "mkhom", @ref prog_decomp "decomp").

# Condensation Programs {#sec_progs_cond}
These programs are used to condense representations. @ref prog_zkd "zkd"
performs a fixed-point condensation of permutation representations.
It can be used after the orbits have been calculated with @ref prog_zmo "zmo".
@ref prog_zuk "zuk" uncondenses vectors.

@ref prog_precond "precond" and @ref prog_tcond "tcond" are used to
condense tensor products or matrix representations, when the restriction to
the condensation subgroup is semisimple.
The algorithm assumes that the irreducible constituents
of the restriction, and corresponding peak words are known, so you must
run @ref prog_chop "chop" and @ref prog_pwkond "pwkond" before.


@page pg_using Using the MeatAxe Programs: General Remarks

Each of the programs is self-contained, reading its input from files, and writing
its output to files. 

# Setting up the environment

It is recommended that you include the MeatAxe installation directory, where
all the programs are installed, in your path. You may also define the following
environment variables which are recognized by the MeatAxe programs:

Variable  | Description
:-------- | :-------------
MTXLIB    | Name of the MeatAxe library directory, see @ref sec_libdir.


# Arithmetic Tables
All programs use lookup tables for row operations and finite field arithmetic.
These tables are read from a file named `pXXX.zzz`, or, if the MeatAxe
was build with the "large fields" arithmetic kernel, from `pXXXXX.zzz`.
For example, tables for GF(25) are read from the file `p025.zzz` or
`p00025.zzz`, respectively.

The table file can be in the library directory (see below) or in the current
working directory. If it is not found in either location, a new file is created
in the library directory.

# File Types
All MeatAxe programs operate on files and do not require user interaction while running.
There are two kinds of files:

- Binary (internal format) data files. Programs usually read
  and write data in binary format.

- Text files, used by some programs as @ref prog_zcv "zcv" and @ref prog_zpr "zpr".
  The format of text files is described with the @ref prog_zcv "zcv" program.

Typically, a program reads a number of input files and produces one or
more output files. For example, the @ref prog_zmu "zmu" (multiply) program expects two
input and one output file. Thus,
```  
zmu mat1 mat2 result
```
reads matrices from "mat1" and "mat2", and writes the product to "result".
To find out which
files are used by a specific program, look up the program description in
section @ref pg_progs_index. You can also run any MeatAxe program
with the "--help" option to get an on-line help.

A detailed description of the MeatAxe file formats can be found in @ref pg_file_formats.

# The Library Directory  {#sec_libdir}
When given a file name, all programs look for the file in the current directory first.
If the file is not found there, some programs try to find it in the library directory.
This extended search applies, for example, to @ref prog_zev "zev" and
@ref prog_zcv "zcv" input files, and to the arithmetic table files.

By default, the library is derived from the executable name by replacing "bin" with "lib".
For the program @c /usr/local/mtx/bin/zad will by default expect library files in
@c /usr/local/mtx/lib.
The default may be overridden at run-time by defining the environment variable @c MTXLIB.
Both the default and @c MTXLIB can be overridden with -L command line option, which is
supported by all programs.


@page pg_programming Programmer's Guide
This section is intended for users who want to add new programs to the MeatAxe or
use MeatAxe functions in their own programs. 

Internally, the MeatAxe is built in several layers where each layer uses
services provided by the lower layers. Layers 1 to 3 constitute the MeatAxe Library.

# Layer 1 - Kernel
The MeatAxe kernel provides the finite field arithmetic,
some low-level vector operations,
and a platform-independent interface to operating
system services, including the interaction with the user.
All kernel functions and global variables have names beginning with "Ff"
(finite field), "Mtx" (general functions), or "Sys" (OS interface).
Kernel functions generally do not allocate memory dynamically and do not
completetely check their arguments.
Most kernel services are simple,
but there are some more complex functions like Gaussian elimination,
which have been implemented in the kernel for performance reasons.

- @ref os
- @ref ff
- @subpage pg_error_handling

# Layer 2 - Objects
This layer encapsulates low-level functionality in "objects" 
like matrices, permutations, and polynomials. It also provides the
elementary operations with these objects like, for example, matrix
multiplication.
- @ref mat
- @ref poly
- @ref perm
- @ref bs
- @ref mf
- @ref mrep
- @ref imat
- @ref matset

# Layer 3 - Algorithms
This layer provides more complex operations which typically involve
several objects. The spin-up and split functions are examples for this
layer.
- @ref charpol
- @ref spinup
- @ref wgen

Throughout this documentation, layers 1 to 3 are referred to as the
"MeatAxe library". Indeed, these layers are implemented as a static library.

# Layer 4 - Programs
The top layer includes the MeatAxe programs and shell scripts like
@ref prog_zad "zad", @ref prog_zmu "zmu", @ref prog_chop "chop", ...
See @ref pg_progs_index for detailed descriptions.





# Compiling and Linking with the MeatAxe Library

To use the MeatAxe library in your programs you need two files:

@par meataxe.h
(in the @c src directory) contains declarations for all MeatAxe library functions.
This header must be included in all source files using MeatAxe functions.

@par libmtx.a
The (static) MeatAxe library. This file can be found in the @c tmp directory after
the MeatAxe was built. Add the library when linking your programs. For example:

Example:
```
cc -c -I/usr/joe/mtx/src -o myprog.o myprog.c
cc -L/usr/joe/mtx/tmp -o myprog myprog.o libmtx.a
```

@page prog_stdopts Standard Command Line Options
All MeatAxe programs expect one or more arguments, ususally the names
of input and output files. Arguments may be preceeded by program options
which always start with a minus sign. Options which don't need arguments
can be grouped together. For example, the two commands
```
chop -G -V -V -V test
chop -GVVV test
```
are equivalent. Some options require an additional argument which must be
separated from the option by a space. Here is an example:
```
chop -g 3 -T 3600 -G -VV module
```
Some options may have an alternative, long name starting with `--'.
These options cannot be grouped together with other options.

@section co Universal Options
The following command line options are understood by all MeatAxe programs:

@par \--help
Print a help message and exit.

@par -Q
Quiet, fewer messages. A single -Q is equivalent to --log=:warning (warnings
and error messages), -QQ is equivalent to --log=stdout:error (only error messages)

@par -V
Verbose, more messages. A single -V turns on debug messages (equivalent to 
--log=:debug), -VV enables verbose debug messages (equivalent to --log=:debug2).

@par \--log=[FILE]:LEVEL[:FORMAT]
Configures message output. The argument consists of up to three parts separated
by colons.<br>
<i>FILE</i> is is the name of the log file or 'stdout' to write messages
to the standard output. The special syntax +<i>FILE</i> means that messages
are appended to the file instead of overwriting its contents.<br>
<i>LEVEL</i> is one of the following keywords: 'error', 'warning', 'info', 'debug', 'debug2'.<br>
<i>FORMAT</i> determines the output format and must be one of the following keywords.
'none' (message text only), 'short' (level and thread ID), 'full' (level, thread ID,
and timestamp). The thread ID is only printed if threads are enabled (see the explanation
of CFLAGS_THREADS in @ref pg_building). The default is "none".<br><br>
<i>FILE</i> and <i>FORMAT</i> can be omitted and default to "stdout" and "none",
respectively. Here are some examples:
```
--log=meataxe.log:
--log=:debug
--log=:debug2:full
```

@par -T @em Time
Limit CPU time to \e Time seconds.

@par -L @em LibDir
Used to specify a different library directory. This overrides both
the directory name compiled into the programs and the value of
MTXLIB.

@par \--threads[=N]
Enables multithreading. This option is currently implemented for
@ref prog_pwkond "pwkond" only. Other program accept but ignore `--threads`.
The value defines the number of worker threads and should normally be equal
to the number of available CPU cores. If not specified, N takes its default
value from the build-time option MTX_DEFAULT_THREADS.
If the MeatAxe was build without thread support,
this option is always ignored, even for programs that support multithreading.

@par -o @em Options
This option is used to set various options for testing and debugging
the MeatAxe. Even if documented in the program descriptions, -o
options may be changed or removed in the next release, so you should
not rely on them.

@section oo Other Commonly Used Options
The following options are not supported by all programs, but their
meaning is always the same.

@par -g @em NGen
  Set the number of generators. This option must be followed
  by a positive integer.
@par -G
  Generate output in GAP format.


@page pg_file_formats File Formats

@section sec_fileformats_binary Binary Data Files

Binary data files contain a sequence of objects.
Each object consists of a 12 byte header, followed by the data part.
The object header consists of three 32-bit integers, T, R, and C,
which are are stored in little-endian format.
The first number, T, specifies which kind of object follows.
R and C carry additional information about the object.
Their exact meaning depends on the object type:
@verbatim
T   | Type of object   | R             |  C 
----|------------------|---------------|-----------------------
> 1 | Matrix over GF(T)| no. of rows   | no. of columns
-1  | Permutation      | no. of points | no. of permutations
-2  | Polynomial       | field         | degree
-3  | Bit string       | size          | 0
-8  | Integer matrix   | no. of rows   | no. of columns
@endverbatim

Note: most MeatAxe programs expect their input objects in separate files
rather than reading multiple objects from a single file. This is true,
for example, for the set of generators that is used by programs like
@ref prog_zsp or @ref prog_chop.
For this reason, most data files contain only a single object.


@subsection sec_fileformats_binary_mat Matrices
Matrices are stored as a sequence of rows in the internal format
(see @ref ff), but the rows are truncated to their actual length. 
For example, a row with 17 entries over GF(125) on a machine with
sizeof(long)=4 would occupy 20 bytes in memory but only 17 bytes in
a data file. Depending on the arithmetic kernel used (see @ref ff)
the may also be unused bits in each byte. These unused bits are always
set to zero.

@subsection sec_fileformats_binary_perm Permutations
A permutation of degree N is stored as N 32-bit integers in little-endian
format, the i-th number being the image of i (for i = 0...N-1).

For example, the permutation (1 2)(3 4 5) on five points is represented
as (0 1)(2 3 4). The following hexdump shows the corresponding data file:
@verbatim
Offset  Content
0000000 ff ff ff ff     Header T = -1 (Permutation)
0000004 05 00 00 00     Header R = 5 (degree)
0000008 01 00 00 00     Header C = 1 (number of permutations)
000000c 01 00 00 00     p(0) = 1
0000010 00 00 00 00     p(1) = 0
0000014 03 00 00 00     p(2) = 3
0000018 04 00 00 00     p(3) = 4
000001c 02 00 00 00     p(4) = 2
@endverbatim

For compatibility with older MeatAxe versions, the program can also read
a different format where the points are numbered from 1 to N.

@subsection sec_fileformats_binary_poly Polynomials
A polynomial of degree n,
p(x)=a<sub>n</sub>x<sup>n</sup>+...+a<sub>1</sub>x+a<sub>0</sub>,
over GF(q) is stored with the header (-2,n,q). 
The header is followd by one row of size n+1 containing the coefficients
in ascendig order: (a<sub>0</sub>,a<sub>0</sub>,...a<sub>n</sub>).

@subsection sec_fileformats_binary_imat Integer matrices
An r-by-c matrix is stored with the header (-8,r,c), followed by the marks
in theit "natural" order, i.e., starting with a<sub>11</sub> and ending
with a<sub>rc<</sub>. The marks are stored in 32-bit little-endian format.

@subsection sec_fileformats_binary_bs Bit strings
A bit string of length n is stored with the header (-3,n,0), followed by
k 32-bit little-endian numbers, where k is chosen such that
32k ≤ n < 32k+1. The least significant bit of the first number is the
first bit of the string.





@section sec_fileformats_text Text Files
@todo old format
@todo new format


@page pg_error_handling Error Handling

Error Handling
--------------

Most MeatAxe functions perform additional checks to detect errors caused
by invalid arguments, corrupted files or wrong assumptions made by the
caller (e.g. trying to load a permutation from a file which contains a matrix).
If an error occurs, a message is written to stderr and the program is terminated.
A custom error handler can be used to implement an alternative error reporting
method. However, the process will always be terminated.

Raising Errors
--------------

@ref mtxAbort
@c MTX_HERE
@ref MtxSourceLocation


Providing Context Information
-----------------------------

@ref mtxBegin
@ref mtxEnd


Custom Error Handlers
---------------------


**/
