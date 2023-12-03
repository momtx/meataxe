@page pg_changelog Version History and Release Notes

Most of the programs (including their documentation) in this package
are based on the FORTRAN code written by Richard A. Parker.
The translation from FORTRAN to C was done in 1989 and since then both
versions have developed independently.


# Release 2.5.X {#r2_5_x}

- **NOTE**: version 2.5. includes major incompatible API changes
  and will probably be renamed to 3.0
- Removed suport for old text file formats.
- Added "large field" arithmetic (q < 65535).
- Added support for multithreading (experimental)
- Program messages can be controlled with --log. -Q and -V still exist but
  should not be used together with --log
- GAP output is not written to a file by default to avoid mixing with
  normal messages. See the documentaion of -G / --gap


# Release 2.4.X {#r2_4_x}
@par Jan-2015
Refactoring of tests (work in progress).

@par Jan-2015
The documentation was updated ad extended (work in progress).

@par Aug-2012
The documentation was updated ad extended (but is still far from complete).

@par 21-May-2012
Bugfix: pwkond can now print words in GAP format if there are more than 9 generators.

@par 31-Jul-2011
- Improved GAP support in PWKOND. Output now includes peak words and associated polynomials.
  Note that the output format has changed and is not compatible with previous releases.
- Added functions for handling of dynamic strings.
- Documentation enhancements.

# Release 2.4.13 {#r2_4_13}
@par 11-May-2009 
- Makefile reorganization and SVN integration. As a consequence of the latter, the third
  component of the version number is now the subversion revision.
@par 23-Apr-2009 
- genmod: fixed broken bit string reading code.
- pseudochop: integrated changes by Klaus Lux to get it working under 2.4.8.
  Now using simpler integer mat functions.

# Release 2.4.8 {#r2_4_8}
@par 21-Apr-2009 
- Minor corrections.
- ZSC: accept spin-up scripts starting with a seed vector differnent from 1.

# Release 2.4.7 {#r2_4_7}
@par 09-Sep-2007
- More bugfixes for 64-bit platforms.
@par 03-Sep-2007
- New: pseudochop
@par 02-Sep-2007
- Various corrections for 64-bit platforms
@par 15-Jul-2007
- Initialize all Makefile variables
@par 29-Apr-2007
- Fixed a bug in ZNU, leading to wrong output if the input is not a square matrix and
  the null-space is not requested (i.e., when using the "znu matrix" syntax,
  "znu matrix nsp" was working correctly).
@par 19-Feb-2007
- Fixed various compiler errors and warnings with newer GCC versions.
@par 12-Jul-2005
- Bugfix in zsp (option handling).  Working on documention.

@par 01-Jun-2004
- Bugfix in @c FfAlloc() for large allocations > 4GB.

@par 14-Apr-2004
- Documentation

@par 10-Apr-2004
- Added -n option to ZTS
- Started to convert documentation do doxygen.

@par ?
- Various changes for 64-bit systems
- Bit String file format is now machine independent

@par 09-Dec-2003
- Renamed restruct to @c ffrestrict to avoid name collision on HP machines
- mkhom: Updated online help.
- mkhom updated test scripts
- Test scripts modified to work with Bourne shell (Solaris)
- Automatic platform detection (config.h)

