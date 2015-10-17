# The MeatAxe

The MeatAxe is a set of programs for working with matrices over finite fields.
Its primary purpose is the calculation of modular character tables, although
it can be used for other purposes, such as investigating subgroup structure,
module structure etc.


## Building

To build with the default settings, just type

    make clean test

This will build the programs in `./bin` and run the tests.

Some aspects of the build can be configured by setting variables. See Makefile
for details. You should write your custom settings to `Makefile.conf` rather than
modifying `Makefile`. This prevents your settings from being overwritten when you
unpack a new version of the MeatAxe.


## Documentation

The current documentation for version 2.4 is available under http://momtx.github.io/meataxe
