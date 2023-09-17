# Helper functions and common initialization for tests

set -eu

error()
{
   echo "TEST FAILED: $*"
   exit 1
}

MTX_TESTCASE_DIR="${0%/run}"
MTX_TESTCASE_NO="${MTX_TESTCASE_DIR##*/}"
MTX_TESTCASE_NO="${MTX_TESTCASE_NO%-*}"

[ -n "$MTX_ZZZ" ] || error "\$MTX_ZZZ is not defined"

# Create test directory
rm -rf "T$MTX_TESTCASE_NO"
mkdir "T$MTX_TESTCASE_NO"
cd "T$MTX_TESTCASE_NO"
PATH="../../bin:$PATH"
MTXLIB="../../lib"
export MTXLIB
MTX_TEST_DATA_DIR="../../tests/data/$MTX_ZZZ"
MTX_TESTCASE_DIR="../$MTX_TESTCASE_DIR"
[ -d "$MTX_TEST_DATA_DIR" ] \
   || error "Invalid data directory \"$MTX_TEST_DATA_DIR\" (\$MTX_ZZZ must be 0 or 1)"

TEST_FIELDS="2 5 9 25 67 125 256"
if [ $MTX_ZZZ -eq 1 ]; then TEST_FIELDS="$TEST_FIELDS 625"; fi

# Prefix for temporary files
TF="$MTX_TESTCASE_NO"

# compareWithReference <file> [<reference>]
# 
# Compares a file created by the test with a reference file from the test case directory.
# If the files are not equal the test fails.
# <reference> defaults to <file>.expected.

compareWithReference()
{
   if [ $# -lt 2 ]; then
      cmp "$1" "${MTX_TESTCASE_DIR}/$1.expected" || exit 1
   else
      cmp "$1" "${MTX_TESTCASE_DIR}/$2" || exit 1
   fi
}

# compareBinaryWithReference <file>
# 
# Converts a binary file to text format and compares the result with a reference file.

compareBinaryWithReference()
{
   zpr "$1" "$1.txt"
   compareWithReference "$1.txt"
}

