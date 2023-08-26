# Helper functions and common initialization for tests

set -eu

error()
{
   echo "TEST FAILED: $*"
   exit 1
}

MTX_TESTCASE_DIR="${0%/run}"


[ -n "$MTX_ZZZ" ] || error "\$MTX_ZZZ is not defined"
MTX_TEST_DATA_DIR="../tests/data/$MTX_ZZZ"
[ -d "$MTX_TEST_DATA_DIR" ] \
   || error "Invalid data directory \"$MTX_TEST_DATA_DIR\" (\$MTX_ZZZ must be 0 or 1)"

   
TEST_FIELDS="2 5 9 25 67 125 256"
if [ $MTX_ZZZ -eq 1 ]; then TEST_FIELDS="$TEST_FIELDS 625"; fi


# compareWithReference <file>
# 
# Compares a file created by the test with a reference file from the test case directory.
# If the files are not equal the test fails.
# For example, "compareWithReference xyz" would compare these files:
#    tmp/tmp.xyz
#    tests/TTTT/xyz.expected

compareWithReference()
{
   cmp "tmp.$$.$1" "${MTX_TESTCASE_DIR}/$1.expected" || exit 1
}


# compareBinaryWithReference <file>
# 
# Converts a binary file to text format and compares the result with a reference file.
# The "tmp.$$." prefix is prepended to <file>, and ".text" is appended to the text file.
# For example, "compareBinaryWithReference xyz" is equivalent to 
#
#    zpr tmp/tmp.$$.<file> tmp/tmp.$$.<file>.txt
#    cmp tmp/tmp.$$.<file>.txt tests/TTTT/<file>.txt.expected
#
# where TTTT is the test id.

compareBinaryWithReference()
{
   zpr "tmp.$$.$1" "tmp.$$.$1.txt"
   compareWithReference "$1.txt"
}

