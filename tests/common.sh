# helper functions for tests

TEST_FIELDS="2 5 9 25 67 125 256"

error()
{
   echo "TEST FAILED: $*"
   exit 1
}
