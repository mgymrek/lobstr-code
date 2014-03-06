#!/bin/sh

##
## Test harness for lobSTR
## Prepares few environment variables and verify reference/binary locations.
##

## Default values
SCRIPT=$(perl -MCwd -e '$f=Cwd::realpath($ARGV[0]); die unless -e $f; print $f' "$0") || exit 1
TESTDIR=$(dirname "$SCRIPT") || exit 1
BASEDIR="$TESTDIR/../"
LOBSTR_BIN="$BASEDIR/src/lobSTR"
LOBSTR_INDEX_PREFIX="$TESTDIR/smallref/small_lobstr_ref/lobSTR_"
INPUT_FILE="$TESTDIR/tiny.fq"
USE_GDB=
USE_VALGRIND=
LOBSTR_THREADS=
SHOW_HELP=

## Read options from the command line
set -- $(getopt hdvp:i:f: "$@")
while [ $# -gt 0 ]
do
    case "$1" in
    (-h) SHOW_HELP=yes;;
    (-d) USE_GDB=yes;;
    (-v) USE_VALGRIND=yes;;
    (-p) LOBSTR_THREADS="-p $2"; shift;;
    (-i) LOBSTR_INDEX_PREFIX="$2"; shift;;
    (-f) INPUT_FILE="$2"; shift;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*)  break;;
    esac
    shift
done

if [ "x$SHOW_HELP" = "xyes" ]; then
BASE=$(basename "$0")
cat<<EOF
lobSTR test script.

Usage:
  $BASE [OPTIONS] -- [LOBSTR-OPTIONS]

Options:
 -d        =  Run under GDB
 -v        =  Run under valgrind
 -p N      =  Use N threads (default = single-threaded)
 -i PREFIX =  Use index PREFIX (default = '$LOBSTR_INDEX_PREFIX')
 -f FILE   =  Use input FILE (default = '$INPUT_FILE').
              By default, file must be FASTQ. see below for examples below
              for other types.

To pass extra parameters to directly to lobSTR, specify them AFTER a double-dash.

Examples:

  # Default, use provided input file and the provided (very small) reference
  \$ $BASE

  # Same as above, but run with 5 threads
  \$ $BASE -p 5

  # Same as above, run under GDB
  # (if lobSTR crashes, GDB will interject and allow some debugging)
  \$ $BASE -d -p 5

  # Use a custom input file
  \$ $BASE -f /data/test.fq

  # Use a custom lobSTR reference index
  \$ $BASE -i /data/lobSTR/hg19/lobSTR_

  # Use a BAM file as input
  # (all options after the double-dash are passed directly to lobSTR)
  \$ $BASE -f /data/test.bam -- --bam

EOF
    exit 0
fi


## Varify parameters (after command-line processing),
## Set some more ENV values
[ -d "$TESTDIR" ] || { echo "Error: can't detect test-files directory" >&2 ; exit 1 ; }
[ -d "$BASEDIR" ] || { echo "Error: lobSTR base directory ($BASEDIR) doesn't exist." >&2 ; exit 1 ; }
[ -x "$LOBSTR_BIN" ] || { echo "Error: can't find lobSTR executable ($LOBSTR_BIN)" >&2 ; exit 1 ; }
LOBSTR_INDEX_DIR=$(dirname "$LOBSTR_INDEX_PREFIX")
[ -d "$LOBSTR_INDEX_DIR" ] || { echo "Error: can't find lobSTR reference directory ($LOBSTR_INDEX_DIR)" >&2 ; exit 1 ; }
LOBSTR_INDEX_REF_FILE="${LOBSTR_INDEX_PREFIX}ref.fa"
[ -e "$LOBSTR_INDEX_REF_FILE" ] || { echo "Error: Invalid lobSTR reference index (expecting file '$LOBSTR_INDEX_REF_FILE' based on index prefix '$LOBSTR_INDEX_PREFIX'" >&2 ; exit 1 ; }

[ "x$USE_GDB" = "xyes" ] && GDB_COMMAND="gdb --eval-command=run --eval-command=quit --args"
[ "x$USE_VALGRIND" = "xyes" ] && VALGRIND_COMMAND="valgrind --leak-check=yes"


##
## Setup output directory, and stdout/stderr logging to a file.
##
OUTDIR=$(mktemp -d) || exit 1
STDOUTPIPE=$OUTDIR/stdout.pipe
STDERRPIPE=$OUTDIR/stderr.pipe
mkfifo "$STDOUTPIPE" || exit 1
mkfifo "$STDERRPIPE" || exit 1
trap 'rm "$STDOUTPIPE" "$STDERRPIPE"' EXIT
tee $OUTDIR/stdout.log < "$STDOUTPIPE" &
tee $OUTDIR/stderr.log < "$STDERRPIPE" >&2 &


$VALGRIND_COMMAND \
  $GDB_COMMAND \
    $LOBSTR_BIN --index-prefix "$LOBSTR_INDEX_PREFIX" \
	-f "$INPUT_FILE" -q \
	$LOBSTR_THREADS \
	--out $OUTDIR/test \
	--rg-sample test \
	--rg-lib test \
	-v "$@" >"$STDOUTPIPE" 2>"$STDERRPIPE"

BAM=$OUTDIR/test.aligned.bam
[ -e "$BAM" ] || { echo "test failed: BAM file ($BAM) was not created" >&2 ;  }
[ -s "$BAM" ] || { echo "test failed: BAM file ($BAM) is empty" >&2 ; }
samtools view "$BAM" > /dev/null || { echo "test failed: samtools detected invalid BAM file ($BAM)">&2 ; }


echo ""
echo "Test complete"
echo "   INPUT_FILE: $INPUT_FILE"
echo "   LOBSTR_INDEX_PREFIX: $LOBSTR_INDEX_PREFIX"
if [ -z "$LOBSTR_THREADS" ]; then
  echo "   THREADS: Single-Threaded"
else
  echo "   THREADS: $LOBSTR_THREADS"
fi
if [ ! -z "$@" ]; then
echo "     Extra parameters to lobSTR: $@"
fi
echo " Output directory: $OUTDIR"
echo "   STDERR: $OUTDIR/stderr.log"
echo "   STDOUT: $OUTDIR/stdout.log"
echo ""
echo " Files in output directory:"
find "$OUTDIR/" -type f -printf "   %p\n"
echo ""
