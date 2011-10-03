#include "Genotyper.h"
#include "runtime_parameters.h"

using namespace std;

// keep track of genotypes
Genotyper genotyper;

void show_help(){
  const char* help = "lobstr_genotyper [OPTIONS] --aligned <aligned.tab file>\n" \
"-f,--files     file or comma-separated list of files containing reads in fasta or fastq format\n" \
"-t,--table     file containing table of markers to test for. See README for table format\n" \
"-g,--genome    fasta file containing the genome (e.g. hg18.fa)\n" \
"-o,--out       prefix for out put files. will output:\n" \
"                      <prefix>.<filename>.fast(a/q): file of raw reads with detected reads modified for each <filename>\n" \
"                      <prefix>.msalign: file of alignments\n" \
"\n\nOptions:\n" \
"-h,--help                  display this help screen\n" \
"-v,--verbose               print out useful progress messages\n" \
"-q,--fastq                 reads are in fastq format (default: fasta)\n" \
"-p,--threads <int>         number of threads (default: 1)\n" \
"-m,--mismatch <int>        number of mismatches to allow in each flanking region (defult: 0). An alignment is reported if there is a unique best alignment.\n" \
"-s,--sam <filename>        output aligned reads in .sam format\n" \
"--rmdup                    remove PCR duplicates when reporting genotype counts\n" \
"\n\nAdvanced options:\n" \
"--min-read-length    minimum number of nucleotides for a read to be processed. This should be at least two times fft-window-size (default: 48)\n"
"--max-read-length    maximum number of nucleotides for a read to be processed. (default: 200)\n"
"--fft-window-size    size of fft window (default: 24)\n" \
"--fft-window-step    step size of sliding window (default: 12)\n" \
"--lobe-threshold     threshold score to call a window periodic (defualt: 3)\n" \
"--extend             length of flanking regions in the genome to align against (default: 100)\n" \
"--minperiod          minimum period to attempt to detect (default: 1)\n" \
"--maxperiod          maximum period to attempt to detect (default: 8)\n" \
"--minflank           minimum length of flanking region to try to align (default: 10)\n" \
"--maxflank           length to trim flanking regions to if they exceed that length (default: 1000)\n" \
    "--max-diff-ref       maximum difference in length from the reference sequence to allow for alignment (default 50) (will take the absolute value)\n" \
"--pi<pi0, pi1, pi2>  priors for genotype values aa, ab, bb, where a = reference and b = non-reference\n" \
"--mu<mu0, mu1, mu2>  probability of read coming from reference allele for genotypes aa, ab, and bb\n" \
"--female             The sample is from a female. Used in genotyping steps for genotyping STRs from the X and Y chromosomes (default false)\n" \
"--min-coverage       Minimum required coverage of an STR locus to return a genotyper prediction. (default 2)\n" \
"This program takes in raw reads, detects and aligns reads containing microsatellites, and genotypes STR locations.";
	cout << help;
	exit(1);
}

