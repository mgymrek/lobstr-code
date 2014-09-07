# make small test reference w/ only chr22

# Get reference set of STRs in tiny region
cat /san/melissa/workspace/lobstr-code-production/resource_bundles/hg19/lobstr_v2.0.3_hg19_ref.bed | awk '($1=="chr22")' | awk '($2 >=38100651 && $2 <= 38120651)' > lobstr_test_ref.bed

# Run the indexer
python /san/melissa/workspace/lobstr-code-production/scripts/lobstr_index.py \
  --str lobstr_test_ref.bed \
  --ref /data/dbase/human/hg19/fasta/hg19.fa \
  --out_dir small_lobstr_ref_v2 -v
