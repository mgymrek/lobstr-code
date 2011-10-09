#include "BWAReadAligner.h"
#include <algorithm>
#include "common.h"
#include "runtime_parameters.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern unsigned char nst_nt4_table[256];

BWAReadAligner::BWAReadAligner(map<std::string, BWT>* bwt_references,
			       map<std::string, BNT>* bnt_annotations,
			       gap_opt_t *opts) {
  bwase_initialize();
  _bwt_references = bwt_references;
  _bnt_annotations = bnt_annotations;
  _opts = opts;
}

bool BWAReadAligner::ProcessRead(MSReadRecord* read) {
  const string& repseq = read->repseq;
  int is_comp = _opts->mode&BWA_MODE_COMPREAD;

  // get bwa_seq_t *seqs
  bwa_seq_t *seqs, *seq_left, *seq_right;
  seqs=(bwa_seq_t*)calloc(2, sizeof(bwa_seq_t));
  seq_left = &seqs[0];
  seq_right= &seqs[1];

  // LEFT FLANKING REGION
  const char* left_flank = read->left_flank_nuc.c_str();
  const char* left_qual =  read->quality_scores.substr(0, strlen(left_flank)).c_str();
  seq_left->bc[0] = 0;
  seq_left->tid = -1;
  seq_left->qual = 0;
  seq_left->full_len = seq_left->clip_len = seq_left->len = strlen(left_flank);
  seq_left->seq = (ubyte_t*)calloc(seq_left->len, 1);
  for (int i = 0; i != seq_left->full_len; ++i)
    seq_left->seq[i] = nst_nt4_table[(int)left_flank[i]];
  seq_left->qual = (ubyte_t*)strdup(left_qual);
  seq_left->rseq = (ubyte_t*)calloc(seq_left->full_len,1);
  memcpy(seq_left->rseq, seq_left->seq, seq_left->len);
  seq_reverse(seq_left->len, seq_left->seq, 0);
  seq_reverse(seq_left->len, seq_left->rseq, is_comp);
  seq_left->name = strdup((const char*)read->ID.c_str());
				   
  
  // RIGHT FLANKING REGION
  const char* right_flank = read->right_flank_nuc.c_str();
  const char* right_qual = read->quality_scores.substr(read->quality_scores.size() -
						       strlen(right_flank),
						       strlen(right_flank)).c_str();
  seq_right->bc[0] = 0;
  seq_right->tid = -1;
  seq_right->qual = 0;
  seq_right->full_len = seq_right->clip_len = seq_right->len = strlen(right_flank);
  seq_right->seq = (ubyte_t*)calloc(seq_right->len, 1);
  for (int i = 0; i != seq_right->full_len; ++i)
    seq_right->seq[i] = nst_nt4_table[(int)right_flank[i]];
  seq_right->qual = (ubyte_t*)strdup((char*)right_qual);
  seq_right->rseq = (ubyte_t*)calloc(seq_right->full_len,1);
  memcpy(seq_right->rseq, seq_right->seq, seq_right->len);
  seq_reverse(seq_right->len, seq_right->seq, 0);
  seq_reverse(seq_right->len, seq_right->rseq, is_comp);
  seq_right->name = strdup((const char*)read->ID.c_str());

  // call bwa with appropriate options
  bwa_cal_sa_reg_gap(0, _bwt_references->at(repseq).bwt, 2, seqs, _opts); 

  // Find matches between them (need same STR, same direction)
  int num_left_aln = seq_left->n_aln;
  int num_right_aln = seq_right->n_aln;
  if (num_left_aln == 0 || num_right_aln == 0) return false;
  // get ref seqs and index of alignment for each
  list<ALIGNMENT> left_alignments;
  list<ALIGNMENT> right_alignments;

  // get coordinates
  if (!GetAlignmentCoordinates(seq_left, repseq, &left_alignments)) return false;
  if (!GetAlignmentCoordinates(seq_right, repseq, &right_alignments)) return false;
  
  // get shared keys
  ALIGNMENT left_refid;
  ALIGNMENT right_refid;
  if (!GetSharedAln(left_alignments, right_alignments,
  		    &left_refid, &right_refid)) return false;

  // get the desired info out of seqs and put it back into the MSReadRecord
  read->chrom = left_refid.chrom;
  read->msStart = left_refid.start+extend;
  read->msEnd = left_refid.end-extend;
  read->msRepeat = repseq;
  read->refCopyNum = left_refid.copynum;
  read->name = left_refid.name;
  read->lDist = left_refid.score;
  read->rDist = right_refid.score;
  read->reverse = left_refid.strand;
  read->lStart = left_refid.pos;
  read->lEnd = left_refid.pos + strlen(left_flank);
  read->rStart = right_refid.pos;
  read->rEnd = right_refid.pos + strlen(right_flank);
  read->name = left_refid.name;
  int read_len_at_locus = read->nucleotides.size();
  int ref_len_at_locus = abs(read->rEnd - read->lStart);
  if (left_refid.strand) {
    ref_len_at_locus = abs(read->lEnd - read->rStart);
  }
  read->diffFromRef = read_len_at_locus - ref_len_at_locus;
  return true;
}

bool BWAReadAligner::GetSharedAln(const list<ALIGNMENT>& map1,
				  const list<ALIGNMENT>& map2,
				  ALIGNMENT* left_refid,
				  ALIGNMENT* right_refid) {
  // get keys that are shared and consistent
  map<int,ALIGNMENT> left_id_to_ref;
  bool found = false;
  for (list<ALIGNMENT>::const_iterator it = map1.begin();
       it != map1.end(); ++it) {
    left_id_to_ref[(*it).id] = *it;
  }
  for (list<ALIGNMENT>::const_iterator it = map2.begin();
       it != map2.end(); ++it) {
    int ref_key = (*it).id;
    if (left_id_to_ref.find(ref_key) != left_id_to_ref.end()) { // found match!
      // check strand, L, R
      if (left_id_to_ref[ref_key].strand == (*it).strand &&
	left_id_to_ref[ref_key].left != (*it).left) {
	if (found) return false; // multiple
	found = true;
	*left_refid = left_id_to_ref[ref_key];
	*right_refid = *it;
      }
    }
  }
  return found;
 }

void BWAReadAligner::ParseRefid(const string& refstring, ALIGNMENT* refid) {
  vector<string> items;
  split(refstring, '_', items);
  refid->id = atoi(items.at(0).c_str());
  refid->left = (items.at(1) == "L");
  refid->chrom = items.at(2);
  refid->start = atoi(items.at(3).c_str());
  refid->end = atoi(items.at(4).c_str());
  refid->copynum = atof(items.at(5).c_str());
  if (items.size() > 6)
    refid->name = items.at(6);
}

int64_t pos_end_multi(const bwt_multi1_t *p, int len) // analogy to pos_end()
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + len;
}

bool BWAReadAligner::GetAlignmentCoordinates(bwa_seq_t* aligned_seqs,
					     const std::string&repseq,
					     list<ALIGNMENT>* alignments) {
  // fill in alignment properties
  bwa_aln2seq_core(aligned_seqs->n_aln, aligned_seqs->aln, aligned_seqs,0,1);
  // forward
  if (aligned_seqs[0].strand) {
    // TODO change max_mm, fnr, once I know what they do..
    bwa_cal_pac_pos_core(_bwt_references->at(repseq).bwt[0], 0, &aligned_seqs[0], 10, 0.04);
  }
  for (int j = 0; j < aligned_seqs[0].n_multi; ++j) {
    bwt_multi1_t *p = aligned_seqs[0].multi + j;
    if (p->strand) p->pos = bwt_sa(_bwt_references->at(repseq).bwt[0], p->pos);
  }

  // reverse
  if (!aligned_seqs[0].strand) {
    // TODO change max_mm, fnr, once I know what they do..
    bwa_cal_pac_pos_core(0, _bwt_references->at(repseq).bwt[1], &aligned_seqs[0], 10, 0.04);
  }
  for (int j = 0; j < aligned_seqs[0].n_multi; ++j) {
    bwt_multi1_t *p = aligned_seqs[0].multi + j;
    if (!p->strand) p->pos = _bwt_references->at(repseq).bwt[1]->seq_len -
		      bwt_sa(_bwt_references->at(repseq).bwt[1], p->pos) +
		      aligned_seqs[0].len;
  }

  // get coords
  for (int i = 0; i < aligned_seqs->n_multi; ++i) {
    bwt_multi1_t *q = aligned_seqs->multi + i;
    int k, j, nn, seqid; 
    j = pos_end_multi(q, aligned_seqs->len) - q->pos;
    nn = bns_coor_pac2real(_bnt_annotations->at(repseq).bns, q->pos, j, &seqid);
    ALIGNMENT refid;
    ParseRefid(_bnt_annotations->at(repseq).bns->anns[seqid].name, &refid);
    refid.strand = q->strand;
    refid.score = (aligned_seqs->aln +i)->score;
    if (q->strand) {
      refid.pos = refid.start + (int)(q->pos - _bnt_annotations->at(repseq).bns->anns[seqid].offset) + 1;
    } else {
      refid.pos = (int)(q->pos - _bnt_annotations->at(repseq).bns->anns[seqid].offset - 2*j)
	+ refid.start + 1;
    }
    alignments->push_back(refid);
  }
  return true;
}



BWAReadAligner::~BWAReadAligner(){}
