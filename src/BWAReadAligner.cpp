#include "BWAReadAligner.h"
#include <algorithm>
#include "common.h"
#include "runtime_parameters.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* 
TODO!

- cigar strings incorrect for flanking regions with indels
- check strands for complimentary matches
 */

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
  if (align_debug)
    cout << endl << "processing " << read->ID << " " << read->nucleotides << " repseq " << read->repseq << endl;
  const string& repseq = read->repseq;
  if (_bwt_references->find(repseq) == _bwt_references->end()) {
    return false;
  }

  int is_comp = _opts->mode&BWA_MODE_COMPREAD;

  // get bwa_seq_t *seqs
  bwa_seq_t *seqs, *seq_left, *seq_right;
  seqs=(bwa_seq_t*)calloc(2, sizeof(bwa_seq_t)); // TODO memory leak here?
  seq_left = &seqs[0];
  seq_right= &seqs[1];

  // For some reason if I don't do it it doesn't work??
  read->left_flank_nuc = reverseComplement(read->left_flank_nuc);
  read->right_flank_nuc = reverseComplement(read->right_flank_nuc);

  // LEFT FLANKING REGION
  
  const char* left_flank = read->left_flank_nuc.c_str();
  const char* left_qual =  read->quality_scores.substr(0, strlen(left_flank)).c_str();
  seq_left->bc[0] = 0;
  seq_left->tid = -1;
  seq_left->qual = 0;
  seq_left->full_len = seq_left->clip_len = seq_left->len = strlen(left_flank);
  seq_left->seq = (ubyte_t*)calloc(seq_left->len, 1);
  seq_left->qual = (ubyte_t*)calloc(seq_left->len, 1);
  for (int i = 0; i != seq_left->full_len; ++i) {
    seq_left->seq[i] = nst_nt4_table[(int)left_flank[i]];
    if (fastq) {
      seq_left->qual[i] = left_qual[i]+33<126?left_qual[i]+33:126;
    } else {
      seq_left->qual[i] = 126;
    }
  }
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
  seq_right->qual = (ubyte_t*)calloc(seq_right->len, 1);
  for (int i = 0; i != seq_right->full_len; ++i) {
    seq_right->seq[i] = nst_nt4_table[(int)right_flank[i]];
    if (fastq) {
      seq_right->qual[i] = right_qual[i]+33<126?right_qual[i]+33:126;
    } else {
      seq_right->qual[i] = 126;
    }
  }
  seq_right->rseq = (ubyte_t*)calloc(seq_right->full_len,1);
  memcpy(seq_right->rseq, seq_right->seq, seq_right->len);
  seq_reverse(seq_right->len, seq_right->seq, 0);
  seq_reverse(seq_right->len, seq_right->rseq, is_comp);
  seq_right->name = strdup((const char*)read->ID.c_str());

  if (align_debug)
    cout << "left flank " << left_flank << " right flank " << right_flank << endl;

  // call bwa with appropriate options
  bwa_cal_sa_reg_gap(0, _bwt_references->at(repseq).bwt, 2, seqs, _opts); 

  // Find matches between them (need same STR, same direction)
  int num_left_aln = seq_left->n_aln;
  int num_right_aln = seq_right->n_aln;
  if (num_left_aln == 0 || num_right_aln == 0) {
    bwa_free_read_seq(2, seqs);
    return false;
  }
  // get ref seqs and index of alignment for each
  list<ALIGNMENT> left_alignments;
  list<ALIGNMENT> right_alignments;

  // get coordinates
  //cout << "num left " << num_left_aln << " num right " << num_right_aln << endl;

  if (align_debug)
    cout <<"coords left " << endl;
  if (!GetAlignmentCoordinates(seq_left, repseq, &left_alignments)) {
    bwa_free_read_seq(2, seqs);
    return false;
  }
  if (align_debug)
    cout << "coords right " << endl;
  if (!GetAlignmentCoordinates(seq_right, repseq, &right_alignments)) {
    bwa_free_read_seq(2, seqs);
    return false;
  }

  // get shared keys
  ALIGNMENT left_refid;
  ALIGNMENT right_refid;
  if (!GetSharedAln(left_alignments, right_alignments,
  		    &left_refid, &right_refid)) {
    bwa_free_read_seq(2, seqs);
    return false;
  }
  if (align_debug)
    cout << "found shared aln " <<endl;


  // get the desired info out of seqs and put it back into the MSReadRecord
  read->chrom = left_refid.chrom;
  read->msStart = left_refid.start;
  read->msEnd = left_refid.end;
  read->msRepeat = repseq;
  read->refCopyNum = left_refid.copynum;
  read->lDist = seq_left->score;
  read->rDist = seq_right->score;
  read->reverse = !left_refid.strand;
  read->lStart = left_refid.pos;
  read->lEnd = left_refid.pos + strlen(left_flank);
  read->rStart = right_refid.pos;
  read->rEnd = right_refid.pos + strlen(right_flank);
  read->name = left_refid.name;
  int read_len_at_locus = read->nucleotides.size();
  int ref_len_at_locus = abs(read->rEnd - read->lStart);
  if (!left_refid.strand) {
    ref_len_at_locus = abs(read->lEnd - read->rStart);
  }
  read->diffFromRef = read_len_at_locus - ref_len_at_locus;
  bwa_free_read_seq(2, seqs);

  // complement back again
  read->left_flank_nuc = reverseComplement(read->left_flank_nuc);
  read->right_flank_nuc = reverseComplement(read->right_flank_nuc);
  return abs(read->diffFromRef) <= max_diff_ref;
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
      // check strand, L, R are compatible
      /* Need:
	 1. Strand the same
	 2. One right and one left
	 3. If L, R -> repeat is in forward dir. 
	    If R, L -> repeat is in backward dir 
	    (check #3 in calling function) TODO
       */
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
  if (left) {
    refid->start = atoi(items.at(3).c_str())+extend;
  } else {
    refid->start = atoi(items.at(3).c_str());
  }
  if (left) {
    refid->end = atoi(items.at(4).c_str());
  } else {
    refid->end = atoi(items.at(4).c_str())-extend;
  }
  refid->copynum = atof(items.at(6).c_str());
  if (items.size() > 7)
    refid->name = items.at(7);
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
  bwa_aln2seq_core(aligned_seqs->n_aln, aligned_seqs->aln, aligned_seqs,0,1000);

  if (align_debug)
    cout << "nmulti " << aligned_seqs->n_multi << endl;

  // forward
  if (aligned_seqs->strand) {
    bwa_cal_pac_pos_core(_bwt_references->at(repseq).bwt[0], 0, &aligned_seqs[0],
			 _opts->max_diff, _opts->fnr);
  }
  for (int j = 0; j < aligned_seqs[0].n_multi; ++j) {
    bwt_multi1_t *p = aligned_seqs[0].multi + j;
    if (p->strand) p->pos = bwt_sa(_bwt_references->at(repseq).bwt[0], p->pos);
  }

  // reverse
  if (!aligned_seqs[0].strand) {
    bwa_cal_pac_pos_core(0, _bwt_references->at(repseq).bwt[1], &aligned_seqs[0], 
			 _opts->max_diff, _opts->fnr);
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
    if (q->pos >= _bnt_annotations->at(repseq).bns->l_pac) {
      cout << "ERROR" << endl;
      continue;
    }
    nn = bns_coor_pac2real(_bnt_annotations->at(repseq).bns, q->pos, j, &seqid);
    ALIGNMENT refid;
    ParseRefid(_bnt_annotations->at(repseq).bns->anns[seqid].name, &refid);

    if (align_debug) {
      cout << "adding " << _bnt_annotations->at(repseq).bns->anns[seqid].name << endl;
    }    
    refid.strand = q->strand;
    if (q->strand) {
      refid.pos = (refid.start-extend) + (int)(q->pos - _bnt_annotations->at(repseq).bns->anns[seqid].offset) + 1;
    } else {
      refid.pos = (int)(q->pos - _bnt_annotations->at(repseq).bns->anns[seqid].offset - 2*j)
	+ (refid.start-extend);
    }
    alignments->push_back(refid);
  }
  return true;
}



BWAReadAligner::~BWAReadAligner(){}
