#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cassert>
#include <ctype.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include "ds.h"

using namespace std;

/**
 * Implementations of the various passes that Qsim makes over SAM files.
 *
 * NOTE: uses strtok so not multithreaded.
 */

/* 64K buffer for all input and output */
const static size_t BUFSZ = 65536;

struct OpRunOffset {
	char op;
	int run;
	int offset;
	
	void init(char _op, int _run, int _offset) {
		op = _op;
		run = _run;
		offset = _offset;
	}
};

struct Alignment {
	
	Alignment() { clear(); }
	
	~Alignment() { }
	
	void clear() {
		rest_of_line = NULL;
		valid = false;
		qname = NULL;
		flag = 0;
		rname = NULL;
		pos = 0;
		mapq = 0;
		cigar = NULL;
		rnext = NULL;
		pnext = 0;
		tlen = 0;
		seq = NULL;
		len = 0;
		qual = NULL;
		mdz = NULL;
		cigar_equal_x = false;
		best_score = 0;
		left_clip = 0;
		right_clip = 0;
		rf_aln_buf.clear();
		rd_aln_buf.clear();
		edit_xscript.clear();
		cigar_ops.clear();
		cigar_run.clear();
		mdz_char.clear();
		mdz_oro.clear();
	}
	
	inline bool is_aligned() const {
		return (flag & 4) == 0;
	}
	
	inline bool is_fw() const {
		return (flag & 16) == 0;
	}
	
	inline bool is_concordant() const {
		return (flag & 2) != 0;
	}
	
	inline bool is_paired() const {
		return (flag & 1) != 0;
	}
	
	inline char mate_flag() const {
		return ((flag & 64) != 0) ? '1' : (((flag & 128) != 0) ? '2' : '0');
	}
	
	/**
	 * Return value is first pre-comma token of ZT:Z strings.s
	 */
	char * parse_extra(char *extra) {
		char *ztz = NULL;
		extra = strtok(extra, "\t");
		bool found_ztz = false, found_mdz = false;
		while(extra != NULL && (!found_mdz || !found_ztz)) {
			if(strncmp(extra, "ZT:Z:", 5) == 0) {
				ztz = extra + 5;
				found_ztz = true;
			}
			if(strncmp(extra, "MD:Z:", 5) == 0) {
				assert(mdz == NULL);
				mdz = extra + 5;
				mdz_to_list();
				found_mdz = true;
			}
			extra = strtok(NULL, "\t");
		}
		if(cigar != NULL && mdz != NULL && !cigar_equal_x) {
			cigar_and_mdz_to_edit_xscript();
		}
		assert(ztz != NULL);
		char *ztz_tok = strtok(ztz, ",");
		assert(ztz_tok != NULL);
		return ztz_tok;
	}
	
	/**
	 * CIGAR string to list.
	 */
	void parse_cigar() {
		assert(cigar_ops.empty());
		assert(cigar_run.empty());
		assert(cigar != NULL);
		const size_t clen = strlen(cigar);
		int i = 0;
		while(i < clen) {
			assert(isdigit(cigar[i]));
			int run = 0;
			do {
				run *= 10;
				run += ((int)cigar[i] - (int)'0');
				i++;
			} while(isdigit(cigar[i]));
			assert(isalpha(cigar[i]));
			if(cigar_ops.empty() && cigar[i] == 'S') {
				left_clip = run;
			} else if(i+1 >= clen && cigar[i] == 'S') {
				right_clip = run;
			}
			if(cigar[i] == 'X' || cigar[i] == '=') {
				cigar_equal_x = true;
			}
			cigar_ops.push_back(cigar[i]);
			cigar_run.push_back(run);
			i++;
		}
		if(cigar_equal_x) {
			cigar_to_edit_xscript();
		}
	}
	
	/**
	 * MD:Z to list.
	 */
	void mdz_to_list() {
		assert(mdz_char.empty());
		assert(mdz_oro.empty());
		assert(mdz != NULL);
		const size_t mlen = strlen(mdz);
		int i = 0;
		while(i < mlen) {
			if(isdigit(mdz[i])) {
				// Matching stretch
				int run = 0;
				do {
					run *= 10;
					run += ((int)mdz[i] - (int)'0');
					i++;
				} while(i < mlen && isdigit(mdz[i]));
				if(run > 0) {
					mdz_oro.expand();
					mdz_oro.back().init(0, run, -1);
				}
			} else if(isalpha(mdz[i])) {
				// Mismatching stretch
				int run = 0;
				do {
					mdz_char.push_back(mdz[i++]);
					run++;
				} while(i < mlen && isalpha(mdz[i]));
				mdz_oro.expand();
				mdz_oro.back().init(1, run, (int)(mdz_char.size() - run));
			} else if(mdz[i] == '^') {
				i++; // skip the ^
				int run = 0;
				while(i < mlen && isalpha(mdz[i])) {
					mdz_char.push_back(mdz[i++]);
					run++;
				}
				mdz_oro.expand();
				mdz_oro.back().init(2, run, (int)(mdz_char.size() - run));
			} else {
				fprintf(stderr, "Unexpected character at position %d of MD:Z string '%s'\n", i, mdz);
			}
		}
		assert(i == mlen);
	}
	
	void cigar_to_edit_xscript() {
		assert(cigar_equal_x);
		assert(edit_xscript.empty());
		for(size_t i = 0; i < cigar_run.size(); i++) {
			char cop = cigar_ops[i];
			int crun = cigar_run[i];
			assert(cop != 'M');
			if(cop == 'S') {
				// TODO: insert something else here
			}
			for(int j = 0; j < crun; j++) {
				edit_xscript.push_back(cop);
			}
		}
		edit_xscript.push_back(0);
	}

	/**
	 * Turn the CIGAR and MD:Z fields into an edit transcript.
	 *
	 * MODIFIES mdz_run
	 */
	void cigar_and_mdz_to_edit_xscript() {
		assert(!cigar_equal_x);
		assert(edit_xscript.empty());
		size_t mdo = 0;
		size_t rdoff = 0;
		for(size_t i = 0; i < cigar_run.size(); i++) {
			char cop = cigar_ops[i];
			int crun = cigar_run[i];
			assert(cop != 'X' && cop != '=');
			if(cop == 'M') {
				int mdrun = 0;
				int runleft = crun;
				while(runleft > 0 && mdo < mdz_oro.size()) {
					char op_m = mdz_oro[mdo].op;
					int run_m = mdz_oro[mdo].run;
					int run_comb = min<int>(runleft, run_m);
					runleft -= run_comb;
					assert(op_m == 0 || op_m == 1);
					if(op_m == 0) {
						for(size_t j = 0; j < run_comb; j++) {
							edit_xscript.push_back('=');
						}
					} else {
						assert(run_m == run_comb);
						for(size_t j = 0; j < run_m; j++) {
							edit_xscript.push_back('X');
						}
					}
					mdrun += run_comb;
					rdoff += run_comb;
					if(run_comb < run_m) {
						assert(op_m == 0);
						mdz_oro[mdo].run -= run_comb;
					} else {
						mdo++;
					}
				}
			} else if(cop == 'I') {
				for(size_t j = 0; j < crun; j++) {
					edit_xscript.push_back('I');
				}
				rdoff += crun;
			} else if(cop == 'D') {
				char op_m = mdz_oro[mdo].op;
				int run_m = mdz_oro[mdo].run;
				assert(op_m == 2);
				assert(crun == run_m);
				assert(run_m == crun);
				mdo++;
				for(size_t j = 0; j < run_m; j++) {
					edit_xscript.push_back('D');
				}
			} else if(cop == 'N') {
				for(size_t j = 0; j < crun; j++) {
					edit_xscript.push_back('N');
				}
			} else if(cop == 'S') {
				// TODO: insert something else hwere
				rdoff += crun;
			} else if(cop == 'H') {
				// pass
			} else if(cop == 'P') {
				throw 1;
			} else if(cop == '=') {
				throw 1;
			} else if(cop == 'X') {
				throw 1;
			} else {
				throw 1;
			}
		}
		assert(mdo == mdz_oro.size());
		edit_xscript.push_back(0);
	}

	/**
	 * MODIFIES mdz_run
	 */
	void cigar_and_mdz_to_stacked() {
		size_t mdo = 0;
		size_t rdoff = 0;
		for(size_t i = 0; i < cigar_run.size(); i++) {
			char cop = cigar_ops[i];
			int crun = cigar_run[i];
			if(cop == 'M') {
				int mdrun = 0;
				int runleft = crun;
				while(runleft > 0 && mdo < mdz_oro.size()) {
					char op_m = mdz_oro[mdo].op;
					int run_m = mdz_oro[mdo].run;
					int st_m = mdz_oro[mdo].offset;
					int run_comb = min<int>(runleft, run_m);
					runleft -= run_comb;
					assert(op_m == 0 || op_m == 1);
					for(size_t j = rdoff; j < rdoff + run_comb; j++) {
						rd_aln_buf.push_back(seq[j]);
					}
					if(op_m == 0) {
						for(size_t j = rdoff; j < rdoff + run_comb; j++) {
							rf_aln_buf.push_back(seq[j]);
						}
					} else {
						assert(run_m == run_comb);
						for(size_t j = st_m; j < st_m + run_m; j++) {
							rf_aln_buf.push_back(mdz_char[j]);
						}
						// TODO: assert that all added characters mismatch
					}
					mdrun += run_comb;
					rdoff += run_comb;
					if(run_comb < run_m) {
						assert(op_m == 0);
						mdz_oro[mdo].run -= run_comb;
					} else {
						mdo++;
					}
				}
			} else if(cop == 'I') {
				for(size_t j = rdoff; j < rdoff + crun; j++) {
					rd_aln_buf.push_back(seq[j]);
					rf_aln_buf.push_back('-');
				}
				rdoff += crun;
			} else if(cop == 'D') {
				char op_m = mdz_oro[mdo].op;
				int run_m = mdz_oro[mdo].run;
				int st_m = mdz_oro[mdo].offset;
				assert(op_m == 2);
				assert(crun == run_m);
				assert(run_m == crun);
				mdo++;
				for(size_t j = st_m; j < st_m + run_m; j++) {
					rd_aln_buf.push_back('-');
					rf_aln_buf.push_back(mdz_char[j]);
				}
			} else if(cop == 'N') {
				for(size_t j = 0; j < crun; j++) {
					rd_aln_buf.push_back('-');
					rf_aln_buf.push_back('-');
				}
			} else if(cop == 'S') {
				rdoff += crun;
			} else if(cop == 'H') {
				// pass
			} else if(cop == 'P') {
				throw 1;
			} else if(cop == '=') {
				throw 1;
			} else if(cop == 'X') {
				throw 1;
			} else {
				throw 1;
			}
		}
		assert(mdo == mdz_oro.size());
		rd_aln_buf.push_back(0);
		rf_aln_buf.push_back(0);
	}
	
	char *rest_of_line;
	bool valid;
	char *qname;
	int flag;
	char *rname;
	int pos;
	int mapq;
	char *cigar;
	char *rnext;
	int pnext;
	int tlen;
	char *seq;
	size_t len;
	char *qual;
	char *mdz;
	bool cigar_equal_x;
	int best_score;
	int left_clip;
	int right_clip;
	
	// For holding stacked alignment result
	EList<char> rf_aln_buf;
	EList<char> rd_aln_buf;
	
	// For holding edit transcript
	EList<char> edit_xscript;

	// For holding cigar parsing info
	EList<int> cigar_run;
	EList<char> cigar_ops;
	
	// For holding MD:Z parsing info
	EList<OpRunOffset> mdz_oro;
	EList<char> mdz_char;
};

/**
 * strtok already used to parse up to rname.  Parse the rest and return char *
 * to the extra flags.
 *
 * This function uses strtok so it will interrrupt any strtoks already in
 * progress.  But when it returns, it's finished with its stateful strtok
 * session, so there's nothing keeping the caller from starting another
 * strtok.
 */
static char * parse_from_rname_on(Alignment& al) {
	assert(al.rest_of_line != NULL);
	al.rname = strtok(al.rest_of_line, "\t"); assert(al.rname != NULL);
	char *pos_str = strtok(NULL, "\t"); assert(pos_str != NULL);
	al.pos = atoi(pos_str);
	char *mapq_str = strtok(NULL, "\t"); assert(mapq_str != NULL);
	al.mapq = atoi(mapq_str);
	assert(al.mapq < 256);
	al.cigar = strtok(NULL, "\t"); assert(al.cigar != NULL);
	al.parse_cigar();
	al.rnext = strtok(NULL, "\t"); assert(al.rnext != NULL);
	char *pnext_str = strtok(NULL, "\t"); assert(pnext_str != NULL);
	al.pnext = atoi(pnext_str);
	char *tlen_str = strtok(NULL, "\t"); assert(tlen_str != NULL);
	al.tlen = atoi(tlen_str);
	al.seq = strtok(NULL, "\t"); assert(al.seq != NULL);
	al.len = strlen(al.seq);
	al.qual = strtok(NULL, "\t"); assert(al.qual != NULL);
	al.rest_of_line = al.qual + strlen(al.qual) + 1;
	return al.rest_of_line;
}

/**
 * No guarantee about state of strtok upon return.
 */
static void print_unpaired(Alignment& al, size_t ordlen, FILE *fh_model, FILE *fh_recs) {
	assert(al.is_aligned());
	char *extra = parse_from_rname_on(al);
	char *ztz_tok = al.parse_extra(extra);
	al.best_score = atoi(ztz_tok);
	char fw_flag = al.is_fw() ? 'T' : 'F';
	
	/* TODO: add correct/incorrect info */
	
	if(fh_model != NULL) {
		/* Output information relevant to input model */
		fprintf(fh_model, "%d,%c,%s,%u,%c,%u,%s\n",
				al.best_score,
				fw_flag,
				al.qual,
				(unsigned)al.len,
				al.mate_flag(),
				(unsigned)ordlen,
				al.edit_xscript.ptr());
	}
	
	if(fh_recs != NULL) {
		/* Output information relevant to MAPQ model */
		fprintf(fh_recs, "%u,%d,%d",
				(unsigned)al.len,
				al.mapq,
				al.tlen);
		
		/* ... including all the ZT:Z fields */
		while(ztz_tok != NULL) {
			fprintf(fh_recs, ",%s", ztz_tok);
			ztz_tok = strtok(NULL, ",");
		}
	}
}

/**
 * No guarantee about state of strtok upon return.
 */
static void print_paired(Alignment& al1, struct Alignment& al2, FILE *fh_model, FILE *fh_recs) {
	assert(al1.is_aligned());
	assert(al2.is_aligned());
	
	char *extra1 = parse_from_rname_on(al1);
	char *extra2 = parse_from_rname_on(al2);
	
	char *ztz_tok1 = al1.parse_extra(extra1);
	al1.best_score = atoi(ztz_tok1);
	char fw_flag1 = al1.is_fw() ? 'T' : 'F';
	
	int fraglen = 0;
	int upstream1 = 0;
	
	/* TODO: add correct/incorrect info */
	
	if(fh_recs != NULL) {
		
		/*
		 * Mate 1
		 */
		
		/* Output information relevant to input model */
		fprintf(fh_recs, "%u,%d,%u,%d,%d",
				(unsigned)al1.len,
				al1.mapq,
				(unsigned)al2.len,
				al2.mapq,
				fraglen);
		/* ... including all the ZT:Z fields */
		while(ztz_tok1 != NULL) {
			size_t toklen = strlen(ztz_tok1);
			/* remove trailing whitespace */
			while(ztz_tok1[toklen-1] == '\n' || ztz_tok1[toklen-1] == '\r') {
				ztz_tok1[toklen-1] = '\0';
				toklen--;
			}
			fprintf(fh_recs, ",%s", ztz_tok1);
			ztz_tok1 = strtok(NULL, ",");
		}
	}
	
	char *ztz_tok2 = al2.parse_extra(extra2);
	al2.best_score = atoi(ztz_tok2);
	char fw_flag2 = al2.is_fw() ? 'T' : 'F';
	
	if(fh_recs != NULL) {
		
		/*
		 * Mate 2
		 */
		
		/* Output information relevant to input model */
		fprintf(fh_recs, ",%u,%d,%u,%d,%d",
				(unsigned)al2.len,
				al2.mapq,
				(unsigned)al1.len,
				al1.mapq,
				fraglen);
		/* ... including all the ZT:Z fields */
		while(ztz_tok2 != NULL) {
			fprintf(fh_recs, ",%s", ztz_tok2);
			ztz_tok2 = strtok(NULL, ",");
		}
	}
	
	if(fh_model != NULL) {
		/* Output information relevant to input model */
		fprintf(fh_model, "%d,%c,%s,%d,%u,%s,%c,%s,%d,%u,%s,%d,%d\n",
				al1.best_score + al2.best_score,
				fw_flag1,
				al1.qual,
				al1.best_score,
				(unsigned)al1.len,
				al1.edit_xscript.ptr(),
				fw_flag2,
				al2.qual,
				al2.best_score,
				(unsigned)al2.len,
				al2.edit_xscript.ptr(),
				fraglen,
				upstream1);
	}
}

/**
 * Read the input SAM file while simultaneously writing out records used to
 * train a MAPQ model as well as records used to build an input model.
 */
static int sam_pass1(FILE *fh,
					 FILE *orec_u_fh, FILE *omod_u_fh,
					 FILE *orec_b_fh, FILE *omod_b_fh,
					 FILE *orec_c_fh, FILE *omod_c_fh,
					 FILE *orec_d_fh, FILE *omod_d_fh,
					 bool quiet)
{
	/* Advise the kernel of our access pattern.  */
	/* posix_fadvise(fd, 0, 0, 1); */ /* FDADVICE_SEQUENTIAL */
	
	char linebuf1[BUFSZ], linebuf2[BUFSZ];
	int line1 = 1;
	
	Alignment al1, al2;
	
	int al_cur1 = 1;
	
	int nline = 0, nignored = 0, npair = 0, nunp = 0;
	int nunp_al = 0, nunp_unal = 0, npair_badend = 0, npair_conc = 0, npair_disc = 0, npair_unal = 0;
	
	while(1) {
		char *line = line1 ? linebuf1 : linebuf2;
		if(fgets(line, BUFSZ, fh) == NULL) {
			break; /* done */
		}
		if(line[0] == '@') {
			continue; // skip header
		}
		nline++;
		char *qname = strtok(line, "\t"); assert(qname != NULL);
		assert(qname == line);
		char *flag_str = strtok(NULL, "\t"); assert(flag_str != NULL);
		int flag = atoi(flag_str);
		if((flag & 2048) != 0) {
			nignored++;
			continue;
		}
		
		/* switch which buffer "line" points to */
		line1 = !line1;
		
		/* somehow switch between alignments? */
		Alignment& al_cur  = al_cur1 ? al1 : al2;
		assert(!al_cur.valid);
		al_cur.clear();
		Alignment& al_prev = al_cur1 ? al2 : al1;
		al_cur1 = !al_cur1;
		
		al_cur.rest_of_line = flag_str + strlen(flag_str) + 1; /* for re-parsing */
		al_cur.qname = qname;
		al_cur.flag = flag;
		
		/* If we're able to mate up ends at this time, do it */
		Alignment *mate1 = NULL, *mate2 = NULL;
		if(al_cur.mate_flag() != '0' && al_prev.valid) {
			if(al_cur.mate_flag() == '1') {
				assert(al_prev.mate_flag() == '2');
				mate1 = &al_cur;
				mate2 = &al_prev;
			} else {
				assert(al_cur.mate_flag() == '2');
				assert(al_prev.mate_flag() == '1');
				mate1 = &al_prev;
				mate2 = &al_cur;
			}
			mate1->valid = mate2->valid = false;
			npair++;
		}
		
		if(al_cur.mate_flag() == '0') {
			nunp++;
			
			/* Case 1: Current read is unpaired and unlineigned, we can safely skip */
			if(!al_cur.is_aligned()) {
				nunp_unal++;
				continue;
			}
			
			/* Case 2: Current read is unpaired and aligned */
			else {
				nunp_al++;
				print_unpaired(al_cur, 0, omod_u_fh, orec_u_fh);
			}
		}
		
		else if(mate1 != NULL) {
			/* Case 3: Current read is paired and unlineigned, opposite mate is also unlineigned; nothing more to do! */
			assert(mate2 != NULL);
			if(!mate1->is_aligned() && !mate2->is_aligned()) {
				npair_unal++;
				continue;
			}
			
			/* Case 4: Current read is paired and aligned, opposite mate is unlineigned */
			/* Case 5: Current read is paired and unlineigned, opposite mate is aligned */
			/* we handle both here */
			else if(mate1->is_aligned() != mate2->is_aligned()) {
				npair_badend++;
				print_unpaired(
							   mate1->is_aligned() ? *mate1 : *mate2,
							   mate1->is_aligned() ? mate2->len : mate1->len,
							   omod_b_fh, orec_b_fh);
			}
			
			else {
				assert(mate1->is_concordant() == mate2->is_concordant());
				
				if(mate1->is_concordant()) {
					/* Case 6: Current read is paired and both mates aligned, concordantly */
					npair_conc++;
					print_paired(*mate1, *mate2, omod_c_fh, orec_c_fh);
				}
				
				else {
					/* Case 7: Current read is paired and both mates aligned, not condordantly */
					npair_disc++;
					print_paired(*mate1, *mate2, omod_d_fh, orec_d_fh);
				}
			}
		}
		
		else {
			/* This read is paired but we haven't seen the mate yet */
			assert(al_cur.mate_flag() != '0');
			al_cur.valid = true;
		}
	}
	
	if(!quiet) {
		fprintf(stderr, "%d lines\n", nline);
		fprintf(stderr, "%d ignored b/c secondary\n", nignored);
		fprintf(stderr, "%d unpaired\n", nunp);
		fprintf(stderr, "    %d aligned\n", nunp_al);
		fprintf(stderr, "    %d unaligned\n", nunp_unal);
		fprintf(stderr, "%d paired-end\n", npair);
		fprintf(stderr, "    %d concordant\n", npair_conc);
		fprintf(stderr, "    %d discordant\n", npair_disc);
		fprintf(stderr, "    %d bad-end\n", npair_badend);
		fprintf(stderr, "    %d unaligned\n", npair_unal);
	}
	
	fclose(fh);
	
	return 0;
}

/**
 * Caller gives path to one or more SAM files, then the final argument is a prefix where all the 
 */
int main(int argc, char **argv) {
	
	string orec_u_fn, omod_u_fn;
	string orec_b_fn, omod_b_fn;
	string orec_c_fn, omod_c_fn;
	string orec_d_fn, omod_d_fn;
	
	vector<string> sams;

	for(int i = 1; i < argc; i++) {
		if(i == argc-1) {
			string prefix = argv[i];
			
			orec_u_fn = prefix + string("_rec_u.csv");
			orec_b_fn = prefix + string("_rec_b.csv");
			orec_c_fn = prefix + string("_rec_c.csv");
			orec_d_fn = prefix + string("_rec_d.csv");

			omod_u_fn = prefix + string("_mod_u.csv");
			omod_b_fn = prefix + string("_mod_b.csv");
			omod_c_fn = prefix + string("_mod_c.csv");
			omod_d_fn = prefix + string("_mod_d.csv");
		} else {
			sams.push_back(string(argv[i]));
		}
	}
	
	if(sams.empty()) {
		cerr << "Usage: qsim_parse_input [sam]* [output prefix]" << endl;
		return 1;
	}
	
	//
	// Unpaired
	//
	
	char orec_u_buf[BUFSZ];
	FILE *orec_u_fh = fopen(orec_u_fn.c_str(), "wb");
	if(orec_u_fh == NULL) {
		cerr << "Could not open output record file \"" << orec_u_fn << "\"" << endl;
		return -1;
	}
	setvbuf(orec_u_fh, orec_u_buf, _IOFBF, BUFSZ);
	
	char omod_u_buf[BUFSZ];
	FILE *omod_u_fh = fopen(omod_u_fn.c_str(), "wb");
	if(omod_u_fh == NULL) {
		cerr << "Could not open output model file \"" << omod_u_fn << "\"" << endl;
		return -1;
	}
	setvbuf(omod_u_fh, omod_u_buf, _IOFBF, BUFSZ);
	
	//
	// Bad-end
	//
	
	char orec_b_buf[BUFSZ];
	FILE *orec_b_fh = fopen(orec_b_fn.c_str(), "wb");
	if(orec_b_fh == NULL) {
		cerr << "Could not open output record file \"" << orec_b_fn << "\"" << endl;
		return -1;
	}
	setvbuf(orec_b_fh, orec_b_buf, _IOFBF, BUFSZ);
	
	char omod_b_buf[BUFSZ];
	FILE *omod_b_fh = fopen(omod_b_fn.c_str(), "wb");
	if(omod_b_fh == NULL) {
		cerr << "Could not open output model file \"" << omod_b_fn << "\"" << endl;
		return -1;
	}
	setvbuf(omod_b_fh, omod_b_buf, _IOFBF, BUFSZ);
	
	//
	// Concordant
	//
	
	char orec_c_buf[BUFSZ];
	FILE *orec_c_fh = fopen(orec_c_fn.c_str(), "wb");
	if(orec_c_fh == NULL) {
		cerr << "Could not open output record file \"" << orec_c_fn << "\"" << endl;
		return -1;
	}
	setvbuf(orec_c_fh, orec_c_buf, _IOFBF, BUFSZ);
	
	char omod_c_buf[BUFSZ];
	FILE *omod_c_fh = fopen(omod_c_fn.c_str(), "wb");
	if(omod_c_fh == NULL) {
		cerr << "Could not open output model file \"" << omod_c_fn << "\"" << endl;
		return -1;
	}
	setvbuf(omod_c_fh, omod_c_buf, _IOFBF, BUFSZ);
	
	//
	// Discordant
	//
	
	char orec_d_buf[BUFSZ];
	FILE *orec_d_fh = fopen(orec_d_fn.c_str(), "wb");
	if(orec_d_fh == NULL) {
		cerr << "Could not open output record file \"" << orec_d_fn << "\"" << endl;
		return -1;
	}
	setvbuf(orec_d_fh, orec_d_buf, _IOFBF, BUFSZ);
	
	char omod_d_buf[BUFSZ];
	FILE *omod_d_fh = fopen(omod_d_fn.c_str(), "wb");
	if(omod_d_fh == NULL) {
		cerr << "Could not open output model file \"" << omod_d_fn << "\"" << endl;
		return -1;
	}
	setvbuf(omod_d_fh, omod_d_buf, _IOFBF, BUFSZ);

	char buf_input_sam[BUFSZ];
	for(size_t i = 0; i < sams.size(); i++) {
		FILE *fh = fopen(sams[i].c_str(), "rb");
		if(fh == NULL) {
			cerr << "Could not open input SAM file \"" << sams[i] << "\"" << endl;
			return -1;
		}
		setvbuf(fh, buf_input_sam, _IOFBF, BUFSZ);
		sam_pass1(fh,
				  orec_u_fh, omod_u_fh,
				  orec_b_fh, omod_b_fh,
				  orec_c_fh, omod_c_fh,
				  orec_d_fh, omod_d_fh,
				  false); // not quiet
		fclose(fh);
	}
	
	fclose(omod_u_fh);
	fclose(orec_u_fh);
	fclose(omod_b_fh);
	fclose(orec_b_fh);
	fclose(omod_c_fh);
	fclose(orec_c_fh);
	fclose(omod_d_fh);
	fclose(orec_d_fh);
}
