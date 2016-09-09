//
//  fasta.cpp
//  qtip
//
//  Created by Ben Langmead on 9/4/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#include "fasta.h"
#include <iostream>
#include <cassert>

using namespace std;

/**
 * Mapping from ASCII characters for ambiguous nucleotides into masks:
 */
int dna_upper[] = {
	/*   0 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/*  16 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/*  32 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/*  48 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/*  64 */ 'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/*  80 */ 'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/*  96 */ 'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 112 */ 'N', 'N', 'N', 'N', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 128 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 144 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 160 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 176 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 192 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 208 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 224 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	/* 240 */ 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

/**
 * Returns pointer to sequence.  Return value can be NULL.  If
 * done() is also true, then we're done.  Otherwise, caller can
 * call again.
 */
const char *FastaChunkwiseParser::next(
	std::string& refid, // out: name of reference buffer is from
	std::string& refid_full, // out: name of reference buffer is from
	size_t& refoff, // out: reference offset of first character in buffer
	size_t& retsz) // out: number of characters in returned buffer
{
	// fni_ and foff_ are set according to where we want to read our next
	// chunk, so foff_ == 0 means we're starting a new file
	if(done()) {
		return NULL;
	}
	if(foff_ == 0) {
		fh_ = fopen(fns_[fni_].c_str(), "rb");
		if(fh_ == NULL) {
			cerr << "Could not open FASTA file \"" << fns_[fni_] << "\"" << endl;
			throw 1;
		}
		setvbuf(fh_, fabuf_, _IOFBF, FASTA_BUFSZ);
	}
	
	// Move over
	char *buf = buf_;
	if(bufcur_ >= olap_) {
		for(size_t i = 0; i < olap_; i++) {
			buf_[i] = buf_[bufcur_-olap_+i];
		}
		buf += olap_;
		bufcur_ = olap_;
	}
	
	bool first = true;
	while(true) {
		//assert(fh_ != NULL);
		int c = get(fh_);
		if(c == EOF) {
			if(!first) {
				pushback(EOF);
				retsz = bufcur_ = buf - buf_;
				refoff = refoff_ - retsz;
				return buf_; // End of file
			}
			if(ferror(fh_)) {
				throw 1;
			}
			assert(feof(fh_));
			bufcur_ = foff_ = 0;
			fni_++;
			fclose(fh_); fh_ = NULL;
			return NULL;
		} else if(c == '>') {
			if(!first) {
				pushback('>');
				retsz = bufcur_ = buf - buf_;
				refoff = refoff_ - retsz;
				return buf_; // Beginning of new FASTA record
			}
			bufcur_ = 0;
			buf = buf_;
			refoff_ = 0;
			refid.clear();
			refid_full.clear();
			c = getc_unlocked(fh_);
			while(!isspace(c)) {
				refid.push_back(c);
				refid_full.push_back(c);
				c = getc_unlocked(fh_);
			}
			while(c != '\n' && c != '\r') {
				refid_full.push_back(c);
				c = getc_unlocked(fh_);
			}
		} else if(!isspace(c)) {
			first = false;
			*buf++ = (char)dna_upper[c];
			refoff_++;
			foff_++;
			if(buf - buf_ == chunksz_) {
				retsz = bufcur_ = chunksz_;
				refoff = refoff_ - retsz;
				return buf_; // Chunk buffer full
			}
			assert(buf - buf_ < chunksz_);
		}
	}
	assert(false);
	return NULL;
}

#ifdef FASTA_MAIN

#include <fstream>

static void test1() {
	string fn1 = ".test1.1.fa";
	string fn2 = ".test1.2.fa";
	string fn3 = ".test1.3.fa";
	
	ofstream ofs1(fn1, ofstream::out);
	ofs1 << ">record1 ok" << endl;
	ofs1 << "AAAACCCCGGGG" << endl;
	ofs1 << "TTTT" << endl;
	ofs1 << ">record2 mk" << endl;
	ofs1 << "A" << endl;
	ofs1 << "T" << endl;
	ofs1 << ">record3\tblah" << endl;
	ofs1 << "A";
	ofs1.close();

	ofstream ofs2(fn2, ofstream::out);
	ofs2 << "" << endl;
	ofs2 << ">record4 ok" << endl;
	ofs2 << "TG" << endl;
	ofs2.close();

	ofstream ofs3(fn3, ofstream::out);
	ofs3 << "" << endl;
	ofs3 << ">record5 ok" << endl;
	ofs3 << "CA" << endl;
	ofs3.close();

	vector<string> fns;
	fns.push_back(fn1);
	fns.push_back(fn2);
	fns.push_back(fn3);
	
	FastaChunkwiseParser fa(fns, 2, 1);
	string refid, refid_full;
	size_t refoff = std::numeric_limits<size_t>::max();
	size_t retsz = std::numeric_limits<size_t>::max();
	
	const char *buf;
	for(int i = 0; i < 3; i++) {
		buf = fa.next(refid, refid_full, refoff, retsz);
		assert(retsz == 2);
		assert(refid == string("record1"));
		assert(refid_full == string("record1 ok"));
		assert(refoff == i);
		assert(buf[0] == 'A');
		assert(buf[1] == 'A');
	}

	buf = fa.next(refid, refid_full, refoff, retsz);
	assert(retsz == 2);
	assert(refoff == 3);
	assert(buf[0] == 'A');
	assert(buf[1] == 'C');

	for(int i = 0; i < 3; i++) {
		buf = fa.next(refid, refid_full, refoff, retsz);
		assert(retsz == 2);
		assert(refoff == 4 + i);
		assert(buf[0] == 'C');
		assert(buf[1] == 'C');
	}

	buf = fa.next(refid, refid_full, refoff, retsz);
	assert(retsz == 2);
	assert(refoff == 7);
	assert(buf[0] == 'C');
	assert(buf[1] == 'G');

	for(int i = 0; i < 3; i++) {
		buf = fa.next(refid, refid_full, refoff, retsz);
		assert(retsz == 2);
		assert(refoff == 8 + i);
		assert(buf[0] == 'G');
		assert(buf[1] == 'G');
	}

	buf = fa.next(refid, refid_full, refoff, retsz);
	assert(retsz == 2);
	assert(refoff == 11);
	assert(buf[0] == 'G');
	assert(buf[1] == 'T');

	for(int i = 0; i < 3; i++) {
		buf = fa.next(refid, refid_full, refoff, retsz);
		assert(retsz == 2);
		assert(refoff == 12 + i);
		assert(buf[0] == 'T');
		assert(buf[1] == 'T');
	}
	
	assert(refid == string("record1"));
	assert(refid_full == string("record1 ok"));
	
	buf = fa.next(refid, refid_full, refoff, retsz);
	assert(refid == string("record2"));
	assert(refid_full == string("record2 mk"));
	assert(retsz == 2);
	assert(refoff == 0);
	assert(buf[0] == 'A');
	assert(buf[1] == 'T');

	buf = fa.next(refid, refid_full, refoff, retsz);
	assert(refid == string("record3"));
	assert(refid_full == string("record3\tblah"));
	assert(retsz == 1);
	assert(refoff == 0);
	assert(buf[0] == 'A');

	do {
		assert(!fa.done());
		buf = fa.next(refid, refid_full, refoff, retsz);
	} while(buf == NULL);
	
	assert(refid == string("record4"));
	assert(refid_full == string("record4 ok"));
	assert(retsz == 2);
	assert(refoff == 0);
	assert(buf[0] == 'T');
	assert(buf[1] == 'G');

	do {
		assert(!fa.done());
		buf = fa.next(refid, refid_full, refoff, retsz);
	} while(buf == NULL);

	assert(refid == string("record5"));
	assert(refid_full == string("record5 ok"));
	assert(retsz == 2);
	assert(refoff == 0);
	assert(buf[0] == 'C');
	assert(buf[1] == 'A');
}

int main(void) {
	test1();
	cout << "ALL TESTS PASSED" << endl;
}
#endif

