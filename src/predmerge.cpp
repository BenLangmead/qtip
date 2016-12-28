#include "predmerge.h"
#include <iostream>
#include <cassert>

using namespace std;

#define BUFSZ (64 * 1024)

/**
 * Construct new prediction merger; open files and read first prediction from
 * each.
 */
PredictionMerger::PredictionMerger(const vector<string>& in_fns) :
    in_fns_(in_fns), next_(-1)
{
    in_.resize(in_fns.size(), NULL);
    bufs_.resize(in_fns.size());
    preds_.resize(in_fns.size());
    done_.resize(in_fns.size(), false);
    for(size_t i = 0; i < in_fns.size(); i++) {
        in_[i] = fopen(in_fns[i].c_str(), "rb");
        bufs_[i].resize(BUFSZ);
        if(in_[i] == NULL) {
			cerr << "Could not open output file \"" << in_fns[i] << "\"" << endl;
            throw 1;
        }
        setvbuf(in_[i], &(bufs_[i][0]), _IOFBF, BUFSZ);
        advanceFile(i);
    }
}

/**
 * Close all open files.
 */
PredictionMerger::~PredictionMerger() {
	for(size_t i = 0; i < in_.size(); i++) {
		if(in_[i] != NULL) {
			fclose(in_[i]);
		}
		in_[i] = NULL;
	}
}

/**
 * Get and return next prediction.
 */
Prediction PredictionMerger::next() {
    if(next_ >= 0) {
        // Next file is known
        Prediction next_pred = preds_[next_];
        if(!advanceFile(next_)) {
			assert(preds_[next_].line >= next_pred.line);
			if(preds_[next_].line != next_pred.line + 1) {
				next_ = -1;
			}
		}
        return next_pred;
    } else {
        // Next file not known; have to check
        assert(next_ == -1);
        uint64_t min_line = numeric_limits<uint64_t>::max();
        int argmin = -1;
        const size_t insz = in_.size();
        for(size_t i = 0; i < insz; i++) {
            if(!done_[i]) {
                assert(preds_[i].valid());
                assert(preds_[i].line != numeric_limits<uint64_t>::max());
                assert(preds_[i].line != min_line);
                if(preds_[i].line < min_line) {
                    argmin = (int)i;
                    min_line = preds_[i].line;
                }
            }
        }
        if(argmin == -1) {
            // All input files exhausted
            return Prediction();
        }
        Prediction next_pred = preds_[argmin];
        if(!done_[argmin]) {
            if(!advanceFile(argmin)) {
				if(preds_[argmin].line == next_pred.line + 1) {
					next_ = argmin;
				}
			}
        }
        return next_pred;
    }
}

/**
 * Read the next prediction from one of the files.
 */
bool PredictionMerger::advanceFile(size_t i) {
    assert(!done_[i]);
    assert(!feof(in_[i]));
    double line = 0.0, mapq = 0.0;
    size_t ret = fread(&line, 8, 1, in_[i]);
	if(ret == 0 && feof(in_[i])) {
		done_[i] = true;
		preds_[i].reset();
		return false;
	}
    if(ret != 1) {
        cerr << "Could not read line id from prediction file \"" << in_fns_[i] << "\"; read " << ret << " items" << endl;
        throw 1;
    }
    assert(!feof(in_[i]));
    ret = fread(&mapq, 8, 1, in_[i]);
    if(ret != 1) {
        cerr << "Could not read MAPQ from prediction file \"" << in_fns_[i] << "\"; read " << ret << " items" << endl;
        throw 1;
    }
    preds_[i].line = (uint64_t)line;
    preds_[i].mapq = mapq;
    assert(mapq >= 0.0);
    assert(mapq <= 100.0);
    return true;
}

#ifdef PREDMERGE_MAIN

#include <fstream>

static void write2_or_throw(double d1, double d2, FILE *fh) {
    if(fwrite(&d1, 8, 1, fh) != 1) {
        cerr << "error writing d1" << endl;
        throw 1;
    }
    if(fwrite(&d2, 8, 1, fh) != 1) {
        cerr << "error writing d2" << endl;
        throw 1;
    }
}

static void write_file_a(string fn) {
	FILE *fh = fopen(fn.c_str(), "wb");
	if(fh == NULL) {
		cerr << "could not write test file a" << endl;
		throw 1;
	}
	write2_or_throw(0.0, 10.0, fh);
	write2_or_throw(2.0, 20.0, fh);
	write2_or_throw(3.0, 30.0, fh);
	write2_or_throw(10.0, 11.0, fh);
	write2_or_throw(12.0, 1.0, fh);
	fclose(fh);
}

static void write_file_b(string fn) {
	FILE *fh = fopen(fn.c_str(), "wb");
	if(fh == NULL) {
		cerr << "could not write test file b" << endl;
		throw 1;
	}
	write2_or_throw(1.0, 17.0, fh);
	write2_or_throw(4.0, 27.0, fh);
	write2_or_throw(6.0, 37.0, fh);
	write2_or_throw(11.0, 47.0, fh);
	write2_or_throw(14.0, 17.0, fh);
	write2_or_throw(15.0, 18.0, fh);
	fclose(fh);
}

static void write_file_c(string fn) {
	FILE *fh = fopen(fn.c_str(), "wb");
	if(fh == NULL) {
		cerr << "could not write test file c" << endl;
		throw 1;
	}
	write2_or_throw(5.0, 15.0, fh);
	write2_or_throw(7.0, 13.0, fh);
	write2_or_throw(8.0, 13.0, fh);
	write2_or_throw(9.0, 13.0, fh);
	fclose(fh);
}

static void test1() {
    string fn(".predmerge.test1.npy");
    write_file_a(fn);

    vector<string> fns;
    fns.push_back(fn);
    PredictionMerger m(fns);
    Prediction pred = m.next();
    assert(pred.line == 0);
    assert(pred.mapq == 10.0);
    pred = m.next();
    assert(pred.line == 2);
    assert(pred.mapq == 20.0);
    pred = m.next();
    assert(pred.line == 3);
    assert(pred.mapq == 30.0);
    pred = m.next();
    assert(pred.line == 10);
    assert(pred.mapq == 11.0);
    pred = m.next();
    assert(pred.line == 12);
    assert(pred.mapq == 1.0);
    pred = m.next();
    assert(!pred.valid());
}

static void test2() {
    string fn_a(".predmerge.test2.1.npy");
    string fn_b(".predmerge.test2.2.npy");
	write_file_a(fn_a);
	write_file_b(fn_b);

    vector<string> fns;
    fns.push_back(fn_a);
    fns.push_back(fn_b);
    PredictionMerger m(fns);
    Prediction pred = m.next();
    assert(pred.line == 0);
    assert(pred.mapq == 10.0);
    pred = m.next();
    assert(pred.line == 1);
    assert(pred.mapq == 17.0);
    pred = m.next();
    assert(pred.line == 2);
    assert(pred.mapq == 20.0);
    pred = m.next();
    assert(pred.line == 3);
    assert(pred.mapq == 30.0);
    pred = m.next();
    assert(pred.line == 4);
    assert(pred.mapq == 27.0);
    pred = m.next();
    assert(pred.line == 6);
    assert(pred.mapq == 37.0);
    pred = m.next();
    assert(pred.line == 10);
    assert(pred.mapq == 11.0);
    pred = m.next();
    assert(pred.line == 11);
    assert(pred.mapq == 47.0);
    pred = m.next();
    assert(pred.line == 12);
    assert(pred.mapq == 1.0);
    pred = m.next();
    assert(pred.line == 14);
    assert(pred.mapq == 17.0);
    pred = m.next();
    assert(pred.line == 15);
    assert(pred.mapq == 18.0);
    pred = m.next();
    assert(!pred.valid());
}

static void test3() {
    string fn_a(".predmerge.test3.1.npy");
    string fn_b(".predmerge.test3.2.npy");
    string fn_c(".predmerge.test3.3.npy");
	write_file_a(fn_a);
	write_file_b(fn_b);
	write_file_c(fn_c);

    vector<string> fns;
    fns.push_back(fn_a);
    fns.push_back(fn_b);
    fns.push_back(fn_c);
    PredictionMerger m(fns);
    Prediction pred = m.next();
    assert(pred.line == 0);
    assert(pred.mapq == 10.0);
    pred = m.next();
    assert(pred.line == 1);
    assert(pred.mapq == 17.0);
    pred = m.next();
    assert(pred.line == 2);
    assert(pred.mapq == 20.0);
    pred = m.next();
    assert(pred.line == 3);
    assert(pred.mapq == 30.0);
    pred = m.next();
    assert(pred.line == 4);
    assert(pred.mapq == 27.0);
    pred = m.next();
    assert(pred.line == 5);
    assert(pred.mapq == 15.0);
    pred = m.next();
    assert(pred.line == 6);
    assert(pred.mapq == 37.0);
    pred = m.next();
    assert(pred.line == 7);
    assert(pred.mapq == 13.0);
    pred = m.next();
    assert(pred.line == 8);
    assert(pred.mapq == 13.0);
    pred = m.next();
    assert(pred.line == 9);
    assert(pred.mapq == 13.0);
    pred = m.next();
    assert(pred.line == 10);
    assert(pred.mapq == 11.0);
    pred = m.next();
    assert(pred.line == 11);
    assert(pred.mapq == 47.0);
    pred = m.next();
    assert(pred.line == 12);
    assert(pred.mapq == 1.0);
    pred = m.next();
    assert(pred.line == 14);
    assert(pred.mapq == 17.0);
    pred = m.next();
    assert(pred.line == 15);
    assert(pred.mapq == 18.0);
    pred = m.next();
    assert(!pred.valid());
}

int main(void) {
	test1();
	test2();
	test3();
	cout << "ALL TESTS PASSED" << endl;
}
#endif

