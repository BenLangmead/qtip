//
//  qsim_rewrite.cpp
//  qsim
//
//  Created by Ben Langmead on 9/25/15.
//  Copyright (c) 2015 JHU. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include "qsim_rewrite.h"

using namespace std;

/**
 * Re-write a SAM file to include our new MAPQ prediction.
 */

bool write_orig_mapq = false;
const char *orig_mapq_flag = "Zm:Z:";

bool write_precise_mapq = false;
const char *precise_mapq_flag = "Zp:Z:";

const static size_t BUFSZ = 65536;

#define FILEDEC(fn, fh, buf, typ, do_open) \
	char buf [BUFSZ]; \
	FILE * fh = NULL; \
	if(do_open) { \
		fh = fopen( fn .c_str(), "wb"); \
		if(fh == NULL) { \
			cerr << "Could not open output " << typ << " file \"" << fn << "\"" << endl; \
			return -1; \
		} \
		setvbuf(fh, buf, _IOFBF, BUFSZ); \
	}

// Iterate through each line of the input SAM file in tandem with the lines of
// the prediction file.

int main(int argc, char **argv) {

	if(argc == 1) {
		// print which arguments from ts.py should pass through to here
		cout << "orig-mapq-flag "
		     << "precise-mapq-flag "
		     << "write-orig-mapq "
		     << "write-precise-mapq"
		     << endl;
		return 0;
	}

	string fn;
	string prefix;
	vector<string> sams;
	vector<string> preds;
	string osam_fn;
	
	// All arguments except last are SAM files to parse.  Final argument is
	// prefix for output files.
	{
		int section = 0, prefix_set = 0;
		for(int i = 1; i < argc; i++) {
			if(strcmp(argv[i], "--") == 0) {
				section++;
				continue;
			}
			if(section == 0) {
				if(strcmp(argv[i], "orig-mapq-flag") == 0) {
					if(i >= argc-1) {
						cerr << "Error: orig-mapq-flag requires argument" << endl;
						throw 1;
					}
					orig_mapq_flag = argv[++i];
				}
				if(strcmp(argv[i], "precise-mapq-flag") == 0) {
					if(i >= argc-1) {
						cerr << "Error: precise_mapq_flag requires argument" << endl;
						throw 1;
					}
					precise_mapq_flag = argv[++i];
				}
				if(strcmp(argv[i], "write-orig-mapq") == 0) {
					write_orig_mapq = true;
				}
				if(strcmp(argv[i], "write-precise-mapq") == 0) {
					write_precise_mapq = true;
				}
			} else if(section == 1) {
				sams.push_back(string(argv[i]));
			} else if(section == 2) {
				preds.push_back(string(argv[i]));
			} else {
				prefix = argv[i];
				prefix_set++;
				osam_fn = prefix + string(".sam");
			}
			if(prefix_set > 1) {
				cerr << "Warning: More than output prefix specified; using last one: \"" << prefix << "\"" << endl;
			}
		}
		if(sams.empty() || !prefix_set) {
			cerr << "Usage: qsim_rewrite" << endl;
		}
	}

	// Output SAM file
	FILEDEC(osam_fn, osam_fh, osam_buf, "SAM", true);

	char buf_input_sam[BUFSZ];
	char buf_input_pred[BUFSZ];

	for(size_t i = 0; i < sams.size(); i++) {
		// Input SAM file
		cerr << "Parsing SAM file \"" << sams[i] << "\"" << endl;
		FILE *fh_sam = fopen(sams[i].c_str(), "rb");
		if(fh_sam == NULL) {
			cerr << "Could not open input SAM file \"" << sams[i] << "\"" << endl;
			return -1;
		}
		setvbuf(fh_sam, buf_input_sam, _IOFBF, BUFSZ);
		
		// Input prediction file
		cerr << "Parsing prediction file \"" << preds[i] << "\"" << endl;
		FILE *fh_pred = fopen(preds[i].c_str(), "rb");
		if(fh_pred == NULL) {
			cerr << "Could not open input prediction file \"" << preds[i] << "\"" << endl;
			return -1;
		}
		setvbuf(fh_pred, buf_input_pred, _IOFBF, BUFSZ);
		
		// Do your thing
		
		fclose(fh_sam);
		fclose(fh_pred);
	}
}
