#include "api/BamReader.h"
#include <iostream>

using namespace BamTools;
using namespace std;

int main(void) {

    // at some point, start our merge operation
    string inputFilename = "../input.bam";

    // provide some input & output filenames
    // attempt to open our BamMultiReader
    BamReader reader;
    if(!reader.Open(inputFilename)) {
        cerr << "Could not open input BAM files." << endl;
        return -1;
    }

    // retrieve 'metadata' from BAM files, these are required by BamWriter
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();

    // iterate through all alignments
    BamAlignment al;
    while(reader.GetNextAlignment(al)) {
        string ztz;
        al.GetTag<string>("ZT", ztz);
        cerr << al.MapQuality << ":" << ztz << endl;
    }

    // close the reader
    reader.Close();
}
