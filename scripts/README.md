These scripts are used from various other places in the repo, including from the `Makefile`s under `data/simulated_data` and from the SNAP test scripts under `software/snap`.

* `mason_convert.py`: Convert Mason-formatted FASTQ files to the augmented `wgsim`-like formatting used by the simulation scripts
* `fastq_interleave.py`: Interleave two paired-end FASTQ files.  Sometimes useful for tools like BWA-MEM and SNAP that take interleaved FASTQ.
