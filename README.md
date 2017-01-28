Qtip
====

[![Build Status](https://travis-ci.org/BenLangmead/qtip.svg?branch=master)](https://travis-ci.org/BenLangmead/qtip)

### Background

Read alignment is the first step in most sequencing data analyses. It is also a source of errors and interpretability problems. Repetitive genomes, algorithmic shortcuts, and genetic variation impede the aligner's ability to find a read's true point of origin. Aligners therefore report a mapping quality: the probability the reported point of origin for a read is incorrect.

Qtip is an accurate, aligner-agnostic tool for predicting mapping qualities that works by simulating a set of _tandem_ reads, similar to the input reads in important ways, but for which the true point of origin is known. Alignments of tandem reads are used to build a model for predicting mapping quality, which is then applied to the input-read alignments.
The model is automatically tailored to the alignment scenario at hand, allowing it to make accurate mapping-quality predictions across a range of read lengths, alignment parameters, genomes, and read aligners.

### Using Qtip

Qtip runs alongside an existing aligner, though the aligner requires modifications for Qtip to obtain the feature data it needs to make predictions.  We have already made these modifications for the popular Bowtie 2, BWA-MEM and SNAP tools.  See the `software` subdirectory for details.

### Building Qtip

If Qtip was cloned/extracted to a directory `$QTIP_HOME`, then:

    make -C $QTIP_HOME/src

### Testing Qtip

    make -C $QTIP_HOME/test

### Running Qtip

    $QTIP_HOME/qtip

### Qtip architecture

![Qtip flow diagram](images/qtip_flow.png)

More extensive documentation is coming soon.  For now, please see the usage message for `qtip`
