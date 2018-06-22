![PANDASEQ](https://rawgithub.com/neufeld/pandaseq/master/pandaseq.svg)
========

PANDASEQ is a program to align Illumina reads, optionally with PCR primers embedded in the sequence, and reconstruct an overlapping sequence.

INSTALLATION
------------

[![Build Status](https://travis-ci.org/neufeld/pandaseq.png?branch=master)](https://travis-ci.org/neufeld/pandaseq)  

Binary packages are available for recent versions of Windows, MacOS and Linux. Source code is also available. See [Installation instructions](https://github.com/neufeld/pandaseq/wiki/Installation) for details.

Development packages for zlib and libbz2 are needed, as well as a standard compiler environment. On Ubuntu, this can be installed via:

	sudo apt-get install build-essential libtool automake zlib1g-dev libbz2-dev pkg-config

On MacOS, the Apple Developer tools and Fink (or MacPorts or Brew) must be installed, then:

	sudo fink install bzip2-dev pkgconfig

The newer AutoTools from Fink are needed over the ones provided by Apple, so ensure that Fink's `bin` directory precedes `/usr/bin` in the `$PATH`.

After the support packages are installed, one should be able to do:

	./autogen.sh && ./configure && make && sudo make install

If you receive an error that `libpandaseq.so.[number]` is not found on Linux, try running:

	sudo ldconfig

USAGE
-----

Please consult the manual page by invoking:

	man pandaseq

or visiting [online PANDAseq manual page](https://storage.googleapis.com/pandaseq/pandaseq.html).

The short version is:

	pandaseq -f forward.fastq -r reverse.fastq

REPORTING BUGS
--------------

Before filing a bug, consult [how to file a bug](https://github.com/neufeld/pandaseq/wiki/Filing-Bugs).

Please run:

	curl https://raw.github.com/neufeld/pandaseq/master/pandabug | sh

or

	wget -O- https://raw.github.com/neufeld/pandaseq/master/pandabug | sh

to create a header with basic details about your system. Please include:

1. The output of the above script.
2. The exact error message. If this is a compilation error, do not truncate the output. If this is a problem when assembling, keep the `INFO ARG` lines, and the last few lines, but you may truncate the middle.
3. If you have tried multiple different things, please list them all.
4. Your sequencing data may be requested. This usually does not necessitate all the reads.

BINDING
-------

PANDAseq may be used in other programs via a programmatic interface. Consult the header file `pandaseq.h` for more details. The C interface is pseudo-object oriented and documented in the header. The library provides `pkg-config` information, so compiling against it can be done using something like:

	cc mycode.c `pkg-config --cflags --libs pandaseq-2`

or using, in `configure.ac`:

	PKG_CHECK_MODULES(PANDASEQ, [ pandaseq-2 >= 2.5 ])

Other language bindings are welcome.

FAQ
---

### Can I insist that PANDAseq only assemble perfect sequences?
Yes, but you shouldn't want to do it. The whole point is to fix sequences which are probably good. There is no quality setting that will achieve this effect. You can use the plugin `completely_miss_the_point`, but this really does miss the point. Moreover, assuming that the sequencer is right in the overlap region and in the non-overlapping regions requires an unsound leap in statistics.

### Can I use SAM/BAM files as input without converting them to FASTQ?
Yes. [PANDAseq-sam](https://github.com/neufeld/pandaseq-sam) extends PANDAseq to do this. SAM/BAM files do not guarantee that sequences will be in the right order, so using SAM/BAM files may be slower and PANDAseq will use more memory.

### The scores of the output bases seem really low. What's wrong?
Nothing. The quality scores of the output do not have any similarity to the original quality scores and are not uniform across the sequence (i.e., the overlap is scored differently from the unpaired ends.

In the overlap region where there is a mismatch, it is probable that one base was sequenced correctly and the other was sequenced incorrectly. If both bases have high scores (i.e., are probably correct), the chance of the resulting base is low (i.e., is probably incorrect). For more information, see the paper. Also, remember that the PHRED to probability conversion is not linear, so most scores are relatively high. It's also not uncommon to see the PHRED score `!`, which is zero, but in this context, it means less than `"` (PHRED = 1, P = .20567).

Again, these scores are not meant to be interpreted as regular scores and should not be processed by downstream applications expecting PHRED scores from Illumina sequences.

### The scores of the non-overlapping regions are not the same as the original reads. Why?
The PHRED scores from the input are not copied directly to the output when using FASTQ (`-F`) output. They go through a transformation from PHRED scores into probabilities, which is what PANDAseq uses. When output as FASTQ, the probabilities are converted back to PHRED scores. The rounding error involved can cause a score to jump by one.

### How many sequences should there be in the output?
You should expect that PANDAseq will output fewer sequences than the read pairs given to it. The log contains several `STAT` lines that will help with the analysis. Lines containing `STAT READS` report the number of read pairs in the input. Sequences first go through a number of basic filtering steps and then user-specified filtering steps. If provided, forward and reverse primers are aligned and clipped. The optimal overlap is selected and the sequence is constructed. The quality score is verified and any user-specified filtering is done. Any of these steps might fail and cause the sequence to be rejected. For each of the possible rejection reasons, the log file will contain a `STAT` line reporting the number of sequences filtered, as is described in the _Output Statistics_ section of the manual.

If multiple threads are used, which the default on most platforms, each thread collects this information separately. The output log will output a group of `STAT` lines per thread.

The `STAT SLOW` line is informative; those sequences were not rejected. The other `STAT` lines (i.e., not `READS` or `SLOW`) should sum to the `STAT READS` line.

ALTERNATIVES
------------

Similar algorithms (i.e., determine the overlap, then fuse the reads):
 - [COPE (Connecting Overlapped Pair-End reads)](http://sourceforge.net/projects/coperead/) – Algorithm similar to FLASH
 - FastqJoin in [ea-utils](https://code.google.com/p/ea-utils/wiki/FastqJoin) – Algorithm included in PANDAseq
 - [FLASH (Fast Length Adjustment of SHort reads)](http://ccb.jhu.edu/software/FLASH/) – Algorithm included in PANDAseq
 - [PEAR (Paired-End AssembleR)](http://www.exelixis-lab.org/pear) – Algorithm included in PANDAseq
 - [stitch](https://github.com/audy/stitch)
 - [XORRO (Rapid Pair-end Read Overlapper)](http://arxiv.org/pdf/1304.4620v1.pdf)
 - [leeHom](http://nar.oxfordjournals.org/content/42/18/e141.full)

Completely different methods:
 - [SeqPrep](https://github.com/jstjohn/SeqPrep/) – Uses alignment

CITATION
--------

Andre P Masella, Andrea K Bartram, Jakub M Truszkowski, Daniel G Brown and Josh D Neufeld. PANDAseq: paired-end assembler for illumina sequences. BMC Bioinformatics 2012, 13:31. <http://www.biomedcentral.com/1471-2105/13/31>
