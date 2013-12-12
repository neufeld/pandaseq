![PANDASEQ](http://neufeldserver.uwaterloo.ca/~apmasell/pandaseq.svg)
========

PANDASEQ is a program to align Illumina reads, optionally with PCR primers embedded in the sequence, and reconstruct an overlapping sequence.

INSTALLATION
------------

[![Build Status](https://travis-ci.org/neufeld/pandaseq.png?branch=master)](https://travis-ci.org/neufeld/pandaseq) [![Build Status](https://travis-ci.org/neufeld/pandaseq-sam.png?branch=master)](https://travis-ci.org/neufeld/pandaseq-sam)

Binary packages are available for recent versions of Windows, MacOS and Linux. Installing from source is not too difficult. See [Installation instructions](https://github.com/neufeld/pandaseq/wiki/Installation) for details.

Development packages for zlib and libbz2 are needed, as is a standard compiler environment. On Ubuntu, this can be installed via

	sudo apt-get install build-essentials libtool automake zlib1g-dev libbz2-dev

On MacOS, the Apple Developer tools and Fink (or MacPorts or Brew) must be installed, then

	sudo fink install bzip2-dev

After the support packages are installed, one should be able to do:

	./autogen.sh && ./configure && make && sudo make install

If you receive an error that `libpandaseq.so.[number]` is not found on Linux, try running:

	sudo ldconfig

USAGE
-----

Please consult the manual page by invoking

	man pandaseq

or visiting [online PANDAseq manual page](http://neufeldserver.uwaterloo.ca/~apmasell/pandaseq_man1.html).

The short version is

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

A [Vala binding](http://neufeldserver.uwaterloo.ca/~apmasell/pandaseq-vapi/) is also included.

Other lanugage bindings are welcome.

FAQ
---

### Can I insist that PANDAseq only assembler perfect sequences?
Yes, but you shouldn't want to do it. The whole point is to fix sequences which are probably good. There is no quality setting that will achieve this effect. You can use the plugin `completely_miss_the_point`, but this really does miss the point. Moreover, assuming that the sequencer is right in the overlap region and in the non-overlapping regions requires an unsound leap in statistics.

### Can I use SAM/BAM files as input without converting them to FASTQ?
Yes. [PANDAseq-sam](https://github.com/neufeld/pandaseq-sam) extends PANDAseq to do this. SAM/BAM files do not guarantee that sequences will be in the right order, so files may be slower and PANDAseq will use more memory.

### The scores of the output bases seem really low. What's wrong?
Nothing. The quality scores of the output do not have any similarity to the original quality scores and are not uniform across the sequence (i.e., the overlap is scored differently from the unpaired ends.

In the overlap region where there is a mismatch, it is the probability that one base was sequenced correctly and the other was sequenced incorrectly. If both bases have high scores (i.e., are probably correct), the chance of the resulting base is low (i.e., is probably incorrect). For more information, see the paper. Also, remember that the PHRED to probability conversion is not linear, so most scores are relatively high. It's also not uncommon to see the PHRED score `!`, which is zero, but in this context, it means less than `"` (PHRED = 1, P = .20567).

Again, these scores are not meant to be interpreted as regular scores and should not be processed by downstream applications expecting PHRED scores from Illumina sequences.


### The scores of the non-overlapping regions are not the same as the original reads. Why?
The PHRED scores from the input are not copied directly to the output when using FASTQ (`-F`) output. They go through a transformation from PHRED scores into probabilities, which is how PANDAseq uses them. When output as FASTQ, the probability is converted back to a PHRED scores. The rounding error involved can cause a score to jump by one.

ALTERNATIVES
------------

[PEAR (Paired-End AssembleR)](http://www.exelixis-lab.org/pear)  
[FLASH (Fast Length Adjustment of SHort reads)](http://ccb.jhu.edu/software/FLASH/)  
[COPE (Connecting Overlapped Pair-End reads)](ftp://ftp.genomics.org.cn/pub/cope)  
[XORRO (Rapid Pair-end Read Overlapper)](http://arxiv.org/pdf/1304.4620v1.pdf)  

CITATION
--------

Andre P Masella, Andrea K Bartram, Jakub M Truszkowski, Daniel G Brown and Josh D Neufeld. PANDAseq: paired-end assembler for illumina sequences. BMC Bioinformatics 2012, 13:31. <http://www.biomedcentral.com/1471-2105/13/31>
