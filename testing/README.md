This is an experimental tool for determining if changes to PANDAseq affect the output. You will need Vala to compile it.

1. Install PANDAseq.
2. Install Vala (sudo apt-get install valac or sudo yum install vala vala-tools).
3. Compiled the binary using:

	valac pandaseq-diff.vala dynamic_loader.c --pkg pandaseq-2 -X -lltdl -g

If PANDAseq has been installed somewhere that pkg-config, and Vala are not looking, set the evironment variable:

	export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig

and pass the compiler flag

	--vapidir $PREFIX/share/vala/vapi

where $PREFIX is the prefix passed to `./configure` when building PANDAseq.

4. Modify PANDAseq and compile it.
5. Run PANDAseq-Diff on a sample dataset:

	./pandaseq-diff -l libpandaseq.la -f mcbath_1.fastq.bz2 -r mcbath_2.fastq.bz2

The input files must be BZipped. The supplied `libpandaseq.la` is the file built in the top-level directory. This file will be opened and each read pair will be assembled by both the existing and new assemblers, and the results compared.

The recommended data set is the sample dataset provided at http://neufeldserver.uwaterloo.ca/~apmasell/pandaseq_sampledata.tar or a small subset, http://neufeldserver.uwaterloo.ca/~apmasell/pandaseq_sampledata_small.tar that was used in the original publication.

This code makes use of strange dynamic loading behaviour and so requires that methods are not called on the assembler and that no ABI changes have occured.
