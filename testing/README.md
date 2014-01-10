This is an experimental tool for determining if changes to PANDAseq affect the output. You will need Vala to compile it.

1. Install PANDAseq.
2. Install Vala (`sudo apt-get install valac` or `sudo yum install vala vala-tools`).
3. Modify PANDAseq and compile it, but do not install.
3. Compiled the tester using `make`.

	If PANDAseq has been installed somewhere that pkg-config, and Vala are not looking, set the `PREFIX` in the Makefile.

5. Run PANDAseq-Diff on a sample dataset:

		./pandaseq-diff -f mcbath_1.fastq.bz2 -r mcbath_2.fastq.bz2

Each read pair will be assembled by both the existing and new assemblers, and the results compared.

The recommended data set is the sample dataset provided at [McBath dataset](http://neufeldserver.uwaterloo.ca/~apmasell/pandaseq_sampledata.tar) or [a small subset](http://neufeldserver.uwaterloo.ca/~apmasell/pandaseq_sampledata_small.tar) that was used in the original publication. You can use the `-W` option to download and use these sequences automatically.

This code makes use of strange `objcopy` behaviour and so requires that methods are not called on the assembler and that no ABI changes have occured. It also is probably very non-portable, but does work on Linux.

The tool will compile `access_new_library.vala` in such a way that it will only link against the library just compiled. The main file, `pandaseq-diff.vala` will be compiled against the installed library.
