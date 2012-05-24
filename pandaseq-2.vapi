/**
 * PANDAseq Illumina assembler
 *
 * Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
 */
[CCode(cheader_filename = "pandaseq.h")]
namespace Panda {
	/**
	 * Codes used for logging conditions during the assembly
	 *
	 * Some of these are errors and some are informational.
	 * The variadic arguments provided are listed.
	 */
	[CCode(cname = "PandaCode", has_type_id = false, cprefix = "PANDA_CODE_")]
	public enum Code {
		/**
		 * Show version of the module API.
		 *
		 * (int api_version)
		 */
		API_VERSION,
		/**
		 * Invalid character in nucleotide data
		 *
		 * (panda_seq_identifier *id, int character)
		 */
		BAD_NT,
		/**
		 * The best computed overlap
		 *
		 * (panda_seq_identifier *id, int overlap)
		 */
		BEST_OVERLAP,
		/**
		 * A nucleotide copied during reconstruction of the forward sequence
		 *
		 * (panda_seq_identifier *id, int index_in_assembly, int index_in_reverse, panda_result *nt)
		 */
		BUILD_FORWARD,
		/**
		 * A nucleotide determined during reconstruction of the overlap region
		 *
		 * (panda_seq_identifier *id, int index_in_assembly, int index_in_forward, int index_in_reverse, panda_result *nt, panda_nt *forward, panda_nt *reverse)
		 */
		BUILD_OVERLAP,
		/**
		 * A nucleotide copied during reconstruction of the reverse sequence
		 *
		 * (panda_seq_identifier *id, int index_in_assembly, int index_in_reverse, panda_result *nt)
		 */
		BUILD_REVERSE,
		/**
		 * A //k//-mer found in the forward sequence
		 *
		 * (panda_seq_identifier *id, unsigned int kmer, ssize_t position)
		 */
		FORWARD_KMER,
		/**
		 * An incorrect Illumina FASTQ header
		 *
		 *  (const char *id)
		 */
		ID_PARSE_FAILURE,
		/**
		 * The //k//-mer table is too small for the sequence
		 *
		 * (panda_seq_identifer *id)
		 */
		INSUFFICIENT_KMER_TABLE,
		/**
		 * A //k//-mer is thrown away due to a collision
		 *
		 * (panda_seq_identifier *id, unsigned int kmer, ssize_t position)
		 */
		LOST_KMER,
		/**
		 * A reconstruction is rejected because the quality score is below threshold
		 *
		 * (panda_seq_identifier *id, double quality, double threshold)
		 */
		LOW_QUALITY_REJECT,
		/**
		 * A pair of bases disagree in the reconstruction
		 *
		 * (panda_seq_identifier *id, int index_in_forward, int index_in_reverse, panda_qual *forward, panda_qual *reverse)
		 */
		MISMATCHED_BASE,
		/**
		 * Display information about a module
		 *
		 * (PandaModule module)
		 */
		MOD_INFO,
		/**
		 * The computed sequence length is invalid
		 *
		 * (panda_seq_identifier *id)
		 */
		NEGATIVE_SEQUENCE_LENGTH,
		/**
		 * No sequence data is availble in the file
		 *
		 * (panda_seq_identifier *id)
		 */
		NO_DATA,
		/**
		 * The file could not be opened
		 *
		 * (const char *filename)
		 */
		NO_FILE,
		/**
		 * The forward primer cannot be found in the sequence
		 *
		 * (panda_seq_identifier *id)
		 */
		NO_FORWARD_PRIMER,
		/**
		 * The quality information is missing in the FASTQ file
		 *
		 * (panda_seq_identifier *id)
		 */
		NO_QUALITY_INFO,
		/**
		 * The reverse primer cannot be found in the sequence
		 *
		 * (panda_seq_identifier *id)
		 */
		NO_REVERSE_PRIMER,
		/**
		 * The Illumina headers from the forward and reverse sequences do not match
		 *
		 * (panda_seq_identifier *foward, panda_seq_identifier *reverse)
		 */
		NOT_PAIRED,
		/**
		 * A possible overlap has been determined
		 *
		 * (panda_seq_identifier *id, size_t overlap, size_t matches, size_t mismatches, size_t unknowns, double probability)
		 */
		OVERLAP_POSSIBILITY,
		/**
		 * Error parsing FASTQ data
		 *
		 * (panda_seq_identifier *id)
		 */
		PARSE_FAILURE,
		/**
		 * An input file ended in the middle of a FASTQ file
		 *
		 * (panda_seq_identifier *id)
		 */
		PREMATURE_EOF,
		/**
		 * Reconstruction will commence with provided parameters
		 *
		 * (panda_seq_identifier *id, int overlap, int forward_unpaired, int reverse_unpaired)
		 */
		RECONSTRUCTION_PARAM,
		/**
		 * The number of sequences rejected by a particular module
		 *
		 * (PandaModule module, long rejected_sequence_count)
		 * @see Assembler.module_stats
		 */
		REJECT_STAT,
		/**
		 * A //k//-mer found in the reverse sequence
		 *
		 * (panda_seq_identifier *id, unsigned int kmer, ssize_t position)
		 */
		REVERSE_KMER,
		/**
		 * The reconsructed sequence will exceed the memory buffer
		 *
		 * (panda_seq_identifier *id)
		 */
		SEQUENCE_TOO_LONG,
	}

	/**
	 * A single nucleotide
	 */
	[CCode(cname = "panda_nt", has_type_id = false, cprefix = "PANDA_NT_")]
	[Flags]
	public enum Nt {
		/**
		 * Adenine
		 */
		A,
		/**
		 * Cytosine
		 */
		C,
		/**
		 * Guanine
		 */
		G,
		/**
		 * Thyamine
		 */
		T;
		/**
		 * Is nucleotide degenerate?
		 */
		[CCode(cname = "PANDA_NT_IS_DEGN")]
		public bool is_degenerate();
		/**
		 * Is nucleotide all possible values?
		 */
		[CCode(cname = "PANDA_NT_IS_N")]
		public bool is_n();
		/**
		 * Get the nucleotide code for an ASCII character in IUPAC
		 */
		[CCode(cname = "panda_nt_from_ascii")]
		public static Nt from_ascii(char c);
		/**
		 * Get the complement nucleotide code for an ASCII character in IUPAC
		 */
		[CCode(cname = "panda_nt_from_ascii_complement")]
		public static Nt from_ascii_complement(char c);
		/**
		 * Convert a nucleotide to an IUPAC representation
		 */
		[CCode(cname = "panda_nt_to_ascii")]
		public char to_ascii();
	}

	/**
	 * The manager for an assembly
	 */
	[CCode(cname = "struct panda_assembler", has_type_id = false, ref_func = "panda_assembler_ref", unref_func = "panda_assembler_unref")]
	public class Assembler {
		/**
		 * The number of sequences processed so far.
		 */
		public long count {
			[CCode(cname = "panda_assembler_get_count")]
			get;
		}

		/**
		 * The number of sequences rejected because they contain degenerate (N) bases.
		 */
		public long degenerate_count {
			[CCode(cname = "panda_assembler_get_degenerate_count")]
			get;
		}

		/**
		 * Reject sequences with degenerate (N) bases.
		 */
		public bool disallow_degenerates {
			[CCode(cname = "panda_assembler_get_disallow_degenerates")]
			get;
			[CCode(cname = "panda_assembler_set_disallow_degenerates")]
			set;
		}

		/**
		 * The minimum error estimation in the sequence data (epsilon)
		 */
		public double error_estimation {
			[CCode(cname = "panda_assembler_get_error_estimation")]
			get;
			[CCode(cname = "panda_assembler_set_error_estimation")]
			set;
		}

		/**
		 * The number of sequences rejected because the overlap could not be determined.
		 */
		public long failed_alignment_count {
			[CCode(cname = "panda_assembler_get_failed_alignment_count")]
			get;
		}

		/**
		 * The forward primer sequence to be stripped
		 *
		 * This is mutually exclusive with {@link forward_trim}
		 */
		public Nt[]? forward_primer {
			[CCode(cname = "panda_assembler_get_forward_primer")]
			get;
			[CCode(cname = "panda_assembler_set_forward_primer")]
			set;
		}

		/**
		 * The amount of forward sequence to strip
		 *
		 * This is mutually exclusive with {@link forward_primer}
		 */
		public size_t forward_trim {
			[CCode(cname = "panda_assembler_get_forward_trim")]
			get;
			[CCode(cname = "panda_assembler_set_forward_trim")]
			set;
		}

		/**
		 * The number of sequences rejected because the quality score is too low.
		 */
		public long low_quality_count {
			[CCode(cname = "panda_assembler_get_low_quality_count")]
			get;
		}

		/**
		 * The minimum overlap two sequences must have to be accepted. It must be greater than one.
		 */
		public int minimum_overlap {
			[CCode(cname = "panda_assembler_get_minimum_overlap")]
			get;
			[CCode(cname = "panda_assembler_set_minimum_overlap")]
			set;
		}

		/**
		 * The number of sequences rejected because the forward primer could not be aligned.
		 */
		public long no_forward_primer_count {
			[CCode(cname = "panda_assembler_get_no_forward_primer_count")]
			get;
		}

		/**
		 * The number of sequences rejected because the reverse primer could not be aligned.
		 */
		public long no_reverse_primer_count {
			[CCode(cname = "panda_assembler_get_no_reverse_primer_count")]
			get;
		}

		/**
		 * The number of sequences accepted.
		 */
		public long ok_count {
			[CCode(cname = "panda_assembler_get_ok_count")]
			get;
		}

		/**
		 * The reverse primer sequence to be stripped
		 *
		 * This is mutually exclusive with {@link reverse_trim}
		 */
		public Nt[]? reverse_primer {
			[CCode(cname = "panda_assembler_get_reverse_primer")]
			get;
			[CCode(cname = "panda_assembler_set_reverse_primer")]
			set;
		}

		/**
		 * The amount of reverse sequence to strip
		 *
		 * This is mutually exclusive with {@link reverse_primer}
		 */
		public size_t reverse_trim {
			[CCode(cname = "panda_assembler_get_reverse_trim")]
			get;
			[CCode(cname = "panda_assembler_set_reverse_trim")]
			set;
		}

		/**
		* The minimum quality threshold to have an assembly accepted. Must be between 0 and 1, exclusive.
		*/
		public double threshold {
			[CCode(cname = "panda_assembler_get_threshold")]
			get;
			[CCode(cname = "panda_assembler_set_threshold")]
			set;
		}

		/**
		 * Open a pair of bzipped for assembly.
		 *
		 * @param qualmin the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
		 */
		[CCode(cname = "panda_assembler_open_bz2")]
		public static Assembler? open_bz2(string forward, string reverse, owned Logger logger, uint8 qualmin = 33);

		/**
		 * Open a pair of gzipped (or uncompressed files) for assembly.
		 *
		 * @param qualmin the value to strip from the quality scores. Usually 33 or 64, depending on CASAVA version.
		 */
		[CCode(cname = "panda_assembler_open_gz")]
		public static Assembler? open_gz(string forward, string reverse, owned Logger logger, uint8 qualmin = 33);

		/**
		 * Create a new assembler from a sequence source.
		 *
		 * @param next the function to call to get the next sequence. The assembler does not manage the memory of the returned arrays, but assume it may use them until the next call or destruction. If not provided, you may use {@link assemble} but not {@link next}
		 * @param logger the function to call to report information to the user.
		 */
		[CCode(cname = "panda_assembler_new")]
		public Assembler(owned NextSeq? next, owned Logger logger);

		/**
		 * Create a new assembler for given to FASTQ streams.
		 * @see create_fastq_reader
		 */
		[CCode(cname = "panda_assembler_new_fastq_reader")]
		public Assembler.fastq (owned NextChar forward, owned NextChar reverse, owned Logger logger, uint8 qualmin = 33);

		/**
		 * Assemble a single sequence pair not drawn from the sequence stream.
		 *
		 * This works exactly like {@link next}, but instead of asking the {@link NextSeq} for the data, it expects this information to be provided.
		 */
		[CCode(cname = "panda_assembler_assemble")]
		public unowned result_seq? assemble(identifier id, qual[] forward, qual[] reverse);

		/**
		 * Add a module to this assembly process.
		 *
		 * Sequences will be checked using this module.
		 */
		[CCode(cname = "panda_assembler_add_module")]
		public void add_module(Module module);

		/**
		 * Log the number of sequences rejected by each module.
		 */
		[CCode(cname = "panda_assembler_module_stats")]
		public void module_stats();

		/**
		 * Assemble the next sequence from the input
		 *
		 * This function will process sequences until one is assembled successfully or no more sequences are available from the input stream, after which it will return null.
		 * The returned sequence becomes invalid after the next call.
		 */
		[CCode(cname = "panda_assembler_next")]
		public unowned result_seq? next();

		/**
		 * Increment reference count.
		 */
		[CCode(cname = "panda_assembler_ref")]
		public unowned Assembler @ref();
		/**
		 * Decrement reference count.
		 */
		[CCode(cname = "panda_assembler_unref")]
		public void unref();
	}

	/**
	 * Sequence validity checker
	 */
	[CCode(cname = "struct panda_module", ref_function = "panda_module_ref", unref_function = "panda_module_unref")]
	public class Module {
		/**
		 * The current module API version
		 */
		[CCode(cname = "PANDA_API")]
		public const int API;

		/**
		 * The API version of a module.
		 *
		 * This is only appropriate for loaded modules. Modules constructed by {@link create} will always return {@link API}.
		 */
		public int api {
			[CCode(cname = "panda_module_get_api")]
			get;
		}

		/**
		 * The arguments passed on loading of a module of a module.
		 *
		 * This is only appropriate for loaded modules.
		 */
		public string? args {
			[CCode(cname = "panda_module_get_args")]
			get;
		}

		/**
		 * The description of a module.
		 *
		 * This is only appropriate for loaded modules.
		 */
		public string? description {
			[CCode(cname = "panda_module_get_description")]
			get;
		}

		/**
		 * The name of a module.
		 */
		public string name {
			[CCode(cname = "panda_module_get_name")]
			get;
		}

		/**
		 * The usage information (i.e., help text) of a module.
		 *
		 * This is only appropriate for loaded modules.
		 */
		public string? usage {
			[CCode(cname = "panda_module_get_usage")]
			get;
		}

		/**
		 * The version of a module.
		 *
		 * This is only appropriate for loaded modules.
		 */
		public string? version {
			[CCode(cname = "panda_module_get_version")]
			get;
		}

		/**
		 * Check a sequence after reconstruction for validity.
		 */
		[CCode(cname = "PandaCheck", has_target = false, simple_generics = true)]
		public delegate bool Check<T>(T data, result_seq sequence);
		/**
		 * Check a sequence before reconstruction for validity.
		 */
		[CCode(cname = "PandaPreCheck", has_target = false, simple_generics = true)]
		public delegate bool PreCheck<T>(T data, identifier id, qual[] forward, qual[] reverse);

		/**
		 * Create a module given sequence checking parameters.
		 *
		 * @param name the name of the module, for user interaction
		 * @param check the check function, which must not be null
		 * @param precheck an optional check to be done before the module
		 */
		[CCode(cname = "panda_module_new", simple_generics = true)]
		public static Module create<T>(string name, Check<T> check, PreCheck<T> precheck, owned T user);

		/**
		 * Load a module from a string containg the module name and arguments.
		 *
		 * @param path the name or path to a module separated by LT_PATHSEP_CHAR and any arguments to the initialisation function of that module
		 */
		[CCode(cname = "panda_module_load")]
		public static Module load(string path);

		/**
		 * Increment reference count.
		 */
		[CCode(cname = "panda_module_ref")]
		public unowned Module @ref();
		/**
		 * Decrement reference count.
		 */
		[CCode(cname = "panda_module_unref")]
		public void unref();
	}

	/**
	 * Illumina sequence information from the FASTQ header
	 */
	[CCode(cname = "panda_seq_identifier", has_type_id = false, destroy_function = "")]
	public struct identifier {
		string flowcell;
		string instrument;
		int lane;
		int run;
		string tag;
		int tile;
		int x;
		int y;
		/**
		 * Parse an Illumina header
		 *
		 * Fills `id` with the parse result. The function returns the direction of the sequence (1 for forward, 2 for reverse) or 0 if an error occurs.
		 */
		[CCode(cname = "panda_seqid_parse")]
		public static int parse(out identifier id, string input);
		/**
		 * Compare two Illumina headers
		 */
		[CCode(cname = "panda_seqid_equal")]
		public bool equal(identifier other);
		/**
		 * Write the Illumina header to a printf-like function
		 */
		[CCode(cname = "panda_seqid_xprint")]
		public void print(PrintfFunc func);
		/**
		 * Write an Illumina header for a sequence identifer to a file
		 */
		[CCode(cname = "panda_seqid_print")]
		public void to_file(
#if POSIX
Posix.FILE
#else
GLib.FileStream
#endif
file);
		/**
		 * Create an Illumina header for a sequence identifier
		 *
		 * The return string must not be freed and subsequent calls will obliterate the previously returned string.
		 */
		[CCode(cname = "panda_seqid_str")]
		 public unowned string to_string();
	}

	/**
	 * A single nucleotide with quality information
	 */
	[CCode(cname = "panda_qual", has_type_id = false)]
	public struct qual {
		/**
		 * The nucleotide
		 */
		Nt nt;
		/**
		 * The quality score as a PHRED score
		 */
		char qual;
		/**
		 * Convert the PHRED quality score to a log probability.
		 */
		public double log_probability {
			[CCode(cname = "panda_quality_log_probability")]
			get;
		}
		/**
		 * Convert the PHRED quality score to a probability.
		 */
		public double probability {
			[CCode(cname = "panda_quality_probability")]
			get;
		}
	}

	/**
	 * A reconstructed nucleotide
	 */
	[CCode(cname = "panda_result", has_type_id = false)]
	public struct result {
		/**
		 * The nucleotide
		 */
		Nt nt;
		/**
		 * The quality score as a log probability
		 */
		double p;
	}

	/**
	 * A reconstructed sequence with meta information
	 */
	[CCode(cname = "panda_result_seq", has_type_id = false, destroy_function = "")]
	public struct result_seq {
		/**
		 * Number of uncalled bases in the sequence.
		 */
		size_t degenerates;

		/**
		 * The original forward sequence
		 */
		[CCode(array_length_cname = "forward_length")]
		qual[] forward;

		/**
		 * The number of nucleotides clipped from the forward sequence
		 */
		size_t forward_offset;

		/**
		 * The sequence identification information
		 */
		identifier name;

		/**
		 * Calculated quality score as the geometric mean of the product of the Illumina quality scores of the included bases.
		 *
		 * It will always be between 0 and 1.
		 */
		double quality;

		/**
		 * The reconstructed sequence with quality information
		 */
		[CCode(array_length_cname = "sequence_length")]
		result[] sequence;

		/**
		 * The original reverse sequence
		 */
		[CCode(array_length_cname = "reverse_length")]
		qual[] reverse;

		/**
		 * The number of nucleotides clipped from the reverse sequence
		 */
		size_t reverse_offset;

		/**
		 * The number of mismatches in the overlap region.
		 */
		size_t overlap_mismatches;

		/**
		 * Write an assembly to a FASTA file.
		 */
		[CCode(cname = "panda_output_fasta")]
		public bool write_fasta(
#if POSIX
Posix.FILE
#else
GLib.FileStream
#endif
file);
		/**
		 * Write an assembly to a FASTQ file.
		 */
		[CCode(cname = "panda_output_fastq")]
		public bool write_fastq(
#if POSIX
Posix.FILE
#else
GLib.FileStream
#endif
file);
	}

	/**
	 * Maximum length of a sequence
	 */
	[CCode(cname = "PANDA_MAX_LEN")]
	public const int MAX_LEN;

	/**
	 * Log an error/event
	 *
	 * If the function returns false, assembly will be halted.
	 * The variadic arguments will provide context based on the particular code passed.
	 * @see Code
	 */
	[CCode(cname = "PandaLogger", instance_pos = 0)]
	public delegate bool Logger(Code code, ...);

	/**
	 * Get the next character from a FASTQ file or EOF.
	 *
	 * For assembly from an alternate source of data, this function returns the next character in the stream.
	 */
	[CCode(cname = "PandaNextChar")]
	public delegate int NextChar();

	/**
	 * Get the next sequence pair.
	 *
	 * For assembly from a non-FASTQ text source, this function can provide the next sequence. The function must provide the sequences and metadata for assembly by modifing the values of its parameters.
	 * @param id the identifier information for the sequence pair
	 * @param forward the location of the parsed sequence data of the forward read.
	 * @param reverse the location of the parsed sequence data of the reverse read.
	 */
	[CCode(cname = "PandaNextSeq")]
	public delegate bool NextSeq(out identifier id, out unowned qual[] forward, out unowned qual[] reverse);

	[CCode(cname = "PandaPrintf", instance_pos = 0)]
	[PrintfFormat]
	public delegate void PrintfFunc(string format, ...);

	/**
	 * Create an object to read sequences from two character streams of FASTQ data
	 *
	 * @param forward the stream of forward characters, called every time a new character is required.
	 * @param reverse the stream of reverse characters, called every time a new character is required.
	 * @param logger the logging function to use during assembly.
	 */
	[CCode(cname = "panda_assembler_create_fastq_reader")]
	public static NextSeq create_fastq_reader(owned NextChar forward, owned NextChar reverse, Logger logger, uint8 qualmin = 33);

	[CCode(cname = "panda_version")]
	public unowned string get_version();

	/**
	 * Write errors and information to a file.
	 */
	[CCode(cname = "panda_logger_file")]
	public bool file_logger(
#if POSIX
Posix.FILE
#else
GLib.FileStream
#endif
file, Code code, ...);

	/**
	 * Writes an error message using the supplied printf-like function.
	 */
	[CCode(cname = "panda_logger_v")]
	bool logger_v(PrintfFunc func, Code code, va_list va);
	/**
	 * Create a file logger
	 */
	[CCode(cname = "PANDA_LOGGER")]
	public Logger create_logger(
#if POSIX
Posix.FILE
#else
GLib.FileStream
#endif
file);
}
