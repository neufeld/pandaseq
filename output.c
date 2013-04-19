/* PANDAseq -- Assemble paired FASTQ Illumina reads and strip the region between amplification primers.
     Copyright (C) 2011-2012  Andre Masella

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */
#include "config.h"
#include <math.h>
#include "pandaseq.h"
#include "buffer.h"
#include "table.h"

bool
panda_logger_file(
	PandaCode code,
	panda_seq_identifier *id,
	const char *message,
	FILE *file) {
#if HAVE_PTHREAD
	fprintf(file, "%p\t", static_buffer());
#endif
	(void) fputs(panda_code_str(code), file);
	if (id != NULL) {
		(void) fputc('\t', file);
		panda_seqid_xprint(id, (PandaPrintf) fprintf, file);
	}
	if (message != NULL) {
		(void) fputc('\t', file);
		(void) fputs(message, file);
	}
	(void) fputc('\n', file);
	if (code == PANDA_CODE_ID_PARSE_FAILURE && message != NULL) {
		(void) fprintf(file, "* * * * * Something is wrong with this ID. If tags are absent, try passing the -B option.\n* * * * * Consult `pandaseq-checkid \"%s\"` to get an idea of the problem..\n", message);
	} else if (code == PANDA_CODE_PHRED_OFFSET) {
		(void) fprintf(file, "* * * * * Using the default PHRED+33 offset, but no sequences had quality data under PHRED+64.\n* * * * * This is probably not what you want. Consult the manual about the -6 option.\n");
	} else if (code == PANDA_CODE_READ_TOO_LONG) {
		(void) fprintf(file, "* * * * * The input reads are longer than this version of PANDAseq can handle. Currently %zd nucleotides.\n", PANDA_MAX_LEN);
	}
	return true;
}

const char *const
panda_code_str(
	PandaCode code) {
	switch (code) {
	case PANDA_CODE_BAD_NT:
		return "ERR\tBADNT";
	case PANDA_CODE_ID_PARSE_FAILURE:
		return "ERR\tBADID";
	case PANDA_CODE_MOD_INFO:
		return "INFO\tMOD";
	case PANDA_CODE_NO_DATA:
		return "ERR\tNODATA";
	case PANDA_CODE_NO_FILE:
		return "ERR\tNOFILE";
	case PANDA_CODE_NO_QUALITY_INFO:
		return "ERR\tNOQUAL";
	case PANDA_CODE_NOT_PAIRED:
		return "ERR\tNOTPAIRED";
	case PANDA_CODE_PARSE_FAILURE:
		return "ERR\tBADSEQ";
	case PANDA_CODE_READ_TOO_LONG:
		return "ERR\tREADLEN";
	case PANDA_CODE_PREMATURE_EOF:
		return "ERR\tEOF";
	case PANDA_CODE_REJECT_STAT:
		return "STAT";
	case PANDA_CODE_INSUFFICIENT_KMER_TABLE:
		return "ERR\tKLNG";
	case PANDA_CODE_FORWARD_KMER:
		return "DBG\tFMER";
	case PANDA_CODE_REVERSE_KMER:
		return "DBG\tRMER";
	case PANDA_CODE_LOST_KMER:
		return "DBG\tFML";
	case PANDA_CODE_OVERLAP_POSSIBILITY:
		return "INFO\tOLD";
	case PANDA_CODE_BEST_OVERLAP:
		return "INFO\tBESTOLP";
	case PANDA_CODE_NO_FORWARD_PRIMER:
		return "ERR\tNOFP";
	case PANDA_CODE_NO_REVERSE_PRIMER:
		return "ERR\tNORP";
	case PANDA_CODE_LOW_QUALITY_REJECT:
		return "ERR\tLOWQ";
	case PANDA_CODE_NEGATIVE_SEQUENCE_LENGTH:
		return "ERR\tNEGS";
	case PANDA_CODE_SEQUENCE_TOO_LONG:
		return "ERR\tOOM";
	case PANDA_CODE_BUILD_FORWARD:
	case PANDA_CODE_BUILD_REVERSE:
	case PANDA_CODE_BUILD_OVERLAP:
		return "INFO\tBUILD";
	case PANDA_CODE_RECONSTRUCTION_PARAM:
		return "INFO\tRECR";
	case PANDA_CODE_MISMATCHED_BASE:
		return "INFO\tMISM";
	case PANDA_CODE_PHRED_OFFSET:
		return "INFO\tPHRED OFFSET";
	default:
		return "ERR\tUNKNOWN ERROR";
	}
}

bool
panda_output_fasta(
	const panda_result_seq *sequence,
	FILE *file) {
	size_t it;
	(void) fputc('>', file);
	panda_seqid_print(&sequence->name, file);
	(void) fputc('\n', file);
	for (it = 0; it < sequence->sequence_length; it++) {
		(void) fputc(panda_nt_to_ascii(sequence->sequence[it].nt), file);
	}
	(void) fputc('\n', file);
	return true;
}

bool
panda_output_fastq(
	const panda_result_seq *sequence,
	FILE *file) {
	size_t it;
	(void) fputc('@', file);
	panda_seqid_print(&sequence->name, file);
	(void) fputc('\n', file);
	for (it = 0; it < sequence->sequence_length; it++) {
		(void) fputc(panda_nt_to_ascii(sequence->sequence[it].nt), file);
	}
	fprintf(file, "\n+\n");
	for (it = 0; it < sequence->sequence_length; it++) {
		(void) fputc(33 - (int) (10 * log10(1 - exp(sequence->sequence[it].p)) - 0.5), file);
	}
	(void) fputc('\n', file);
	return true;
}

void
panda_output_fail(
	PandaAssembler assembler,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	FILE *file) {
	size_t it;
	(void) fputc('>', file);
	panda_seqid_print(id, file);
	(void) fputc('\n', file);
	for (it = 0; it < forward_length; it++) {
		(void) fputc(panda_nt_to_ascii(forward[it].nt), file);
	}
	(void) fputc('-', file);
	for (it = reverse_length; it > 0; it--) {
		(void) fputc(panda_nt_to_ascii(reverse[it - 1].nt), file);
	}
	(void) fputc('\n', file);
}
