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
#include <sys/types.h>
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
	fputs(panda_code_str(code), file);
	if (id != NULL) {
		fputc('\t', file);
		panda_seqid_xprint(id, (PandaPrintf) fprintf, file);
	}
	if (message != NULL) {
		fputc('\t', file);
		fputs(message, file);
	}
	fputc('\n', file);
	if (code == PANDA_CODE_PHRED_OFFSET) {
		fprintf(file, "* * * * * Using the default PHRED+33 offset, but no sequences had quality data under PHRED+64.\n* * * * * This is probably not what you want. Consult the manual about the -6 option.\n");
	}
	return true;
}

const char const *
panda_code_str(
	PandaCode code) {
	switch (code) {
	case PANDA_CODE_BAD_NT:
		return "ERR\tBADNT";
		break;
	case PANDA_CODE_ID_PARSE_FAILURE:
		return "ERR\tBADID";
		break;
	case PANDA_CODE_MOD_INFO:
		return "INFO\tMOD";
		break;
	case PANDA_CODE_NO_DATA:
		return "ERR\tNODATA";
		break;
	case PANDA_CODE_NO_FILE:
		return "ERR\tNOFILE";
		break;
	case PANDA_CODE_NO_QUALITY_INFO:
		return "ERR\tNOQUAL";
		break;
	case PANDA_CODE_NOT_PAIRED:
		return "ERR\tNOTPAIRED";
		break;
	case PANDA_CODE_PARSE_FAILURE:
		return "ERR\tBADSEQ";
		break;
	case PANDA_CODE_PREMATURE_EOF:
		return "ERR\tEOF";
		break;
	case PANDA_CODE_REJECT_STAT:
		return "STAT";
		break;
	case PANDA_CODE_INSUFFICIENT_KMER_TABLE:
		return "ERR\tKLNG";
		break;
	case PANDA_CODE_FORWARD_KMER:
		return "DBG\tFMER";
		break;
	case PANDA_CODE_REVERSE_KMER:
		return "DBG\tRMER";
		break;
	case PANDA_CODE_LOST_KMER:
		return "DBG\tFML";
		break;
	case PANDA_CODE_OVERLAP_POSSIBILITY:
		return "INFO\tOLD";
		break;
	case PANDA_CODE_BEST_OVERLAP:
		return "INFO\tBESTOLP";
		break;
	case PANDA_CODE_NO_FORWARD_PRIMER:
		return "ERR\tNOFP";
		break;
	case PANDA_CODE_NO_REVERSE_PRIMER:
		return "ERR\tNORP";
		break;
	case PANDA_CODE_LOW_QUALITY_REJECT:
		return "ERR\tLOWQ";
		break;
	case PANDA_CODE_NEGATIVE_SEQUENCE_LENGTH:
		return "ERR\tNEGS";
		break;
	case PANDA_CODE_SEQUENCE_TOO_LONG:
		return "ERR\tOOM";
		break;
	case PANDA_CODE_BUILD_FORWARD:
	case PANDA_CODE_BUILD_REVERSE:
	case PANDA_CODE_BUILD_OVERLAP:
		return "INFO\tBUILD";
		break;
	case PANDA_CODE_RECONSTRUCTION_PARAM:
		return "INFO\tRECR";
		break;
	case PANDA_CODE_MISMATCHED_BASE:
		return "INFO\tMISM";
		break;
	case PANDA_CODE_PHRED_OFFSET:
		return "INFO\tPHRED OFFSET";
		break;
	default:
		return "ERR\tUNKNOWN ERROR";
		break;
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
		(void) fputc(33 - (int) (10 * log(1 - exp(sequence->sequence[it].p)) / ln_10), file);
	}
	(void) fputc('\n', file);
	return true;
}
