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

const char *panda_code_str(
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

bool panda_output_fasta(
	const panda_result_seq *sequence,
	PandaWriter writer) {
	size_t it;
	panda_writer_append_c(writer, '>');
	panda_writer_append_id(writer, &sequence->name);
	panda_writer_append(writer, ";%f", exp(sequence->quality));

	panda_writer_append_c(writer, '\n');
	for (it = 0; it < sequence->sequence_length; it++) {
		panda_writer_append_c(writer, panda_nt_to_ascii(sequence->sequence[it].nt));
	}
	panda_writer_append_c(writer, '\n');
	panda_writer_commit(writer);
	return true;
}

bool panda_output_fastq(
	const panda_result_seq *sequence,
	PandaWriter writer) {
	size_t it;
	panda_writer_append_c(writer, '@');
	panda_writer_append_id(writer, &sequence->name);
	panda_writer_append(writer, ";%f", exp(sequence->quality));
	panda_writer_append_c(writer, '\n');
	for (it = 0; it < sequence->sequence_length; it++) {
		panda_writer_append_c(writer, panda_nt_to_ascii(sequence->sequence[it].nt));
	}
	panda_writer_append(writer, "\n+\n");
	for (it = 0; it < sequence->sequence_length; it++) {
		panda_writer_append_c(writer, 33 + panda_result_phred(&sequence->sequence[it]));
	}
	panda_writer_append_c(writer, '\n');
	panda_writer_commit(writer);
	return true;
}

void panda_output_fail(
	PandaAssembler assembler,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	PandaWriter writer) {
	size_t it;
	(void) assembler;
	panda_writer_append_c(writer, '>');
	panda_writer_append_id(writer, id);
	panda_writer_append_c(writer, '\n');
	for (it = 0; it < forward_length; it++) {
		panda_writer_append_c(writer, panda_nt_to_ascii(forward[it].nt));
	}
	panda_writer_append_c(writer, '-');
	for (it = reverse_length; it > 0; it--) {
		panda_writer_append_c(writer, panda_nt_to_ascii(reverse[it - 1].nt));
	}
	panda_writer_append_c(writer, '\n');
	panda_writer_commit(writer);
}

void panda_output_fail_qual(
	PandaAssembler assembler,
	const panda_seq_identifier *id,
	const panda_qual *forward,
	size_t forward_length,
	const panda_qual *reverse,
	size_t reverse_length,
	PandaWriter writer) {
	size_t it;
	(void) assembler;
	panda_writer_append_c(writer, '@');
	panda_writer_append_id(writer, id);
	panda_writer_append_c(writer, '\n');
	for (it = 0; it < forward_length; it++) {
		panda_writer_append_c(writer, panda_nt_to_ascii(forward[it].nt));
	}
	panda_writer_append_c(writer, '-');
	for (it = reverse_length; it > 0; it--) {
		panda_writer_append_c(writer, panda_nt_to_ascii(reverse[it - 1].nt));
	}
	panda_writer_append(writer, "\n+\n");
	for (it = 0; it < forward_length; it++) {
		panda_writer_append_c(writer, 33 + forward[it].qual);
	}
	panda_writer_append_c(writer, '!');
	for (it = reverse_length; it > 0; it--) {
		panda_writer_append_c(writer, 33 + reverse[it - 1].qual);
	}
	panda_writer_append_c(writer, '\n');
	panda_writer_commit(writer);
}
