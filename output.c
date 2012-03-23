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
#include<math.h>
#include<sys/types.h>
#include"pandaseq.h"
#include"config.h"
bool panda_logger_file(FILE *file, PandaCode code, ...) {
	va_list args;
	bool ret;
	va_start(args, code);
	ret = panda_logger_v((PandaPrintf)fprintf, file, code, args);
	va_end(args);
	return ret;
}

bool panda_logger_v(PandaPrintf xprintf, void *x, PandaCode code, va_list va)
{
	double d;
	PandaModule m;
	int i;
	int i2;
	int i3;
	panda_qual *q1;
	panda_qual *q2;
	panda_result *r;
	size_t s1;
	size_t s2;
	size_t s3;
	size_t s4;
	switch (code) {
		case PANDA_CODE_API_VERSION:
			xprintf(x, "INFO\tAPI\t%d\n", PANDA_API);
			break;
		case PANDA_CODE_BAD_NT:
			xprintf(x, "ERR\tBADNT\t");
			panda_seqid_xprint(va_arg(va, panda_seq_identifier*), xprintf, x);
			i = va_arg(va, int);
			xprintf(x, "\t%c@%d\n", i, va_arg(va, int));
			break;
		case PANDA_CODE_ID_PARSE_FAILURE:
			xprintf(x, "ERR\tBADID\t%s\n", va_arg(va, const char*));
			break;
		case PANDA_CODE_MOD_INFO:
			m = va_arg(va, PandaModule);
			xprintf(x, "INFO\tMOD\t%s(%s:%d)\t%s\n", panda_module_get_name(m), panda_module_get_version(m), panda_module_get_api(m), panda_module_get_args(m));
			break;
		case PANDA_CODE_NO_DATA:
			xprintf(x, "ERR\tNODATA\t");
			panda_seqid_xprint(va_arg(va, panda_seq_identifier*), xprintf, x);
			xprintf(x, "\n");
			break;
		case PANDA_CODE_NO_FILE:
			xprintf(x, "ERR\tNOFILE\t%s\n", va_arg(va, const char*));
			break;
		case PANDA_CODE_NO_QUALITY_INFO:
			xprintf(x, "ERR\tNOQUAL\t");
			panda_seqid_xprint(va_arg(va, panda_seq_identifier*), xprintf, x);
			xprintf(x, "\n");
			break;
		case PANDA_CODE_NOT_PAIRED:
			xprintf(x, "ERR\tNOTPAIRED\t");
			panda_seqid_xprint(va_arg(va, panda_seq_identifier*), xprintf, x);
			xprintf(x, "\t");
			panda_seqid_xprint(va_arg(va, panda_seq_identifier*), xprintf, x);
			xprintf(x, "\n");
			break;
		case PANDA_CODE_PARSE_FAILURE:
			xprintf(x, "ERR\tBADSEQ\t");
			panda_seqid_xprint(va_arg(va, panda_seq_identifier*), xprintf, x);
			xprintf(x, "\n");
			break;
		case PANDA_CODE_PREMATURE_EOF:
			xprintf(x, "ERR\tEOF\t");
			panda_seqid_xprint(va_arg(va, panda_seq_identifier*), xprintf, x);
			xprintf(x, "\n");
			break;
		case PANDA_CODE_REJECT_STAT:
			m = va_arg(va, PandaModule);
			xprintf(x, "INFO\t%s\t%ld\n", panda_module_get_name(m), va_arg(va, long));
			break;
		case PANDA_CODE_INSUFFICIENT_KMER_TABLE:
			xprintf(x, "ERR\tKLNG\t");
			panda_seqid_xprint(va_arg(va, panda_seq_identifier*), xprintf, x);
			xprintf(x, "\n");
			break;
		case PANDA_CODE_FORWARD_KMER:
			(void)va_arg(va, panda_seq_identifier*);
			i = va_arg(va, unsigned int);
			xprintf(x, "DBG\tFMER\t%x @ %d\n", i, va_arg(va, ssize_t));
			break;
		case PANDA_CODE_REVERSE_KMER:
			(void)va_arg(va, panda_seq_identifier*);
			i = va_arg(va, unsigned int);
			xprintf(x, "DBG\tRMER\t%x @ %d\n", i, va_arg(va, ssize_t));
			break;
		case PANDA_CODE_LOST_KMER:
			(void)va_arg(va, panda_seq_identifier*);
			i = va_arg(va, unsigned int);
			xprintf(x, "DBG\tFML\t%x @ %d\n", i, va_arg(va, ssize_t));
			break;
		case PANDA_CODE_OVERLAP_POSSIBILITY:
			(void)va_arg(va, panda_seq_identifier*);
			s1 = va_arg(va, size_t);
			s2 = va_arg(va, size_t);
			s3 = va_arg(va, size_t);
			s4 = va_arg(va, size_t);
			xprintf(x, "INFO\tOLD\t%d mat = %d, mismat = %d, unk = %d, prob = %f\n", s1, s2, s3, s4, va_arg(va, double));
			break;
		case PANDA_CODE_BEST_OVERLAP:
			(void)va_arg(va, panda_seq_identifier*);
			xprintf(x,  "INFO\tBESTOLP\t%d\n", va_arg(va, int));
			break;
		case PANDA_CODE_NO_FORWARD_PRIMER:
			xprintf(x, "ERR\tNOFP\t\n");
			break;
		case PANDA_CODE_NO_REVERSE_PRIMER:
			xprintf(x, "ERR\tNORP\t\n");
			break;
		case PANDA_CODE_LOW_QUALITY_REJECT:
			(void)va_arg(va, panda_seq_identifier*);
			d = va_arg(va, double);
			xprintf(x, "ERR\tLOWQ\t%f %f\n", d, va_arg(va, double));
			break;
		case PANDA_CODE_NEGATIVE_SEQUENCE_LENGTH:
			xprintf(x, "ERR\tNEGS\n");
			break;
		case PANDA_CODE_SEQUENCE_TOO_LONG:
			xprintf(x, "ERR\tOOM\n");
			break;
		case PANDA_CODE_BUILD_FORWARD:
		case PANDA_CODE_BUILD_REVERSE:
			(void)va_arg(va, panda_seq_identifier*);
			i = va_arg(va, int);
			i2 = va_arg(va, int);
			r = va_arg(va, panda_result*);
			xprintf(x, "INFO\tBUILD\tr[%d] = %c[%d] | %e\n", i, (int)panda_nt_to_ascii(r->nt), i2, exp(r->p));
			break;
		case PANDA_CODE_BUILD_OVERLAP:
			(void)va_arg(va, panda_seq_identifier*);
			i = va_arg(va, int);
			i2 = va_arg(va, int);
			i3 = va_arg(va, int);
			r = va_arg(va, panda_result*);
			q1 = va_arg(va, panda_qual*);
			q2 = va_arg(va, panda_qual*);
			
			xprintf(x, "INFO\tBUILD\tr[%d] = %c | %e given %c[%d] | %e and %c[%d] | %e\n", i, (int)panda_nt_to_ascii(r->nt), exp(r->p), (int)panda_nt_to_ascii(q1->nt), i2, panda_quality_probability(q1),
			(int)panda_nt_to_ascii(q2->nt), i3, panda_quality_probability(q1));
			break;
		case PANDA_CODE_RECONSTRUCTION_PARAM:
			(void)va_arg(va, panda_seq_identifier*);
			i = va_arg(va, int);
			i2 = va_arg(va, int);
			xprintf(x, "INFO\tRECR\toverlap = %d df = %d dr = %d\n", i, i2, va_arg(va, int));
			break;
		case PANDA_CODE_MISMATCHED_BASE:
			(void)va_arg(va, int);
			(void)va_arg(va, int);
			q1 = va_arg(va, panda_qual*);
			q2 = va_arg(va, panda_qual*);
			xprintf(x, "INFO\tMISM\t%d %d\n", (int)q1->qual, (int)q2->qual);
			break;
		default:
			xprintf(x, "ERR\tUNKNOWN ERROR\t%x\n", code);
			break;
	}
	return true;
}

const char *panda_version(void) {
	return PACKAGE_STRING;
}

bool panda_output_fasta(const panda_result_seq *sequence, FILE *file) {
	size_t it;
	(void)fputc('>', file);
	panda_seqid_print(&sequence->name, file);
	(void)fputc('\n', file);
	for(it = 0; it < sequence->sequence_length; it++) {
		(void)fputc(panda_nt_to_ascii(sequence->sequence[it].nt), file);
	}
	(void)fputc('\n', file);
	return true;
}
bool panda_output_fastq(const panda_result_seq *sequence, FILE *file) {
	static double ln_10 = log(10);
	size_t it;
	(void)fputc('@', file);
	panda_seqid_print(&sequence->name, file);
	(void)fputc('\n', file);
	for(it = 0; it < sequence->sequence_length; it++) {
		(void)fputc(panda_nt_to_ascii(sequence->sequence[it].nt), file);
	}
	fprintf(file, "\n+\n");
	for(it = 0; it < sequence->sequence_length; it++) {
		(void)fputc(33 - (int)(10 * log(1 - exp(sequence->sequence[it].p)) / ln_10), file);
	}
	(void)fputc('\n', file);
	return true;
}

