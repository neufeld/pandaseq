/* PANDAseq-Diff -- Check differences between two versions of PANDAseq.
     Copyright (C) 2013  Andre Masella

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

bool panda_diff(
	PandaNextSeq reader,
	void *reader_data,
	PandaAssemble control,
	void *control_data,
	PandaAssemble experiment,
	void *experiment_data,
	bool suppress_quality_diffs) {

	size_t length_diffs = 0;
	size_t nt_diffs = 0;
	size_t diffs_better_score = 0;
	size_t diffs_worse_score = 0;
	size_t gained = 0;
	size_t lost = 0;
	size_t total = 0;

	const panda_qual *forward;
	size_t forward_length;
	const panda_qual *reverse;
	size_t reverse_length;
	panda_seq_identifier id;

	while (reader(&id, &forward, &forward_length, &reverse, &reverse_length, reader_data)) {
		const panda_result_seq *control_result = control(control_data, &id, forward, forward_length, reverse, reverse_length);
		const panda_result_seq *experiment_result = experiment(experiment_data, &id, forward, forward_length, reverse, reverse_length);
		total++;
		if (control_result == NULL && experiment_result == NULL) {
			/* Both fail. That's a match. */
			continue;
		} else if (control_result == NULL || experiment_result == NULL) {
			panda_seqid_print(&id, stdout);
			if (control_result == NULL) {
				gained++;
				fprintf(stdout, " has been gained.\n");
			} else {
				lost++;
				fprintf(stdout, " has been lost.\n");
			}
			continue;
		}

		if (control_result->quality < experiment_result->quality) {
			diffs_better_score++;
		} else if (control_result->quality > experiment_result->quality) {
			diffs_worse_score++;
		}

		if (experiment_result->sequence_length != control_result->sequence_length) {
			length_diffs++;
			panda_seqid_print(&id, stdout);
			fprintf(stdout, " differ in length %zd → %zd.\n", control_result->sequence_length, experiment_result->sequence_length);
		} else {
			bool nt_diff = false;
			size_t it;
			for (it = 0; it < experiment_result->sequence_length; it++) {
				if (control_result->sequence[it].nt != experiment_result->sequence[it].nt) {
					panda_seqid_print(&id, stdout);
					fprintf(stdout, " differ at nucleotide %zd, %c → %c.\n", it, panda_nt_to_ascii(control_result->sequence[it].nt), panda_nt_to_ascii(experiment_result->sequence[it].nt));
					nt_diff = true;
				} else if (control_result->sequence[it].p != experiment_result->sequence[it].p && !suppress_quality_diffs) {
					panda_seqid_print(&id, stdout);
					fprintf(stdout, " differ at nucleotide %zd (%c), quality %f → %f.\n", it, (int) panda_nt_to_ascii(control_result->sequence[it].nt), exp(control_result->sequence[it].p), exp(experiment_result->sequence[it].p));
					nt_diff = true;
				}
			}
			if (nt_diff) {
				nt_diffs++;
			}
		}
	}
	fprintf(stdout, "%zd sequences compared.\n%zd scored better.\n%zd scored worse.\n%zd changed (%zd length changed, %zd sequence changed).\n%zd gained.\n%zd lost.\n", total, diffs_better_score, diffs_worse_score, nt_diffs + length_diffs, length_diffs, nt_diffs, gained, lost);
	return (total == 0 || length_diffs > 0 || nt_diffs > 0 || gained > 0 || lost > 0);
}
