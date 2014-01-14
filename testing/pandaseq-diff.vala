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
string? forward_file = null;
string? reverse_file = null;
bool web = false;
bool supress_quality_diffs = false;
const string URL = "http://neufeldserver.uwaterloo.ca/~apmasell/mcbath-small_%d.fastq.bz2";

const OptionEntry[] options = {
	{ "forward", 'f', 0, OptionArg.FILENAME, ref forward_file, "Forward read FASTQ file.", "forward.fastq.bz2" },
	{ "reverse", 'r', 0, OptionArg.FILENAME, ref reverse_file, "Reverse read FASTQ file.", "reverse.fastq.bz2" },
	{ "suppress-quality", 'Q', 0, OptionArg.NONE, ref supress_quality_diffs, "Ignore differences in quality scores of output bases.", null },
	{ "web", 'W', 0, OptionArg.NONE, ref web, "Get files from the web.", null },
	{ null }
};

[CCode (array_length_type = "size_t")]
unowned Panda.qual[]? forward = null;
[CCode (array_length_type = "size_t")]
unowned Panda.qual[]? reverse = null;

public int main (string[] args) {
	try {
		var opt_context = new OptionContext ("- PANDAseq Diff");
		opt_context.set_help_enabled (true);
		opt_context.add_main_entries (options, null);
		if (!opt_context.parse (ref args)) {
			stdout.printf ("Problem parsing arguments.\n");
			return 1;
		}
	} catch (OptionError e) {
		stdout.printf ("%s\n", e.message);
		stdout.printf ("Run '%s --help' to see a full list of available command line options.\n", args[0]);
		return 1;
	}
	var logger = new Panda.LogProxy (new Panda.Writer.null ());
	Panda.NextSeq reader;
	if (web) {
		var forward = Panda.bz_decompress (Panda.open_url (URL.printf (1), logger));
		var reverse = Panda.bz_decompress (Panda.open_url (URL.printf (2), logger));
		if (forward == null || reverse == null) {
			return 1;
		}
		reader = Panda.create_fastq_reader ((owned) forward, (owned) reverse, logger);
	} else {
		if (forward_file == null) {
			stdout.printf ("You must supply a forward read file.\n");
			return 1;
		}
		if (reverse_file == null) {
			stdout.printf ("You must supply a reverse read file.\n");
			return 1;
		}

		reader = Panda.open_fastq (forward_file, reverse_file, logger);
	}
	if (reader == null) {
		stdout.printf ("Could not open input sequences.\n");
		return 1;
	}

	var old_assembler = new Panda.Assembler (null, logger);
	var assemble = create_assembler (logger);
	if (assemble == null) {
		return 1;
	}

	var length_diffs = 0;
	var nt_diffs = 0;
	var diffs_better_score = 0;
	var diffs_worse_score = 0;
	var gained = 0;
	var lost = 0;
	var total = 0;
	Panda.identifier id;
	stdout.printf ("Old Assembler: %p\n", old_assembler);
	while (reader (out id, out forward, out reverse)) {
		total++;
		unowned Panda.result_seq? old_result = old_assembler.assemble (id, forward, reverse);
		unowned Panda.result_seq? new_result = assemble (id, forward, reverse);
		if (old_result == null && new_result == null) {
			/* Both fail. That's a match. */
			continue;
		} else if (old_result == null || new_result == null) {
			var is_new = old_result == null;
			id.to_file (stdout);
			if (is_new) {
				gained++;
				stdout.printf (" has been gained.\n");
			} else {
				lost++;
				stdout.printf (" has been lost.\n");
			}
			continue;
		}

		if (old_result.quality < new_result.quality) {
			diffs_better_score++;
		} else if (old_result.quality > new_result.quality) {
			diffs_worse_score++;
		}

		if (new_result.sequence.length != old_result.sequence.length) {
			length_diffs++;
			id.to_file (stdout);
			stdout.printf (" differ in length %d → %d.\n", old_result.sequence.length, new_result.sequence.length);
		} else {
			bool nt_diff = false;
			for (var it = 0; it < new_result.sequence.length; it++) {
				if (old_result.sequence[it].nt != new_result.sequence[it].nt) {
					id.to_file (stdout);
					stdout.printf (" differ at nucleotide %d, %c → %c.\n", it, old_result.sequence[it].nt.to_ascii (), new_result.sequence[it].nt.to_ascii ());
					nt_diff = true;
				} else if (old_result.sequence[it].log_probability != new_result.sequence[it].log_probability && !supress_quality_diffs) {
					id.to_file (stdout);
					stdout.printf (" differ at nucleotide %d (%c), quality %f → %f.\n", it, old_result.sequence[it].nt.to_ascii (), old_result.sequence[it].probability, new_result.sequence[it].probability);
					nt_diff = true;
				}
			}
			if (nt_diff) {
				nt_diffs++;
			}
		}
	}
	stdout.printf ("%d sequences compared.\n%d scored better\n%d scored worse.\n%d changed (%d length changed, %d sequence changed).\n%d gained.\n%d lost.\n", total, diffs_better_score, diffs_worse_score, nt_diffs + length_diffs, length_diffs, nt_diffs, gained, lost);
	return (total == 0 || length_diffs > 0 || nt_diffs > 0 || gained > 0 || lost > 0) ? 2 : 0;
}
