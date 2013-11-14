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

const OptionEntry[] options = {
	{ "forward", 'f', 0, OptionArg.FILENAME, ref forward_file, "Forward read FASTQ file.", "forward.fastq.bz2" },
	{ "reverse", 'r', 0, OptionArg.FILENAME, ref reverse_file, "Reverse read FASTQ file.", "reverse.fastq.bz2" },
	{ null }
};

[CCode(array_length_type = "size_t")]
unowned Panda.qual[]? forward = null;
[CCode(array_length_type = "size_t")]
unowned Panda.qual[]? reverse = null;

public int main(string[] args) {
	try {
		var opt_context = new OptionContext("- PANDAseq Diff");
		opt_context.set_help_enabled(true);
		opt_context.add_main_entries(options, null);
		if (!opt_context.parse(ref args)) {
			stdout.printf("Problem parsing arguments.\n");
			return 1;
		}
	} catch (OptionError e) {
		stdout.printf("%s\n", e.message);
		stdout.printf("Run '%s --help' to see a full list of available command line options.\n", args[0]);
		return 1;
	}
	if (forward_file == null) {
		stdout.printf("You must supply a forward read file.\n");
		return 1;
	}
	if (reverse_file == null) {
		stdout.printf("You must supply a reverse read file.\n");
		return 1;
	}

	var logger = new Panda.LogProxy.stderr();
	var fastq_reader = Panda.open_bz2(forward_file, reverse_file, logger);
	if (fastq_reader == null) {
		stdout.printf("Could not open input sequences.\n");
		return 1;
	}

	var old_assembler = new Panda.Assembler(null, logger);
	var assemble = panda_assembler_new_from_file(logger);
	if (assemble == null) {
		return 1;
	}

	var diffs = 0;
	var gained = 0;
	var lost = 0;
	var total = 0;
	Panda.identifier id;
	stdout.printf("Old Assembler: %p\n", old_assembler);
	while (fastq_reader(out id, out forward, out reverse)) {
		total++;
		unowned Panda.result_seq? old_result = old_assembler.assemble(id, forward, reverse);
		unowned Panda.result_seq? new_result = assemble(id, forward, reverse);
		if (old_result == null && new_result == null) {
			/* Both fail. That's a match. */
		} else if (old_result == null || new_result == null) {
			var is_new = old_result == null;
			id.to_file(stdout);
			if (is_new) {
				gained++;
				stdout.printf(" has been gained.\n");
			} else {
				lost++;
				stdout.printf(" has been lost.\n");
			}
		} else if (new_result.sequence.length != old_result.sequence.length) {
			diffs++;
			id.to_file(stdout);
			stdout.printf(" differ in length %d → %d.\n", old_result.sequence.length, new_result.sequence.length);
		} else {
			bool nt_diff = false;
			for (var it = 0; it < new_result.sequence.length; it++) {
				if (old_result.sequence[it].nt != new_result.sequence[it].nt) {
					id.to_file(stdout);
					stdout.printf(" differ at nucleotide %d, %c → %c.\n", it, old_result.sequence[it].nt.to_ascii(), new_result.sequence[it].nt.to_ascii());
					nt_diff = true;
				} else if (old_result.sequence[it].p != new_result.sequence[it].p) {
					id.to_file(stdout);
					stdout.printf(" differ at nucleotide %d, quality %f → %f.\n", it, old_result.sequence[it].p, new_result.sequence[it].p);
					nt_diff = true;
				}
			}
			if (nt_diff)
				diffs++;
		}
	}
	stdout.printf("%d sequences compared.\n%d changed.\n%d gained.\n%d lost.\n", total, diffs, gained, lost);
	return (total == 0 || diffs > 0 || gained > 0 || lost > 0) ? 2 : 0;
}
