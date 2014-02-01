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
bool suppress_quality_diffs = false;
const string URL = "http://neufeldserver.uwaterloo.ca/~apmasell/mcbath-small_%d.fastq.bz2";

const OptionEntry[] options = {
	{ "forward", 'f', 0, OptionArg.FILENAME, ref forward_file, "Forward read FASTQ file.", "forward.fastq.bz2" },
	{ "reverse", 'r', 0, OptionArg.FILENAME, ref reverse_file, "Reverse read FASTQ file.", "reverse.fastq.bz2" },
	{ "suppress-quality", 'Q', 0, OptionArg.NONE, ref suppress_quality_diffs, "Ignore differences in quality scores of output bases.", null },
	{ "web", 'W', 0, OptionArg.NONE, ref web, "Get files from the web.", null },
	{ null }
};

public extern Panda.Assemble? create_assembler_control (Panda.LogProxy logger);
public extern Panda.Assemble? create_assembler_experiment (Panda.LogProxy logger);

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

	var control = create_assembler_control (logger);
	var experimental = create_assembler_experiment (logger);
	if (control == null || experimental == null) {
		return 1;
	}

	return Panda.diff (reader, control, experimental, suppress_quality_diffs) ? 2 : 0;
}
