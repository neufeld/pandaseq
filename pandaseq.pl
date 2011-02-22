#!/usr/bin/perl

# Filename: Illumina_assembly.v.1.0.pl
# Version: 1.0
# Date: 23 September 2010
# Author: Gabriel Moreno-Hagelsieb (gmoreno@wlu.ca), Michael D.J. Lynch (mdjlynch@sciborg.uwaterloo.ca)
# Usage: perl Illumina_assembly.pl <gzip'd Illumina file (forward)> <gzip'd Illumina file (reverse)>
# Availability: Updates to this program can be downloaded directly from http://popolvuh.wlu.ca/applications

####################################################################################################
# Copyright (C) 2010 Gabriel Moreno-Hagelsieb, Wilfrid Laurier University, Waterloo, Ontario, Canada
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################

#### Adapters (multiplexing barcodes). Add any new adapters to this group
@library_adapters = qw(
CTCTCT
ATTGGC
TGGTCA
GATCTG
GCCTAA
ACATCG
CGTGAT
CACTGT
CAAGTG
CGTACT
GTAGCC
GACTGA
AAGCTA
CTGATC
TCAAGT
TACAAG
CCTTGA
ACGGTA
ACAACC
AGTTGG
TGAGGA
TCGCTT
GCTCAA
ACCTCA
CGTGAT
ACATCG
GCCTAA
TGGTCA
CACTGT
ATTGGC
GATCTG
TCAAGT
CTGATC
AAGCTA
GTAGCC
TACAAG
CGTACT
GACTGA
GCTCAA
TCGCTT
TGAGGA
ACAACC
ACCTCA
ACGGTA
AGTTGG
CTCTCT
CAAGTG
CCTTGA
                  );

$middle_adapter = "AAGTCG";

#### Initial parameters
$total_seq_length = 228;#sequence length of each read is 108, plus index (6) and middle adapter (6)
$cut_position = $total_seq_length / 2;
$min_length   = 12;
$uno_length   = length($library_adapters[0]);
$dos_length   = length($middle_adapter);
$last_point   = $cut_position - $min_length - $dos_length;
@points = ( 0 .. $last_point );
@reversed_points = reverse( @points );

$base_name = $ARGV[0];
$reverse_name = $ARGV[1];

#### Results directory
$results_dir = "EXTRACTS_" . $base_name;
mkdir("$results_dir") unless(-d "$results_dir");

#### load all of the sequences into memory for analysis
open($INFILE_FORWARD, "-|","bzip2 -qdc $base_name");
open($INFILE_REVERSE, "-|","bzip2 -qdc $reverse_name");

FILE_LINE:
while($line_forward = <$INFILE_FORWARD>) {
	# get forward sequence name and sequence of index
    $index_string = substr $line_forward,index($line_forward,"#")+1,$uno_length;
	$forward_name = substr $line_forward,0,index($line_forward,"#");
	# get forward sequence
	$seq_f = <$INFILE_FORWARD>;
	chomp($seq_f);
	# get reverse sequence name
	$line_reverse = <$INFILE_REVERSE>;
	$reverse_name = substr $line_reverse,0,index($line_reverse,"#");
	# get reverse sequence
	$seq_r = <$INFILE_REVERSE>;
	chomp($seq_r);
	if($forward_name ne $reverse_name){
		print "Forward and reverse sequence names do not match up. Ending program\n";
		exit;
	}
	$index_string_rc = reverse_complement($index_string);
	$seq = $index_string_rc . $seq_f . $middle_adapter . $seq_r;
	my $to_count = $seq;
    my $count_correct_alphabet = $to_count =~ tr/ACGTN/ACGTN/;
    if( length($seq) > $count_correct_alphabet ) {#keeps track of the sequences with the wrong alphabet (i.e., that Ns)
        $wrong_alphabet++;
        $wrong_alphabet{"$seq"}++;
print STDERR "$seq has wrong alphabet\n";
        next FILE_LINE;
    }
    # now correct length
    if( length($seq) == $total_seq_length ) {
        #print length($seq),"\t",$seq,"\n";
        $count_of{"$seq"}++;
    }
    else {
    	print "Wrong Size\n";
        $count_wrong_size++;
        $wrong_size{"$seq"}++;
    }
}
close($INFILE_FORWARD);

#### Collect and report errors in the infile
open($ERROR,"|-","bzip2 -9 > $results_dir/$base_name.problematic.bz2");
for my $wrong ( sort keys %wrong_alphabet ) {
    print {$ERROR} join("\t",$wrong,"WRONG_ALPHABET"
                            ,$wrong_alphabet{"$wrong"}),"\n";
}
for my $wrong ( sort keys %wrong_size ) {
    print {$ERROR} join("\t",$wrong,"WRONG_SIZE"
                            ,$wrong_size{"$wrong"}),"\n";
}

#### Cut and process sequences
#processes the sequences
open($VERIF,"|-","bzip2 -9 > $results_dir/$base_name.fail_verif.bz2");
open($GVERIF,"|-","bzip2 -9 > $results_dir/$base_name.fine_verif.bz2");
SEQ:
for my $seq ( sort keys(%count_of) ) {
    $total_processed += $count_of{"$seq"};
    my $first  = substr($seq,0,$cut_position);
    my $second = substr($seq,$cut_position);
    my $library_here = substr($first,0,$uno_length);#trim the index (or library adaptor) off
    my $middle_here   = substr($second,0,$dos_length);#trim the index of the middle off
    ### might be reversed
    my $library_option   = substr($second,0,$uno_length);
    my $middle_option    = substr($first,0,$dos_length);
    my $library    = "";
    my $initial    = "";
    my $complement = "";
    my $terminal   = "";
    if( my @begin = grep { /$library_here/ } @library_adapters ) {#if index is in library
        $total_workable += $count_of{"$seq"};
        $complement .= reverse_complement("$second");
        $initial    .= $first;
        $terminal   .= $second;
        if( $middle_here eq $middle_adapter ) {
            $library .= $begin[0];
        }
        else {
            $library .= $begin[0] . "a";
        }
    }
    elsif( my @begin = grep { /$library_option/ } @library_adapters ) {
        $total_workable += $count_of{"$seq"};
        $complement .= reverse_complement("$first");
        $initial    .= $second;
        $terminal   .= $first;
        if( $middle_option eq $middle_adapter ) {
            $library .= $begin[0];
        }
        else {
            $library .= $begin[0] . "a";
        }
    }
    else {
        print {$ERROR} join("\t",$seq,"DOES_NOT_MATCH_LIBRARIES"
                                ,$count_of{"$seq"},$library_here
                                    ,$library_option),"\n";
        $complement .= reverse_complement("$second");
        print {$VERIF} join("\n",$seq,$first,$complement,$second),"\n\n";
        next SEQ;
    }
    
    #### Now to work
    if( length($library) > 0 ) {
        if( my $piece = compound("$initial","$complement") ) {
            $seqs_in_library{"$library"} += $count_of{"$seq"};
            $previously_saw{"$piece"}++;
            $count_library{"$library"}{"$piece"} += $count_of{"$seq"};
            if( $count_library{"$library"}{"$piece"}
                    == $count_of{"$seq"} ) {
                push( @{ $library },$piece );
                print {$GVERIF} join("\n",$library,$seq,$initial
                                         ,$complement,$terminal
                                             ,$piece),"\n\n";
            }
        }
        else {
            print {$ERROR} join("\t",$seq,"DOES_NOT_COMPLEMENT"
                                    ,$count_of{"$seq"}),"\n";
            print {$VERIF}  join("\n",$seq,$initial,$complement,$terminal)
                ,"\n\n";
        }
    }
    else {
        print {$ERROR} join("\t",$seq,"NONEXISTENT_LIBRARY"
                                ,$count_of{"$seq"}),"\n";
    }
}
close($ERROR);

@all_libraries = sort keys %seqs_in_library;
open(my $TBL,"|-","bzip2 -9 > $results_dir/$base_name.16S.tbl.bz2");
print {$TBL} join("\t","Library","Count","ID","Sequence"
                      ,"ReverseComplement"),"\n";
for my $library ( @all_libraries ) {
    my @seqs = @{ $library };
    my $count = @seqs;
    my $n = 1;
    open(my $FOR,"|-","bzip2 -9 > $results_dir/$library.$base_name.16SF.bz2");
    open(my $REV,"|-","bzip2 -9 > $results_dir/$library.$base_name.16SR.bz2");
    open(my $FULL,"|-"
             ,"bzip2 -9 > $results_dir/$library.$base_name.16SRfull.bz2");
    for my $piece ( sort @seqs ) {
        my $zeroes = length($count) - length($n);
        my $id = $library . 0 x $zeroes . $n;
        $n++;
        my $info = "F";
        my $sub_count = $count_library{"$library"}{"$piece"};
        my $rinfo = "R";
        my $comp = reverse_complement("$piece");
        print {$TBL} join("\t",$library,$sub_count,$id,$piece,$comp),"\n";
        print {$FOR} ">" . $id . $info . "|COUNT:" . $sub_count . "\n"
                . $piece . "\n";
        print {$REV} ">" . $id . $rinfo . "|COUNT:" . $sub_count . "\n"
                . $comp . "\n";
        for my $repeat ( 1 .. $sub_count ) {
            my $sub_zeroes = length($sub_count) - length($repeat);
            my $sub_id = 0 x $sub_zeroes . $repeat;
            print {$FULL} ">" . $id . $rinfo . $sub_id . "\n"
                . $comp . "\n";
        }
    }
    close($FOR);
    close($REV);
    close($FULL);
}
close($TBL);

sub reverse_complement {
    my ($in_seq) = @_;
    my $opposite = reverse $in_seq;
    $opposite =~ tr/ACGT/TGCA/;
    return("$opposite");
}

sub compound {
    my ($uno,$dos) = @_;
    $uno =~ s/^\S{$uno_length}//;#strips the beginning (index)
    $dos =~ s/\S{$dos_length}$//;#strips the middle_adapter
    my $complete_seq = "";
    for my $point ( @reversed_points ) {
        if( $point == 0 ) {
            my $matcher = substr($dos,$point,$min_length);#$matcher holds a substring of the the reverse complement sequence from 0 to the minimum length (12)
			my($prev,$good,$post) = $uno =~ /^(\S*)($matcher)(\S*)$/;
            if( my($prev,$good,$post) = $uno =~ /^(\S*)($matcher)(\S*)$/ ) {
                my $rematcher = $good . $post;
                if( my($regood,$repost) = $dos =~ /^($rematcher)(\S*)$/ ) {
                    my $complete = $prev . $rematcher . $repost;
                    $complete_seq .= $complete;
                    last;
                }
			}
        }
        else {
            my $totalmatch = substr($dos,$point);
            if( my($good) = $uno =~ /^($totalmatch)/ ) {
                $complete_seq .= $totalmatch;
                last;
            }
        }
    }
    #### Now report
    if( length($complete_seq) > 0 ) {
        return($complete_seq);
    }
    else {
        return;
    }
}
