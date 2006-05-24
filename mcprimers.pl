#!/usr/bin/perl
$VERSION = '2.0';

# mcprimers.pl - Designs molecular cloning PCR primers.

# Author:    Stephen G. Lenk, November 2005, May 2006.
# Copyright: Stephen G. Lenk (C) 2005, 2006. All rights reserved. 

# This program is free software; you can redistribute it and/or  
# modify it under the same terms as Perl itself.
# Licenced under the Artistic Licence.

#########################################################################
# This software comes with no guarantee of usefulness. 
# Use at your own risk. Check any solutions you obtain.
# Stephen G. Lenk assumes no responsibility for the use of this software.
#########################################################################

# Limitations: Does not account for redundancy codes in FASTA files.

# Note:        These runs use intermediate files. This means that
#              CGI or other use that can overwrite these files is
#              a bad idea. 

# see Bio::MCPrimers for POD

# Use with V2.0 of BIO::MCPrimers.pm

use strict;
use warnings;

use Bio::MCPrimers;


my $input_fh;
my $use_msg = qq/
MCPrimers generates molecular cloning PCR primers.
Copyright (C) 2005, 2006 Stephen G. Lenk.

Use:     perl mcprimers.pl [options] -vector name dna.fasta > result.pr3
Options: [-h] [-shift search_shift] [-clamp (both | 3prime)]

-h        help
-vector   cloning vector (i.e. pET-32a)
-shift    integer value to shift search of left RE match into gene
          use -shift 24 (or other value) if initial search fails
-clamp    GC clamp
          'both' will clamp last NT at both ends to G or C (default)
          '3prime' will clamp last NT on only 3' ends to G or C

Use at your risk. Check any solutions you obtain.
          
/;


my %flag;                       # hash for all flags
$flag{clamp}        = 'both';
$flag{search_shift} = 0;

my $flag;                       # flag currently being checked
my $vector_name;                # name of cloning vector

# command line options
while (@ARGV > 0) { 

    $flag = shift @ARGV;
    my $found = 0;

    if ($flag eq '-h') {
    
        # help
        print STDERR $use_msg;
        exit(1);
    }     
    elsif ($flag eq '-vector') {
        
        # select cloning vector
        if (@ARGV) {
            $vector_name = shift @ARGV;
            if (defined $vector_name and $vector_name ne '') {
                $found = 1;
            }
        }
    }   
    elsif ($flag eq '-shift') {
    
        # search shift
        if (@ARGV) {
            $flag{search_shift} = shift @ARGV;
            if ($flag{search_shift} =~ /^(\d)+$/) {
                $found = 1;
            } else {
                die "\nError: -shift $flag{search_shift} not recognised\n\n";
            }
        } else {
            die "\nError: -s needs an argument\n\n";
        }
    }
    elsif ($flag eq '-clamp') {

        # GC clamping
        if (@ARGV) {
            $flag{clamp} = shift @ARGV;
            if ($flag{clamp} eq 'both' or $flag{clamp} eq '3prime') {
                $found = 1;
            } else {
                die "\nError: -clamp $flag{clamp} not recognised\n\n";
            }
        } else {
            die "\n-Error: clamp needs an argument\n\n";
        }
    }
    
    # unknown flag
    if (not $found and @ARGV) {
        die "\nError: Option $flag not recognised\n\n";
    }
    
    # Should be a FASTA file name left at the end
    if ($found and @ARGV == 0) {
        die "\nError: Specify a FASTA file for input\n\n";  
    }
}


# Check that cloning vector has been defined
if (not defined $vector_name or $vector_name eq '') {
    die "\n\nError: Vector name not defined\n\n";
}


# open input file
my $infile = $flag;
if (not defined $infile) { 
    print STDERR $use_msg;
    exit(0);
}
open $input_fh, "<$infile" or die "\n\nError: Can\'t use $infile for input\n\n";


my $line;               # read lines
my $gene = '';          # then load ONLY nucleotides into $gene

# read fasta file

# toss > annotation
$line = '>';
while (substr($line,0,1) eq '>') {
    $line = <$input_fh>;   
}

# gete sequence data
while (defined $line) {

   chomp $line;
   $gene .= $line;
   $line = <$input_fh>;
}


my @re;           # restriction enzymes pattern
my %re_name;      # names of restriction enzymes
my @ecut_loc;     # cut location in enzyme
my @vcut_loc;     # cut location in vector

use Bio::Data::Plasmid::CloningVector;  # details of the plasmid used as a vector
my $status = Bio::Data::Plasmid::CloningVector::cloning_vector_data($vector_name,
                                                                    \@re, 
                                                                    \%re_name, 
                                                                    \@ecut_loc, 
                                                                    \@vcut_loc);            
if ($status == 0) {
    die "\n\nError: Data not found for cloning vector $vector_name\n\n";
}

my $version;      # version of MCPrimers returned from MCPrimers

# invoke solver
my $answer_ar;
$answer_ar = Bio::MCPrimers::find_mc_primers($gene, 
                                             \%flag, 
                                             \$version, 
                                             \@ecut_loc,
                                             \@vcut_loc,
                                             @re);


# Copyright notices, Primer3 here as MCPrimers uses Primer3
my $copr = 
qq/
|------------------------------------------------------------------|
| MCPrimers V$version                                                   |    
| Copyright (c) 2005,2006 Stephen G. Lenk. All rights reserved     |
| Primer3 Copyright (c) 1996,1997,1998,1999,2000,2001,2004         |  
| Whitehead Institute for Biomedical Research. All rights reserved |
|------------------------------------------------------------------|

Cloning vector    = $vector_name
Clamp flag        = $flag{clamp}
Shift flag        = $flag{search_shift}
Original sequence =
/;

print $copr;

# give them original gene sequence for reference
while ($gene =~ /(.{1,60})/g) {
    print "$1\n"; 
}
print "\n\n";

if (not defined $answer_ar or @{$answer_ar} == 0) { 

    # No solution found
    print "\nNo solution found\n\n";

} else {

    # Dump the solutions
    my $count = 0;
    foreach my $answer_hr (@{$answer_ar}) {

        # handle count
        $count += 1;
        my $count_text = 
            "=========\n=========  Solution # $count\n=========";
        
        # compose result text
        my $result = qq/
Start codon at  $answer_hr->{start}
Stop codon  at  $answer_hr->{stop}
Left RE site  = $re_name{$answer_hr->{left_re}} ($answer_hr->{left_re})
Right RE site = $re_name{$answer_hr->{right_re}} ($answer_hr->{right_re})

Primer3 analysis of PCR primers designed by MCPrimers:

$answer_hr->{primer3}/;

        print "$count_text\n$result\n\n";
    }
}


exit(0);
