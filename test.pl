# test.pl - test Bio::MCPrimers package
$VERSION='2.0';

# Author: Steve Lenk, May 2006

# Copyright (C) Stephen G. Lenk 2006.
# Artistic License applied.

# Use:    perl test.pl
# Output: OK or Not OK
# Return: 0 = OK, 1 = Not OK

use strict;
use warnings;

# use HIcysS.fa to generate a new test file
my $test = 1;
print "Test $test - ";
system('perl -Ilib mcprimers.pl -vector pET-32a -clamp 3prime -shift 24 HIcysS.fa > generated_file.pr3');

my $good_file;
my $generated_file;

open $good_file,      '<HIcysS.pr3'         or die "Can\'t open HIcysS.pr3\n";
open $generated_file, '<generated_file.pr3' or die "Can\'t open generated_file.pr3\n";

my $good_line      = <$good_file>;
my $generated_line = <$generated_file>;

# OK so far
my $status = 1;

my %good;
my %generated;
my $line = 0;

# compare good and generated solutions
CHECK: while (defined $good_line) {

    chomp $good_line;
    chomp $generated_line;

    # retain right primers in hashes
    if ($good_line =~ /RIGHT PRIMER(.*)/) {
        $good{$1} = 1;
    }    
    if ($generated_line =~ /RIGHT PRIMER(.*)/) {
        $generated{$1} = 1;
    }
    
    # compared line-by-line
    $line = $line + 1;
    if ($good_line ne $generated_line) {
        $status = 0;
        print "Difference at line $line\n";
    }

    $good_line      = <$good_file>;
    $generated_line = <$generated_file>;
}

close $good_file;
close $generated_file;

# see if good generated any unique solutions
foreach (keys %good) {
    if (not defined $generated{$_}) {
        $status = 0;
        print "Unique good: $_\n";
    }
}

# see if generated generated any unique solutions
foreach (keys %generated) {
    if (not defined $good{$_}) {
        $status = 0;
        print "Unique generated: $_\n";
    }
}

# see if there are redundant solutions in good
my %regenerated1;
foreach (keys %good) {
    if (defined $regenerated1{$_}) {
        $status = 0;
        print "Regenerated in good: $_\n";
    }
    $regenerated1{$_} = 1;
}

# see if there are redundant solutions in generated
my %regenerated2;
foreach (keys %generated) {
    if (defined $regenerated2{$_}) {
        $status = 0;
        print "Regenerated in generated: $_\n";
    }
    $regenerated2{$_} = 1;
}

# clean up
unlink 'generated_file.pr3';

#simple presentation of results
if ($status == 1) {
    print "OK\n";
    exit 0;
}
else {
    print "Not OK\n";
    exit 1;
}
