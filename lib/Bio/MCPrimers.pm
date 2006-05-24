package Bio::MCPrimers;
our $VERSION = '1.5.1';

use strict;
use warnings;


# Bio::MCPrimers.pm - generates molecular cloning PCR primers

##########################################################################
# This software comes with no guarantee of usefulness.
# Use at your own risk. Check any solutions you obtain.
# Stephen G. Lenk assumes no responsibility for the use of this software.
# Licenced under the Artistic Licence.
##########################################################################

#######################################################
###### See bottom for pod and other information    ####
#######################################################

my $primer3_name;       # set depending upon OS
my $primer3_dir;        # from $ENV or set to default
my $primer3_exe;        # full executable name

BEGIN {

    # check PRIMER3DIR
    if (not (defined $ENV{PRIMER3DIR})) {
        $primer3_dir = '.'
    }   
    
    if ($^O =~ /^MSW/) { 
    
        # Microsoft
        $primer3_name = 'primer3.exe';
        $primer3_exe  = "$primer3_dir\\$primer3_name";
    } 
    else {         
    
        # non-Microsoft
        $primer3_name = 'primer3_core';
        $primer3_exe  = "$primer3_dir/$primer3_name";
    }

    # Is Primer3 executable available?
    if (not -e $primer3_exe) {    
        print STDERR "\n$primer3_exe not available\n\n";
        exit 1;
    }
}



# forward declarations for subs
sub find_mc_primers;
sub _solver;
sub _get_re_patterns;
sub _get_re_matches;
sub _start_codon;
sub _stop_codon;
sub _primers_ok;
sub _primer_pairs;
sub _create_primers;
sub _create_left_primers;
sub _create_right_primers;
sub _handle_reply;
sub _sanity_check_gene;
sub _primer_pairs;
sub _number_dots;
sub _generate_re_patterns;
sub _too_many_substitutions_in_a_row;
sub _check_3prime_gc;

my $CODON_SIZE = 3;   # likely to stay a constant

my %left_bad;         # for tracking bad primers during calculation
my %right_bad;        # and check right as well

my $fh_write;         # write to Primer3
my $fh_read;          # read from Primer3



#### front end solver sets up for OS needs ####
sub find_mc_primers {

    my ($gene,           # ATGC string
        $flag_hr,        # anonymous hash reference to flags
        $version_sr,     # version scalar reference    
        $ecut_loc_ar,    # enzyme cut location array reference
        $vcut_loc_ar,    # vector cut location array reference
        @re              # array of restriction enzyme strings
       ) = @_;
       
    ${$version_sr} = $VERSION;  # return version to caller

    my $search_shift = 0;  # permit the search for left primer
                           # to shift to the right
    my $clamp = 'both';    # GC clamp at both ends  

    # digest flags into local variables
    if (defined $flag_hr->{search_shift}) { 
        $search_shift = $flag_hr->{search_shift} 
    } 
    if (defined $flag_hr->{clamp}) {

        if ($flag_hr->{clamp} eq 'both' or 
            $flag_hr->{clamp} eq '3prime') {

            $clamp = $flag_hr->{clamp}
        }
    }
    
    # invoke solver
    my $answer = _solver($gene, 
                         $search_shift, 
                         $clamp, 
                         $ecut_loc_ar,
                         $vcut_loc_ar,
                         @re);

    # clean up intermediate files for Primer3    
    unlink 'p3in.txt';      
    unlink 'p3out.txt';
       
    return $answer;
}



### solver function ####
sub _solver {

    my ($gene,         # ATGC string
        $search_shift, # extra search into start of gene for left primer
        $clamp,        # both or 3prime
        $ecut_loc_ar,  # enzyme cut location array reference
        $vcut_loc_ar,  # vector cut location array reference
        @re            # array of RE strings
       ) = @_;

    my $MAX_CHANGES  = 3;       # changes to gene for primer
    my $solution_ar  = [];      # solution array reference

    my @re_matches;   # array of anonymous arrays of RE matches
    for (my $i = 0; $i < @re; $i++) { 
        $re_matches[$i] = [];
    }

    my $gene_start = _start_codon($gene); 
    my $gene_stop  = _stop_codon($gene, $gene_start); 

    if (not defined $gene_stop) {
        return $solution_ar;
    }

    # generate RE matches
    for (my $i=0; $i < @re; $i += 1) {

        my @patterns = _get_re_patterns($re[$i], $MAX_CHANGES); 
        my @matches  = _get_re_matches($gene, 
                                       $ecut_loc_ar->[$i],
                                       $vcut_loc_ar->[$i],
                                       @patterns);

        foreach my $match (@matches) { 
            push @{$re_matches[$i]}, $match;
        }
    }
    
    # keep pulling out a potential left primer until exhausted or solved
    while (@re > 1) {

        my $re_l         = shift @re;
        my $re_pos_l_ref = shift @re_matches;

        # position match list for head enzyme
        foreach my $re_pos_l (@{$re_pos_l_ref}) {

            # control where left primer is placed
            if ($re_pos_l > (10 + $gene_start + $search_shift)) {            
                next; 
            }
         
            # left primers
            my $modified_gene = $gene;
            substr($modified_gene, $re_pos_l, length $re_l, $re_l);
            my @left = _create_left_primers($modified_gene, 
                                            $re_pos_l, 
                                            $re_l,
                                            $clamp);
                                            
            # loop across rest of sites
            my $num_sites_left = @re;
            my $site = 0;
            
            # rest of enzymes left after head has been shifted off
            while ($site < $num_sites_left) {
       
                my $re_r = $re[$site];
                my $re_pos_r_ref = $re_matches[$site]; 
                my $right_hr = {};

                # rest of restriction sites going right in vector
                RIGHT: foreach my $re_pos_r (@{$re_pos_r_ref}) {
    
                    if ($re_pos_r < $gene_stop) {  
                        next RIGHT 
                    }
    
                    # right primers
                    my $more_modified_gene = $modified_gene;
                    substr($more_modified_gene, $re_pos_r, length $re_r, $re_r);     
             
                    my @right = _create_right_primers($more_modified_gene, 
                                                      $re_pos_r, 
                                                      $re_r,
                                                      $clamp);

                    # test new primers only, bypass primers already tested
                    foreach (@right) { 
                        if (defined $right_hr->{$_}) { next RIGHT } 
                    }

                    # generate primer pairs and have them checked
                    my $reply = _primer_pairs(\@left, 
                                              \@right, 
                                              $more_modified_gene, 
                                              $gene_start,
                                              $gene_stop,
                                              $search_shift);

                    # solution obtained if $reply is defined
                    if (defined $reply) { 

                        # see if solution is valid, if so then 
                        # append one solution to @{$solution_ar}
                        _handle_reply($reply,
                                      $more_modified_gene,
                                      $re_l,
                                      $re_r,
                                      $gene_start,
                                      $gene_stop,
                                      $solution_ar);
                    } 

                    # keep track of which right primers have been tested
                    foreach (@right) { $right_hr->{$_} = 1 }
                }
                
                # next RE site for possible insertion
                $site = $site + 1;  
            }
        } 
    }

    # all solutions
    return $solution_ar;
}



#### handle reply from solution checker ####
sub _handle_reply {

    my ($reply,                 # solution checker reply
        $more_modified_gene,    # ATGC modified for left and right RE
        $re_l,                  # left RE
        $re_r,                  # right RE
        $gene_start,            # location of START codon
        $gene_stop,             # location of STOP codon
        $solution_ar            # array reference to solution
       ) = @_;

    # check that RE sequences have only one match
    my $left_cnt  = 0;
    my $right_cnt = 0;
    while ($more_modified_gene =~ /$re_l/g) { 
        $left_cnt  += 1 
    }
    while ($more_modified_gene =~ /$re_r/g) { 
        $right_cnt += 1 
    }

    # one match for each?
    if ($left_cnt == 1 and $right_cnt == 1) { 

        # add anonymous hash to anonymous array
        push @{$solution_ar}, { primer3  => $reply,
                                left_re  => $re_l,
                                right_re => $re_r,
                                start    => $gene_start,
                                stop     => $gene_stop };
        return 1; # solution added
    }
    
    return 0; # solution not added
}



#### check primers with Primer3 ####
sub _primers_ok {

    my ($left,          # left primer string
        $right,         # right primer string
        $gene,          # ATGC string modified to match primers
        $gene_start,    # location of start codon
        $gene_stop,     # location of stop codon
        $search_shift   # move search zone into gene from left
       ) = @_;
 
    my $ok = 1;
    
    if (defined $right_bad{$right}) { 
        $ok = 0 
    }
    if (defined $left_bad{$left})   { 
        $ok = 0 
    }
    
    # outta here if one of the primers has been identified as bad
    if ($ok == 0) { 
        return undef 
    }
  
    # create Boulder file text for Primer3
    my $range = length $gene; 
    my $excluded_region_start = $gene_start + 30 + $search_shift;
    my $excluded_length =  $gene_stop - $excluded_region_start + 3;
    
    my @boulder = 
       ( "SEQUENCE=$gene",
         "PRIMER_PRODUCT_MAX_SIZE=$range",
         "PRIMER_PRODUCT_SIZE_RANGE=100-$range",
         "PRIMER_LEFT_INPUT=$left",
         "PRIMER_RIGHT_INPUT=$right",
         "EXCLUDED_REGION=$excluded_region_start,$excluded_length",
         "PRIMER_EXPLAIN_FLAG=1",
         "=" );
   
    # write intermediate file for Primer3
    open  $fh_write, ">p3in.txt";
    foreach (@boulder) { 
        print $fh_write "$_\n"; 
    }
    close $fh_write;
    
    # primer3 call done here
    my $status;    # for system call
    $status = system("$primer3_exe -format_output <p3in.txt >p3out.txt");
    if ($status != 0) {
        print STDERR "\nError: Primer3 error $status\n";
        exit 1;
    }
        
    my $p3_reply;
    
    # go through Primer3 output
    open $fh_read, "<p3out.txt";   
    PRIMER3_READ: while (<$fh_read>) { 

        my $line = $_;       
        $p3_reply .= $line;

        if ($line =~ /NO PRIMERS FOUND/) { 
        
            # solution fails primer3
            $ok = 0; 
        }        
         if ($line =~ /^primer3 release/) {

            # done with primer3 for this primer pair
            last PRIMER3_READ;
        }   
        if ($line =~ /PRIMER_ERROR/) {
        
            # Primer3 found an error
            print STDERR "\nError: Primer3 error: $line\n";
            exit 1;
        }        
        if ($line =~ /PROBLEM/) {

            # done with primer3 for this primer pair
            print STDERR "\n$line\n";
            exit 1;
        } 

        # check left and right
        if ($line =~ /^Left.*0$/)  { 
            $ok = 0; 
            $left_bad{$left} = 1 
        }
        if ($line =~ /^Right.*0$/) { 
            $ok = 0; 
            $right_bad{$right} = 1 
        }
    }
 
    close $fh_read;
    
    # whew! a solution
    if ($ok == 1) { 
        return $p3_reply; 
    }

    # no solution
    return undef;
}



#### create the primers ####
sub _create_primers {

    my ($re,            # restriction enzyme sequence
        $gene,          # ATGC sequence modified to match RE
        $primers_ref,   # reference to primers array
        $clamp          # both left right
       ) = @_;

    my @qs = ( ['', ''], ['?', ''], ['', '?'], ['?', '?'] );
    
    # padding outside of RE
    for my $pad (3 .. 12) {
    
        # ? marks for different types of matching
        foreach my $q (@qs) {
        
            # left and right '?' for regular expression matches
            my $l = $q->[0];
            my $r = $q->[1];
            
            # establish proper pattern
            my $pattern;

            if ($clamp eq 'both')  {            
                $pattern = "[GC](.{3,$pad}$l)$re(.{3,$pad}$r)[GC]" 
            }
            if ($clamp eq 'left')  { 
                $pattern = ".(.{3,$pad}$l)$re(.{3,$pad}$r)[GC]" 
            }
            if ($clamp eq 'right') { 
                $pattern = "[GC](.{3,$pad}$l)$re(.{3,$pad}$r)." 
            }

            # pattern matches
            while ($gene =~ /($pattern)/g) { 
    
                my $location = $-[0]; 
                my $primer   = $1;
                
                # limit primer sizes
                my $l = length $1;
                if ($l >= 18 and $l <= 24) {               
                    $primers_ref->{$primer} = $location; 
                }
            }
        }
    }

    return undef;
}



#### sanity check gene ####
#### see if stop codon has been 'wacked' ####
sub _sanity_check_gene {

    my ($gene,         # ATGC as primers will make
        $gene_start,   # start codon
        $gene_stop     # stop codon
       ) = @_;

    my $stop = _stop_codon($gene, $gene_start);
    
    if (not defined $stop) {     
        return 0 
    }
    
    if ($gene_stop == $stop) {
        return 1;

    } else {
        return 0; 
    }
}



#### generate primer pairs, then process them one by one ####
sub _primer_pairs {

    my ($left_primers_ref,   # array reference
        $right_primers_ref,  # array reference
        $gene,               # ATGC
        $gene_start,         # start codon location
        $gene_stop,          # stop codon location
        $search_shift        # shift search in from right
       ) = @_;

    if (@{$left_primers_ref} == 0 or @{$right_primers_ref} == 0) {

        # need both left and right or go home
        return undef 
    }

    # lefties
    foreach my $left_pr (@{$left_primers_ref}) { 

        # righties
        foreach my $right_pr (@{$right_primers_ref}) { 

            # sequence to be made OK
            if (_sanity_check_gene($gene, $gene_start, $gene_stop) == 1) { 
                
                # primers OK
                my $reply = _primers_ok($left_pr, 
                                        $right_pr, 
                                        $gene, 
                                        $gene_start, 
                                        $gene_stop, 
                                        $search_shift);

                if (defined $reply) { 

                    # reply is OK here
                    return $reply 
                }
            }
        }
    }

    return undef;
}



#### how many '.' in pattern ####
sub _number_dots {

    my (@chars    # array of characters in pattern
       ) = @_;
       
    my $num = 0;

    foreach (@chars) { 
        if ($_ eq '.') { 
            $num += 1 
        } 
    }

    return $num;
}



#### see if there are too many substitutions in a row being requested ####
sub _too_many_substitutions_in_a_row {

    my ($max_dots,    # maximum nuber of dots
        @chars        # array of characters to check
       ) = @_;

    my $n_in_a_row = 0;

    # count '.' in a row, reset where needed
    foreach my $c (@chars) {

        if ($c eq '.') { 
            $n_in_a_row += 1 
        } else { 
            $n_in_a_row = 0 
        }

        if ($n_in_a_row == $max_dots) {   
            return 1 
        }
    }

    return 0;
}



#### recursively generate patterns ####
sub _generate_re_patterns {

    my ($max_dots,     # maximum number of dots also limits recursion
        $pattern_hr,   # pattern hash reference
        @r             # incoming pattern to modify
       ) = @_;  

    my @s;   # next pattern

    # add to hash, keep track of only one time it's found !!!!!!
    $pattern_hr->{ join '', @r } = '';
   
    # limit to annealing capability of primer
    if (_number_dots(@r) == $max_dots)  {   
        return undef 
    }

    if (_too_many_substitutions_in_a_row($max_dots, @r) == 1) {     
        return undef 
    }

    # successively generate next group of patterns
    for my $i (1 .. @r) { 
    
        # already have a '.' here
        if ($r[$i-1] eq '.') {      
            next 
        }
        
        # empty @s, generate clean pattern array
        while (@s) {        
            pop @s 
        }
    
        # build next pattern
        for my $j (1 .. @r) {
        
            if ($i == $j) {             
                push @s, '.'                
            } else {            
                push @s, $r[$j-1] 
            }
        }
            
        # keep going
        _generate_re_patterns( $max_dots, $pattern_hr, @s ); 
    }

    # all patterns stored in hash
    my @k=keys %{$pattern_hr};     
    
    return @k;
}



#### regular expression patterns with '.' generated ####
sub _get_re_patterns {

    my ($re,         # restriction enzyme
        $max_dots    # maximum number of '.'
       ) = @_;     
    
    my @re = split '', $re; # characters in RE
    my %l;                  # hash function with list of generated RE

    # generate patterns for requested enzyme
    my @pats = _generate_re_patterns($max_dots, \%l, @re);

    # sort patterns here
    my @sorted;

    foreach my $n (0 .. $max_dots) {
        foreach my $p (@pats) { 
            if (_number_dots(split '', $p) == $n) { 
                push @sorted, $p 
            }
        }
    }

    return @sorted;
}

    

#### get matches for re pattern in gene ####
sub _get_re_matches {

    my ($gene,       # ATGC sequence
        $ecut_loc,   # enzyme cut location
        $vcut_loc,   # vector cut location
        @patterns    # array of RE patterns
       ) = @_;

    my @positions;
    my %used;

    # loop through patterns
    foreach my $p (@patterns) {

        # loop across gene
        while ($gene =~ /($p)/g) {
            
            # check for in-frame
            if ($vcut_loc == (($-[0] + $ecut_loc) % 3)) {

                # only use a location once
                if (not defined $used{$-[0]}) {
                    push @positions, $-[0]; 
                }
                $used{$-[0]} = 1;
            }
        } 
    }   

    return @positions;
}



#### create left primers ####
sub _create_left_primers {
    
    my ($modified_gene,   # ATGC modified for left RE
        $re_pos_l,        # position of left RE
        $re_l,            # left RE
        $clamp            # type of GC clamp
       ) = @_;
    
    my $left_primers  = {};
    
    if ($clamp eq 'both') {            
        _create_primers($re_l, $modified_gene, $left_primers, 'both');                
    } else {            
        _create_primers($re_l, $modified_gene, $left_primers, 'left');
    }
 
       
    my @left;

    foreach my $l (keys %{$left_primers}) {    
  
        if (_check_3prime_gc($l) == 1) {
            push @left, $l
        }
    }
 
    return @left;
}



#### create right primers ####
sub _create_right_primers {
    
    my ($more_modified_gene,   # ATGC modified for left and right RE
        $re_pos_r,             # position of right RE
        $re_r,                 # right RE
        $clamp                 # type of GC clamp
       ) = @_;

    my $right_primers = {};
                   
    if ($clamp eq 'both') {                    
        _create_primers($re_r, $more_modified_gene, $right_primers, 'both');
    } else {                    
        _create_primers($re_r, $more_modified_gene, $right_primers, 'right');
    }
                    
    my @right;
     
    foreach my $r (keys %{$right_primers}) {  
        push @right, $r 
    }           
    
    # reverse complement right primer
    my @rev_comp;
    foreach my $r (@right) {
        $r =~ tr/ATGC/TACG/; 
        $r = reverse $r;
        if (_check_3prime_gc($r) == 1) {
            push @rev_comp, $r;
        }
    }
    
    return @rev_comp;
}



#### check 3' end for undesirable [GC] run ####
sub _check_3prime_gc {

    my ($primer   # 5' to 3' order for left or right
       ) =@_;

    my $num_at_end = 5;
    my $end = substr($primer, (length $primer) - $num_at_end, $num_at_end);

    if ($end =~ /[GC][GC][GC]/) { 
            
        # undesirable GC run found at 3' end   
        return 0 
        
    } else { 
  
        return 1 
        
    }
}



#### find in-frame start codon ####
sub _start_codon {

    my ($gene      # ATGC sequence
       ) = @_;

    my $gene_start = 0;
    
    if ($gene =~ /^((.{$CODON_SIZE})*?)((ATG)|(GTG))/) { 
        $gene_start = $-[3];
    }

    return $gene_start;
}



#### find in-frame stop codon location ####
sub _stop_codon {

    my ($gene,         # ATGC sequence
        $gene_start    # look for stop after start
       ) = @_;

    my $WAY_TOO_BIG = 100000000; # bigger than any anticipated pattern
    my $gene_stop = $WAY_TOO_BIG;

    # look for stop codon, keep track of first one in sequence after start codon
    foreach my $stop_codon (('TAA', 'TAG', 'TGA')) {

        if (substr($gene, $gene_start) =~ /^((.{$CODON_SIZE})*?)($stop_codon)/) {

            if ($-[3] < $gene_stop) { 
                $gene_stop = $-[3] 
            }
        }
    }

    # sanity check if stop codon was found
    if ($gene_stop == $WAY_TOO_BIG) { 
        return undef 
    } else { 
        return $gene_stop + $gene_start
    }
}


1;

__END__

=head1 NAME

Bio::MCPrimers (and mcprimers.pl)
 
=head1 DESCRIPTION

Creates molecular cloning PCR primer pairs for a given gene so that the
gene can be directionally inserted into a vector. Solver is generic, 
restriction enzymes and their order in the vector are specified in the 
caller.
 
=head1 EXPORT SUBROUTINES

sub find_mc_primers

    $gene,           # ATGC string (use 21 NT upstream, 200 NT downsteam)
    $flag_hr,        # anonymous hash reference to flags from mcprimers.pl
    $version_sr,     # version scalar reference returned to caller   
    $ecut_loc_ar,    # enzyme cut location array reference from caller
    $vcut_loc_ar,    # vector cut location array reference from caller
    @re              # array of restriction enzyme strings from caller
    
Not explicitily exported. Use Bio::MCPrimers::find_mc_primers

See mcprimers.pl for an example of use.
 
=head1 INSTALLATION

    MCPrimers.pm     - place into Perl Bio/MCPrimers.pm
    CloningVector.pm - place into perl Bio/Data/Plasmid/CloningVector.pm
    mcprimers.pl     - place in a directory where it can be accessed by users. 

    PRIMER3DIR       - set environment variable to point to Primer3 
                       executable directory.

    If PRIMER3DIR is not set, MCPrimers will set it to '.' by default.

    MSWindows -    use primer3.exe
    Other     -    use primer3_core
 
=head1 DEPENDENCIES

Primer3 used as primer3.exe on MSWindows and as primer3_core otherwise.

Specify environment variable PRIMER3DIR for path to Primer3 executable directory.

=head1 SYNOPSIS using mcprimers.pl

perl mcprimers.pl -h

perl mcprimers.pl -vector pET-32a ppib.fa > ppib.pr3

perl mcprimers.pl -vector pET-32a -clamp 3prime hemy.fa > hemy.pr3

perl mcprimers.pl -shift 12 -vector pET-32a -clamp 3prime cyss.fa > cyss.pr3

Note: Use -Ilib if modules are still in local lib directory.

See mcprimers.pl for an example of the use of Bio::MCPrimers itself
 
=head1 LIMITATIONS

Limitation: Primer3 does not account for redundancy codes.

Note:       Runs use intermediate files, so don't use this in a
            manner such that intermediate files will overwrite one another. 
            CGI use is probably a bad idea.
            
There is no guarantee this code will find the best solution, or even any 
solution, or that the solutions it finds will be correct or useful to you.

Use at your risk. Check any solutions you obtain.
 
=head1 BUGS

Probably. Use at your own risk.

This software comes with no guarantee of usefulness. 
Use at your own risk. Check any solutions you obtain. 
Stephen G. Lenk assumes no responsibility for the use of this software.
 
=head1 COPYRIGHT

Stephen G. Lenk (C) 2005, 2006. All rights reserved. 

This program is free software; you can redistribute it and/or 
modify it under the same terms as Perl itself.

Primer3 is Copyright (c) 1996,1997,1998,1999,2000,2001,2004

Whitehead Institute for Biomedical Research. All rights reserved.

=head1 AUTHOR

Stephen G. Lenk, November 2005, May 2006. 

slenk@emich.edu
 
==head1 CHANGES

V1.5.0
Solid base code - start change log with next release

V1.5.1
Package updated to beter reflect good CPAN release practices.

=head1 ACKNOWLEDGEMENTS

Primer3 is called by this code to verify that the PCR primers are OK.

Thanks to Tim Wiggin for algorithm suggestions and encouragement. 
The use of direct string comparisons. 
Modify gene to match PCR primer for Primer3 check. 
    
Thanks to Dan Clemans for showing me molecular cloning in the first place. 
I am using Dr. Clemans's ideas about good MC primers in this code. 
Any errors in interpretation or implementation are mine.

Patricia Sinawe found that earlier versions of MCPrimers did not
detect out-of-frame solutions and suggested that extra binding sites
could be included.

Ken Youens-Clark <kyclark@gmail.com> has provided guidance in the proper 
naming of this software so that it functions cooperatively with other 
Perl modules.

Other references:

(1) http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html

(2) http://www.mcb.uct.ac.za/pcroptim.htm


=cut