package Bio::Data::Plasmid::CloningVector; 
our $VERSION = '1.00';

# Stephen G. Lenk (C) 2006. All rights reserved. 

# This program is free software; you can redistribute it and/or 
# modify it under the same terms as Perl itself.

# Licenced under the Artistic Licence.

# This software comes with no guarantee of usefulness. 
# Use at your own risk. Check any solutions you obtain. 
# Stephen G. Lenk assumes no responsibility for the use of this software.


# Use: cloning_vector_data ( $vector_name,  # vector name (i.e. pET-32a)
#                            $re_ra,        # restriction enzyme pattern
#                            $re_name_rh,   # names of restriction enzymes
#                            $ecut_loc_ra,  # cut location in enzyme
#                            $vcut_loc_ra ) # cut location in vector

# Returns: 0=fails (vector not found), 1=success

 
use strict;
use warnings;

sub cloning_vector_data {

    my ( $vector_name,  # vector name (i.e. pET-32a)
         $re_ra,        # restriction enzymes pattern
         $re_name_rh,   # names of restriction enzymes
         $ecut_loc_ra,  # cut location in enzyme
         $vcut_loc_ra   # cut location in vector
       ) = @_;  


    # check name of cloning vector
    if ($vector_name eq 'pET-32a') {

        push @{$re_ra}, 'TTCGAA';      # NspV       (TT^CGAA) (labrat site) 
        push @{$re_ra}, 'AGATCT';      # BglII      (A^GATCT)
        push @{$re_ra}, 'GGTACC';      # KpnI       (GGTAC^C)
        push @{$re_ra}, 'CCATGG';      # NcoI       (C^CATGG)
        push @{$re_ra}, 'GATATC';      # EcoRV      (GAT^ATC)
        push @{$re_ra}, 'GGATCC';      # BamHI      (G^GATCC)
        push @{$re_ra}, 'GAATTC';      # EcoRI      (G^AATTC)
        push @{$re_ra}, 'GAGCTC';      # SacI       (GAGCT^C)
        push @{$re_ra}, 'GTCGAC';      # SalI       (G^TCGAC)
        push @{$re_ra}, 'AAGCTT';      # HindIII    (A^AGCTT)
        push @{$re_ra}, 'GCGGCCGC';    # NotI, EagI (GC^GGCCGC)
        push @{$re_ra}, 'CTCGAG';      # AvaI, XhoI (C^TCGAG)   
    
        $re_name_rh->{TTCGAA}   = 'NspV';
        $re_name_rh->{AGATCT}   = 'BglII';
        $re_name_rh->{GGTACC}   = 'KpnI';
        $re_name_rh->{CCATGG}   = 'NcoI';
        $re_name_rh->{GATATC}   = 'EcoRV';
        $re_name_rh->{GGATCC}   = 'BamHI';
        $re_name_rh->{GAATTC}   = 'EcoRI';
        $re_name_rh->{GAGCTC}   = 'SacI';
        $re_name_rh->{GTCGAC}   = 'SalI';
        $re_name_rh->{AAGCTT}   = 'HindIII';
        $re_name_rh->{GCGGCCGC} = 'NotI, EagI';
        $re_name_rh->{CTCGAG}   = 'AvaI, XhoI'; 
    
        push @{$ecut_loc_ra}, 2;
        push @{$ecut_loc_ra}, 1;
        push @{$ecut_loc_ra}, 5;
        push @{$ecut_loc_ra}, 1;
        push @{$ecut_loc_ra}, 3;
        push @{$ecut_loc_ra}, 1;
        push @{$ecut_loc_ra}, 1;
        push @{$ecut_loc_ra}, 5;
        push @{$ecut_loc_ra}, 1;
        push @{$ecut_loc_ra}, 1;
        push @{$ecut_loc_ra}, 2;
        push @{$ecut_loc_ra}, 1;    
    
        push @{$vcut_loc_ra}, 2;
        push @{$vcut_loc_ra}, 0;
        push @{$vcut_loc_ra}, 2;
        push @{$vcut_loc_ra}, 2;
        push @{$vcut_loc_ra}, 0;
        push @{$vcut_loc_ra}, 1;
        push @{$vcut_loc_ra}, 1;
        push @{$vcut_loc_ra}, 2;
        push @{$vcut_loc_ra}, 2;
        push @{$vcut_loc_ra}, 2;
        push @{$vcut_loc_ra}, 0;
        push @{$vcut_loc_ra}, 2;  

        return 1;

    }
    else {

        return 0;
    }   
}

return 1;
