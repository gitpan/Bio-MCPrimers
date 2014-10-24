#!/usr/bin/perl -w

# Stephen G. Lenk (c) 2006
# LICENCED UNDER PERL ARTISTIC LICENCE
# See Bio/Mcprimers.pm for POD

$VERSION = '2.3'; 

use Tk;
use IPC::Open3;
require Tk::FileSelect;
require Tk::ROText;

use strict;
use warnings;

###################################################################
######                  ######
###### Global Variables ######
######                  ######
##############################

my $mw;                        # main window
my $app_title = 'MCPrimers';

my $menubar_frame;             # frames
my $right_frame;
my $bottom_frame;
my $para_frame;
my $scale_frame;
my $p3_frame;
my $fname_frame;
my $vec_frame;
my $text_frame;
my $seq_frame;
my $reset_frame;

my $text_message;              # widgets
my $sub_button;
my $vec_message;
my $seq_text;
my $filedialog;
my @site_cb;
my $p3_check;
my $fname_message;

my $clamp = 'both';            # parameters
my $shift = 0;                

my $vector_loaded = 0;         # state
my $seq_loaded = 0;
my $p3_used = 0;
my @text;                      # stores text printed to window

my $initial_text='Please load vector file and sequence first';
my $help_text;

my $vector_name = '';          # vector file name
my @re;                        # restriction enzymes pattern
my %re_name;                   # names of restriction enzymes
my @ecut_loc;                  # cut location in enzyme
my @vcut_loc;                  # cut location in vector
my @sites = ();                # recognition sites array
my @result = ();               # from MCPrimers
my $seq_name = '';             # sequence file
my $p3_name = '';              # primer3 file
my $num_solns = 0;             # number of solutions
my @seq = ();                  # sequence (FASTA)
my $mcprimers_dir;             # mcprimers.pl directory

####################################################################

BEGIN {

    # check MCPRIMERS_DIR
    if (defined $ENV{MCPRIMERS_DIR}) {
        $mcprimers_dir = $ENV{MCPRIMERS_DIR};
    }  
    else {
        $mcprimers_dir = '.';
    }
}

END {
    # intentionally empty		
}

###################################################################
######                 ######
###### Main code block ######
######                 ######
#############################

&set_wd();

$mw = MainWindow->new(-title => $app_title);
$filedialog=$mw->FileSelect(-title     => 'Select File', 
                            -width     => 40, 
                            -takefocus => 1);
&define_help_text();
&create_text();
&create_right();
&create_buttons();
&create_seq();
&display_wd();

$mw->resizable(0,0);
$mw->configure(-takefocus);

MainLoop;


###################################################################

sub create_text {
        
    $text_frame = $mw->Frame->pack(-pady=>5, -side=>'top', -anchor=>'nw');
    $text_message = $text_frame->Message(
        -width => 600,
        -text  => "First load vector file and DNA sequence",
        -font  => 'Helvetica 12 bold',
        );
    $text_message->pack(-anchor=>'nw', -padx=>20);
}

###################################################################

sub create_right {
        
    $right_frame = $mw->Frame->pack(-pady=>0, -side=>'right', -anchor=>'nw');
	
 	$fname_frame = $right_frame->Frame->pack(-padx=>5,
                                             -pady=>2,
                                             -anchor=>'nw');
	$fname_message = $fname_frame->Message(
        -width => 200,
        -text  => "FASTA: undefined",
        -font  => 'Helvetica 10 bold',
        )->pack(-anchor=>'nw');    
    $para_frame = $right_frame->Frame->pack(-padx=>5,
                                            -pady=>2,
                                            -anchor=>'nw');
    $para_frame->Message(
        -text  => "GC clamping",
        -font  => 'Helvetica 10 bold',
        -width => 200,
        )->pack(-anchor=>'nw');
    $para_frame->Radiobutton (
        -text     => 'GC clamp 3\' and 5\'',
        -variable => \$clamp,
        -value    => 'both',
        -font => 'SmallItem',)->pack(-padx => 5, -anchor=>'nw');
    $para_frame->Radiobutton (
        -text     => 'GC clamp 3\' only  ',
        -variable => \$clamp,
        -value    => '3prime',
        -font => 'SmallItem',)->pack(-padx => 5, -anchor=>'nw');
 

    $p3_frame = $right_frame->Frame->pack(-padx=>5,
                                          -pady=>2,
                                          -anchor=>'nw');
    $p3_frame->Message(
        -text  => "Primer3 file",
        -font  => 'Helvetica 10 bold',
        -width => 200,
        )->pack(-anchor=>'nw');
    $p3_check = $p3_frame->Checkbutton (
        -text     => 'Unspecified',
        -variable => \$p3_used,
        -state    => 'disabled',
        -font => 'SmallItem',)->pack(-padx => 5, -anchor=>'nw');
       

    $scale_frame = $right_frame->Frame->pack(-padx=>5,
                                             -pady=>2,
                                             -anchor=>'nw');
    $scale_frame->Message(
        -text  => "Search shift",
        -font  => 'Helvetica 10 bold',
        -width => 200,
        )->pack(-anchor=>'nw');  
    $scale_frame->Scale(
        -resolution   => 1,
        -from         => 30,
        -to           => 0,
        -variable     => \$shift,
        -sliderlength => 10,
        -length       => 140,
        -orient       => 'horizontal',
        -showvalue    => 1,
        -state        => 'active',
        -font         => 'Helvetica 10 bold',)->pack(-padx => 5, -anchor=>'nw');										 
											 										 
    $vec_frame = $right_frame->Frame->pack(-padx=>5,
                                           -pady=>2,
                                           -anchor=>'nw');   
    $vec_message = $vec_frame->Message(
        -width => 200,
        -text  => "\n",
        -font  => 'Helvetica 10 bold',
        )->pack(-anchor=>'nw'); 
}

###################################################################

sub create_buttons {
    
    $bottom_frame = $mw->Frame->pack(-pady=>15, -side=>'bottom');
    
    $bottom_frame->Button(
        -text    => 'Load Vector',
        -bg      => 'light green',
        -command => \&load_vector_file,
        )->pack(-anchor=>'w', -side=>'left', -padx=>2);
        
    $bottom_frame->Button(
        -text    => 'Load FASTA',
        -bg      => 'light green',
        -command => \&load_seq_file,
        )->pack(-anchor=>'w', -side=>'left', -padx=>2);
        
    $bottom_frame->Button(
        -text    => 'Select Primer3',
        -bg      => 'light green',
        -command => \&load_primer3,
        )->pack(-anchor=>'w', -side=>'left', -padx=>2);        
          
    $sub_button = $bottom_frame->Button(
        -text    => 'Execute',
        -bg      => 'light gray',
        -command => \&execute_mcprimers,
        )->pack(-anchor=>'w', -side=>'left', -padx=>2);
        
    $bottom_frame->Button(
        -text    => 'Save',
        -bg      => 'light green',
        -command => \&save_window,
        )->pack(-anchor=>'w', -side=>'left', -padx=>2);
                 
    $bottom_frame->Button(
        -text    => 'Help',
        -bg      => 'light green',
        -command => \&help,
        )->pack(-anchor=>'w', -side=>'left', -padx=>2); 
        
    $bottom_frame->Button(
        -text    => 'About',
        -bg      => 'light green',
        -command => \&about,
        )->pack(-anchor=>'w', -side=>'left', -padx=>2);        

    $bottom_frame->Button(
        -text    => 'Exit',
        -bg      => 'light green',
        -command => sub { exit; },
        )->pack(-anchor=>'w', -side=>'left', -padx=>2);            
}

###################################################################

sub create_seq {
 
    $seq_frame = $mw->Frame->pack(-pady=>5, -side=>'bottom');       
    $seq_text = $seq_frame->Scrolled('ROText');
    $seq_text->configure( -scrollbars => 'e',
                          -height     => 42,
                          -width      => 100,
                          -font       => 'Courier 10',
                          -background => 'white');
    $seq_text->pack(-padx=>15);
}

###################################################################

sub execute_mcprimers {
 
	&clear_text();   
    my $excluded = '';
    @result = (); 
    
    if ($vector_loaded and $seq_loaded) {
			
        $text_message->configure(-text=>
          "Executing MCPrimers - please wait (probably several minutes)");
        $mw->update;
        
        # get excluded sites list
        my $i = 0;
        foreach (@sites) { 
            if (${$_} == 0) { 
                $excluded .= (',' . $re_name{$re[$i]});
            }
            $i += 1;    
        }
        
        my $x = '';
        if ($excluded ne '') { 
            $x = "-excludedsites=" . substr($excluded, 1);
        }      
        my $c = "-clamp=$clamp";
        my $s = "-shift=$shift";
        my $p = '';
        if (defined $p3_name and $p3_name ne '' and $p3_used == 1) {
            $p="-primerfile=\"$p3_name\"";
        }
		
		my $cmd;
        if ($^O =~ /^MSW/) { 
				
			# Microsoft
		    my $start = "perl $mcprimers_dir\\mcprimers.pl -stdout";  
            $cmd = "$start $x $c $p $s \"$vector_name\" \"$seq_name\" ";
            my $pid = open3(\*WTRFH, \*RDRFH, \*RDRFH, $cmd);	
		}	
		else {
				
			# non-Microsoft
		    my $start = "perl $mcprimers_dir/mcprimers.pl -filter";  
            $cmd = "$start $x $c $p $s \"$vector_name\" ";
            my $pid = open3(\*WTRFH, \*RDRFH, \*RDRFH, $cmd);
            foreach (@seq) {
                print WTRFH $_;
            }   
            close WTRFH;
        }
		print "$cmd\n";		
        $mw->update;
		
 		my $m;      
        my $line  = <RDRFH>;
		
RESULTS:        
        while (defined $line) {
			
			chop($line);
			if ($^O =~ /^MSW/) { chop $line; }
            push @result, "$line\n";
			
            if ($line =~ /Sorry/   or 
			    $line =~ /Error/   or 
				$line =~ /Problem/ or 
				$line =~ /not available/) { 
						
                $m = $line;
				last RESULTS;
            }
            if ($line =~ /.*Solution # (\d*)/) {
                $m = "Done - $1 solutions found";
            }
            $line  = <RDRFH>;
        }
        close RDRFH;
        
        write_text(@result);
        $text_message->configure(-text=>$m);
        $mw->update;
        
    }
    else {
        $text_message->configure(-text=>$initial_text);
    }
}

###################################################################

sub load_vector_file {
        
	&clear_text();
    if ($vector_loaded == 1) {
    
        # can't figure out how to properly destroy old widgets
        # ->destroy isn't good enough
        $text_message->configure(-text=>"Sorry - can\'t reload vector file");
        return;
    }
    
    $text_message->configure(-text=>"Select vector file");
    $vector_loaded = 0;
	
    my $status;
    
    # popup get file name
    my $old_name = $vector_name;
    $vector_name = $filedialog->Show;
    if (defined $vector_name and $vector_name ne '') { 
    
        # details of the plasmid used as a vector
        use Bio::Data::Plasmid::CloningVector;  
        $status = 
           Bio::Data::Plasmid::CloningVector::cloning_vector_data
             ($vector_name, \@re, \%re_name, \@ecut_loc, \@vcut_loc);
    }
    else {
        $status = 0;
    }
    
    if ($status == 0) {
			
        my $m = '';
        if (defined $vector_name) { 
            $m = "Error: Can\'t load \'$vector_name\'";
        }
        else {
            $m = 'Load vector cancelled';
        }
        $text_message->configure(-text=>$m);
		
        $vector_name = $old_name;
    }
    else {
        foreach (@re) { 
            {
                my $f = 1;
                push @sites, \$f;
                push @site_cb, $vec_frame->Checkbutton
                ( -text     => $re_name{$_},
                  -variable => \$f,
                  -font => 'SmallItem')->pack(-side   =>'top',
                                              -anchor => 'nw',  
                                              -padx   => 5);
            }
        }       

        $vector_loaded = 1;

        $vector_name =~ /.*\/(.+)/;
        my $short_name = $1;

        $text_message->configure(-text=>"Vector file $short_name loaded");
        if ($vector_loaded and $seq_loaded) {
            $sub_button->configure(-bg=>'light green');     
        }
		$short_name =~ /(.+)\..+/;
        $vec_message->configure(-text  => "\n$1 cloning sites",); 
    }
}

###################################################################

sub load_seq_file {
        
    $text_message->configure(-text=>"Select FASTA file");
    my $status;
    my $line;
    my $seq_fh;
	@seq = ();
	&clear_text();  
	
    # popup get file name
    my $old_name = $seq_name;
    $seq_name = $filedialog->Show;
    if (defined $seq_name and $seq_name ne '') { 
    
        # read sequence file
        open $seq_fh, $seq_name;
        my $t = <$seq_fh>;
        while (defined $t) {
            push @seq, $t;
            $t = <$seq_fh>;
        }
        $status = 1;    
    }
    else {
        $status = 0;
    }
    
    if ($status == 0) {
        my $m = '';
        if (defined $seq_name) { 
            $m = "Error: Can\'t load \'$seq_name\'";
        }
        else {
            $m = 'Load sequence cancelled';
        }
        $text_message->configure(-text=>$m);

        $seq_name = $old_name;
    }
    else {
        write_text(@seq);       
        $seq_loaded = 1;

        $seq_name =~ /.*\/(.+)/;
        my $short_name = $1;

        $text_message->configure(-text=>"FASTA file $short_name loaded");
		$fname_message->configure(-text=>"FASTA: $short_name");
    }        

    if ($vector_loaded and $seq_loaded) {
        $sub_button->configure(-bg=>'light green');
    }   
}

####################################################################

sub clear_text {

    $seq_text->selectAll; 
    $seq_text->deleteSelected;
}

####################################################################

sub write_text {

    @text = @_;      

    &clear_text();
    $mw->update;
	
    foreach (@text) {
		my $l = $_;
        $seq_text->Insert($l);
    }

    $mw->update;
    $seq_text->yviewMoveto(0.0);
}

####################################################################

sub about {

    $text_message->configure(-text=>
      "MCPrimers (c) 2006, Steve Lenk (and Tim Wiggin for CloningVector)");
}

####################################################################

sub help {

	&clear_text();
    my @h = ();
    foreach (split "\n", $help_text) {
        push @h, "$_\n";
    }
    &write_text(@h);       
}

###################################################################

sub save_window {

    $text_message->configure(-text=>"Select file for saving window text");

    my $file_name = $filedialog->Show; 

    if (defined $file_name and $file_name ne '') {

        my $fh;
 
        $file_name =~ /.*\/(.+)/;
        my $short_name = $1;

        open $fh, ">$file_name" or 
          $text_message->configure(-text=>"Unable to open $file_name");

        my $status = 1;
        
SAVE:   foreach (@text) {
            if ($status == 1) {
                print $fh $_ or $status = 0;
            }
            else {
                last SAVE;
            }
        }

        close $fh;

        if ($status == 0) {
            $text_message->configure(-text=>"Unable to print to $file_name");
        }
        else {
            $text_message->configure(-text=>"Window text saved to $short_name");
        }
    }
    else {
        $text_message->configure(-text=>"Save cancelled");
    }
}

####################################################################

sub load_primer3 {

    my $fh;
    my @tmp = ();
    my $p3s_found = 0;
        
    $text_message->configure(-text=>"Select Primer3 file");

    my $old_name = $p3_name;
    $p3_name = $filedialog->Show; 

    if (defined $p3_name and $p3_name ne '') {
    
        open $fh, "$p3_name" or 
          $text_message->configure(-text=>"Unable to access $p3_name"); 
     
        my $line  = <$fh>;

        while (defined $line) {

            if ($line =~ /.+=.+/) {
                $p3s_found += 1;
            }
            push @tmp, $line;
            $line  = <$fh>;
        }

        if (@tmp) { 

            $p3_name =~ /.*\/(.+)/;
            my $short_name = $1;
            my $w = '';

            if ($p3s_found == 0) { 
                $w = 'WARNING,'; 
            }

            $text_message->configure(
              -text=>"Primer3 file \'$short_name\' - $w $p3s_found PRIMER tag(s) identified");    
            $p3_used = 1;    
            $p3_check->configure(-text  => $short_name, 
                                 -state => 'active',);
        }
    }
    else {
        $text_message->configure(-text=>"Primer3 file load cancelled");
        $p3_name = $old_name;
    }

    &write_text(@tmp); 
}

####################################################################

sub define_help_text {

$help_text = qq/
MCPrimers GUI Help Text - V2.3

A vector (plasmid) file and a FASTA file must be loaded before execution.

- The Execute button is greyed out until a vector file and a FASTA file are loaded.
  
- Load a vector file with the Vector button.
  This provides information about the vector.
  This will make the cloning site checklist appear.
  Checked sites are used. 
  Unchecked sites are not used.
  All sites are checked initially.

- Load an in-frame FASTA file with the FASTA button.
  Best to use 21 NT (or more) upstream and 200 downstream.
  The FASTA will load into the text window.

- The Primer3 button allows you to specify a Primer3 parameter file.
  When you have specified the file, a CheckButton will appear
  with the file name. You can select if you wish to use the file.

- The Execute button will turn green when a vector file and a plasmid file are specified.

- You can execute MCPrimers with the Execute button.
  After several minutes, the solution will appear in the text window.
  You can scroll through the solution.

- The Save button saves the text window contents.

- The About button gives author-copyright information.

- The Help button displays this text.

- The Exit button exits the GUI.

Enjoy!

Steve Lenk (c) 2006
/;

}

####################################################################

sub set_wd {
		
    if ($^O =~ /^MSW/) { 
    
        # Microsoft
		my $dir = "$ENV{USERPROFILE}/desktop";
        chdir($dir);
    } 
    else {         
    
        # non-Microsoft
        chdir();
    }
}

####################################################################

sub display_wd {

    use Cwd;
	my $dir = &getcwd();
	$text_message->configure(-text=>"Current working directory is $dir");
}

####################################################################
