#!/usr/bin/perl -w
############################## INTERACTION ACCUMULATOR #########################
#
my $versionID = '1.0.5';
#
#
# Uses a downloaded version of a recent biogrid database to calculate
# interaction distances bewteen elements in an ordered list of genes.  The
# order can be defined by almost anything, but this has been used mostly
# for screen results. The results are written to two separate output
# lists. One is a so-called short list with interactions listed by type in
# a matrix following the caluclated numbers. The second so-called long
# list just has the interaction elements, distances and the system and
# author reported for the data.
#
#   CURRENT BIOGRID FILE
#
#     BIOGRID-ORGANISM-Saccharomyces_cerevisiae-2.0.60.tab.txt
#
#		download updated files at http://www.thebiogrid.org/downloads.php
#



use Data::Dumper; # used only if read out-all file fails
use Storable qw(store retrieve freeze thaw dclone);
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use Getopt::Long;

############################     COMMAND LINE OPTIONS       ######################## 

my $writeRank = ''; # Command line option to store the ranked orf list
my $randomize = ''; # Command line option to run randomization of a list
my $strong = ''; # Command line to list reciprocal interactions in output
my $long = ''; # Command line option to print system and source in output
my $millFile = ''; # command line option to get gene list from screenMill out-all file
my $max = ''; # CL switch to trim promiscuous interactions - default is 400 if switch is used.
#                 A number can be provided on CL if a value other than 400 is needed. This 
#                 parameter sets the maximum number of interactions per gene that are considered
#                 in the plots. Genes with more interactions are deleted from the interaction data.
my $maxDefault = 400;  # default value for max interactions
my $help = ''; # print help screen
my $version = ''; # prints version number

# randomization will affect the following: Output file names (_R) will be appended if 
#  $randomize is TRUE, Randomization subroutine will be invoked
# The strong option will cause the interaction strengthening sub to run on the interaction hash
#  This was originally meant to bump up reciprocals when deleting singleton interactions
#  but I havent coded the dismissal of singletons
# The long option will result in long output including source and system

GetOptions ('randomize' => \$randomize,'strong' => \$strong,'long' => \$long,
							'millfile' => \$millFile,'max' => \$max, 'help' => \$help,
							'write' => \$writeRank, 'version' => \$version);

my($OUT2,$DEBUG,$OUT3);

##############################    First, see if help is required  ##################################


if ($help) {

	&printHelp;
	exit;

	}

##############################    Next, see if version # request  ##################################


if ($version) {

	print "MinACC version $versionID\n";
	exit;

	}


################################  I/O Calls ######################################

###########  SCREEN DATA FILES



print "reading screen data...\n";

my $screen_data_file = $ARGV[0];

die "Oops! A file called '$screen_data_file' does not exist.\n" unless -e $screen_data_file;


## grab the root of the input file name to use for naming output files

my @file_name_parts = split /\./ , $screen_data_file;
my $screen_root;

print "$file_name_parts[0]\t $file_name_parts[1]\n";

if ($file_name_parts[0] ne '') {

	$screen_root = $file_name_parts[0];
	
	}
	
else {

	$screen_root = $screen_data_file;
	
	}
	

# The default read data subroutine '&read_ORF_list' is a quick and dirty input for a 1-column list
#  for a 1-column list of ORFS. HEADERS should be preceeded by '!' and the ORFS are in rank 
#  order based on a screen metric.  Sreen data from a ScreenMILL 'out-all' file is read using the
#  '&read_screen_data' subroutine.

my @ranked_ORFs;

if ($millFile) { # if millfile command line switch was used, get ORF list from screenMill out-all file
	@ranked_ORFs = &read_screen_data($screen_data_file);
	}
else {
	@ranked_ORFs = &read_ORF_list($screen_data_file);
	}

print "done reading list of $#ranked_ORFs\n";




###########  LOAD INTERACTION DATA FILES


my $BIOGRID_file = "BIOGRID-ORGANISM-Saccharomyces_cerevisiae-2.0.63.tab.txt";

print "loading BIOGRID data...\n";

my ($counted_interactions_ref, $interactions_ref, $all_sources_ref, $all_systems_ref) = &get_BIOGRID($BIOGRID_file,@ranked_ORFs);


########### PROMPT USER FOR NOISE REDUCTION VALUE IF SWITCH USED

####  Need to add validation for a minimum

if ($max) {
	print "Enter a maximum number of interactions per gene (default $maxDefault):";
	$_ = <STDIN>;
	chomp;
	$max = $_ ? $_ : $maxDefault; # sets to default if no value
	}




###########  OUTPUT FILES

#### Tag output file names based on CL switches

my $outFile = ">".$screen_root."_output";
my $rankFile = ">".$screen_root."_ranks";
if ($randomize) {$outFile .= "_R";$rankFile .= "_R";}
if ($strong) {$outFile .= "_S";$rankFile .= "_S";}
if ($long) {$outFile .= "_L";}
if ($max) {$outFile .= "_NR$max";$rankFile .= "_NR$max";}
$outFile .= ".txt";
$rankFile .= ".txt";

print "Output file name is $outFile\n";


####################################################################################################
#                                                MAIN
####################################################################################################


#***************************************************************************************************
#***************************************** OUTPUT FILES   ******************************************
#***************************************************************************************************

open ($OUT2, $outFile) || die "Could not open output file. $!";


#####  RANDOMIZE ranked list if the 'randomize' switch was used on CL.

if ($randomize) {
	print "ORF list is being randomized for control graph\n";
	@ranked_ORFs = @{&Fisher_Yates_Shuffle(\@ranked_ORFs)};
	@ranked_ORFs = @{&Fisher_Yates_Shuffle(\@ranked_ORFs)};
	}

##### IDENTIFY PROMISCUOUS INTERACTORS


if ($max) {

#	&scrubCasualEncounters ($interactions_ref,\@ranked_ORFs,$max);
	&trimNoisyItems ($interactions_ref,\@ranked_ORFs,$max);
	print "Noise reduction accomplished\n";

	}

##### Dump ranked list (with randomization, noisy interactions removed, etc) to file

if ($writeRank) {
	
	open ($OUT3, $rankFile) || die "Could not write the ranked ORFs to disk\n";
	foreach (@ranked_ORFs) {
		
		print $OUT3 "$_\n";
		
		}
	close $OUT3;
	
}

#####  RANK the ORFs in a numbered list


my %rank = map { $ranked_ORFs[$_], $_ + 1 } 0 .. $#ranked_ORFs;

## The problem with the above statmenent is that if there are repeated
## ORFs in the ranked list, an ORF will take on the highest value of the ranks.
## This is an issue if your affected (interesting) ranks are near 1, and you have an
## affected gene that has a non-phenotypic duplicate further down in the ranked list.
## e.g with TOP1-TA, MMS4 has 3 library copies, one of which is unaffected. Need to fix
## this issue!



########### OUTPUT
#
#  make a list with each interaction printed separately for generation of a histogram
#
#		we can build in an extra print for strong interactions

#exit;

my $header_long = "rank\tORF_A\tORF_B\tdistance\tsigned\trect\tsystem\tsource";
my $header_short = "rank\tORF_A\tORF_B\tdistance\tsigned\trect";

if ($long) {
	print $OUT2 "$header_long\n";
	}
else {
	print $OUT2 "$header_short\n";
	}

foreach my $A (@ranked_ORFs) {#   A's ->
	
	print "$A\t";

	my @list = keys %{${$interactions_ref}{$A}};

	foreach my $B (@list) {#      -> B's
	
		print "$B\t";

		my $difference = abs ($rank{$A} - $rank{$B});# position of A - position of B
		
		my $signed = $rank{$B} - $rank{$A}; # signed value for new graph
		
		my $rect = $signed + $rank{$A} -1;

		if ($difference==0) {next;} #ignore self interactions

		my @system = keys %{${$interactions_ref}{$A}->{$B}};

		foreach my $SYS (@system) { # 	-> systems
		
			print "$SYS\t";

			my $out_SYS = $SYS;
			$out_SYS =~ tr/ -//d;

			my @source = keys %{${$interactions_ref}{$A}->{$B}->{$SYS}};

			foreach my $SOUR (@source)  {#		-> source
			
				print "$SOUR\n";

				my $out_SOUR = $SOUR;
				$out_SOUR =~ tr/ -//d;
				
				if ($long) {
					print $OUT2 "$rank{$A}\t$A\t$B\t$difference\t$signed\t$rect\t$out_SYS\t$out_SOUR\n";
					}
				else {
					print $OUT2 "$rank{$A}\t$A\t$B\t$difference\t$signed\t$rect\n";
					}
				# at some point we need to remove special chars from the system and source for output
				# this would be friendlier to R. use tr/ -//d

			}# end foreach source

		}# end foreach system
		
		if (defined(${$counted_interactions_ref}{$B}->{$A}) && $strong) {# if there is a reciprocal interaction B -> A
		
			my $recipDistance = $signed; # keep the reciprocal value the same 
			if ($long) {
				print $OUT2 "$rank{$A}\t$A\t$B\t$difference\t$recipDistance\treciprocal\treciprocal\n";
				}
			else {
				print $OUT2 "$rank{$A}\t$A\t$B\t$difference\t$signed\t$rect\n";
				}
		}# end if reciprocal

	}# end foreach B

}# end foreach A

close $OUT2;

##################### End of output section #################

exit;

###################################################################################################
#                                     END OF MAIN
###################################################################################################


###################################################################################################
#                                     SUBROUTINES
###################################################################################################

################################### GET_BIOGRID
#
# Biogrid interaction data is acquired from a recent organism-specific set of interaction
# data (download at http://www.thebiogrid.org/downloads/BIOGRID-ORGANISM-x.x.x.x etc). Then
# the interactions for the list of genes in your distribution/screen are identified and
# stored in a hash for the main program.
#
# INPUTS:  $BIOGRIDfilename, @rankedORFs
#
# RETURN VALUES:  $counted_interactions_ref, $interactions_ref, $all_sources_ref, $all_systems_ref

sub get_BIOGRID {

	my $BIOGRID_file = shift; # pull file name off of parameter list
	
	my (%outAll_ORFs,$BIOGRID);  
	
	map {$outAll_ORFs{$_}++} @_; # move remaining parameters into hash for easy lookup below (note this is just an ordered list of ORFs)

	open ($BIOGRID, "<".$BIOGRID_file) || die "Couldn't open BIOGRID data\n";

	my (%biogrid_colindex, $iloop, $head, @headers);
	
	$/ = line_break_check( $BIOGRID ); # check what type of line breaks this file contains

	####### get rid of header in interaction file (Biogrid)
	
	while(! (defined $biogrid_colindex{'interactor_a'} && defined $biogrid_colindex{'interactor_b'})){
		$iloop++;
		if($iloop==100){ die "Your data file is not formatted properly.$!";}
		chomp ($head = <$BIOGRID>);
		#if($iteration=~/1st/){print temp_DATA "$head\n";}
		chomp(@headers = split /\t/, "\L$head"); # want it to be case insensitive 
		%biogrid_colindex=();
		@biogrid_colindex{@headers} = (0..$#headers);
	}
	%biogrid_colindex=();
	@biogrid_colindex{@headers} = (0..$#headers);
	
	####### now iterate over BIOGRID file to accumulate interactions in LIST
	
	my (%interactions,@data,%counted_interactions, $temp_counter);
	my (%all_sources, %all_systems, @all_systems_list, @all_sources_list, %source_rank, %system_rank, %source_colindex,%system_colindex);
	# the above declarations are for variables to catch all sources & systems and apply an ID number to them
	
	foreach(<$BIOGRID>){
	# 	$temp_counter++;
	# 	print "$temp_counter\n";
		chomp;
		@data=split /\t/;
		if($data[$biogrid_colindex{'interactor_a'}]							# if A and B both exist and are not '' and both A & B are in dataset then add to data structure
			 && $data[$biogrid_colindex{'interactor_a'}] ne '' 
			 && $data[$biogrid_colindex{'interactor_b'}] 
			 && $data[$biogrid_colindex{'interactor_b'}] ne ''
			 && $outAll_ORFs{$data[$biogrid_colindex{'interactor_a'}]} 
			 && $outAll_ORFs{$data[$biogrid_colindex{'interactor_b'}]} ){


				$interactions{$data[$biogrid_colindex{'interactor_a'}]}->{$data[$biogrid_colindex{'interactor_b'}]}->{$data[$biogrid_colindex{'experimental_system'}]}->{$data[$biogrid_colindex{'source'}]}++;
				$all_sources{$data[$biogrid_colindex{'source'}]}++;
				$all_systems{$data[$biogrid_colindex{'experimental_system'}]}++;
				$counted_interactions{$data[$biogrid_colindex{'interactor_a'}]}->{$data[$biogrid_colindex{'interactor_b'}]}++;

		}#	end if
	
	}  # end foreach
	close $BIOGRID;

	return (\%counted_interactions, \%interactions, \%all_sources, \%all_systems); #return hash refs



} # end sub get_BIOGRID


################################### RANDOMIZE AN ARRAY
#
#
sub Fisher_Yates_Shuffle {
	my ($list) = @_;
	my @array = @$list;
	my $i;
	for ($i = @array; --$i; ) {
			my $j = int rand ($i+1);
			next if $i == $j;
			@array[$i,$j] = @array[$j,$i];
	}
	return \@array;
}

################################### LINE_BREAK_CHECK
#
# line_break_check receives a file handle as its input and returns the new line character used in the file handle
sub line_break_check{
	my $file = shift;
	local $/ = \1000; # read first 1000 bytes
	local $_ = <$file>; # read
	my ($newline) = /(\015\012?)/ ? $1 : "\012"; # Default to unix.
	seek $file,0,0; # rewind to start of file
 	return $newline;
}

################################### MIN
#
#
sub min{
	my($a,$b)=@_;
	if($a<$b){return $a;}
	return $b;
}

################################### TRAVERSE_SUM
#
# 	This subroutine recursively follows a complex data structure and sums the values
# 	it finds in any arrays or hashes that it encounters. Hash or array limbs to a
# 	data tree are detected by a ref call, but other refs (CODE, GLOB, SCALAR & REF)
# 	are ignored.  May want to build in handling of scalars at some point. This borrows
# 	heavily from Data::Dumper - but I could not figure out how to make Data::Dumper do the
# 	same job.

sub traverse_sum {

#	print "@_\n";
	
	my $y=shift;
	my @list;
	my $sum;
	
# 	my $a =  ref $y;
# 	
# 	print "$a\n";
	
	if (ref $y eq 'HASH') {
	
		@list = keys (%{$y});
		
		foreach my $element (@list) {
		
			if (ref $y->{$element})  {
			
				$sum += traverse_sum ($y->{$element});
				
			}
			elsif (looks_like_number($y->{$element})) {
			
				$sum += $y->{$element};
				
			} # else it is not a ref or a number
			
		} # end foreach
		
	} #end if
	
	elsif (ref $y eq 'ARRAY') {
	
		@list = @{$y};
		
		foreach my $element (@list) {
		
			if (ref $element) {
			
				$sum += traverse_sum ($element);
				
			}
			elsif (looks_like_number($element)) {
		
			$sum += $element;
			
			} 
			
		} # end foreach
		
	} # end elsif
	
	return $sum;
	
} # end sub

################################### READ_SCREEN_DATA
#
#
sub read_screen_data {
	
	my $SCREEN_DATA;
	my $screen_data_file = $_[0]; # note parameters always passed as an array

	open ($SCREEN_DATA, "<".$screen_data_file) || die "Could not open screen data file. $!";

	$/ = line_break_check( $SCREEN_DATA ); # check what type of line breaks this file contains

	my $iloop=0; #used to check if we are in a seemingly infinite loop

	my (%outAll_colindex, @outAll_data, @outAll_ORFs,@headers,$head,$orf_col, $log_ratio_col, $v, @temp);
	
	#
	# get rid of header in an out-all file from a screen
	#
	
#	while(!(defined $outAll_colindex{'description'} && defined $outAll_colindex{'plate #'} && defined $outAll_colindex{'row'} && defined $outAll_colindex{'col'} && defined $outAll_colindex{'p-value'} && (defined $outAll_colindex{'calculated log ratio (control::exp)'} || defined $outAll_colindex{'calculated log ratio (control::exp)'}) && (defined $outAll_colindex{'sgd'} || defined $outAll_colindex{'orf'}) ) ){
	while(!(defined $outAll_colindex{'plate #'} && defined $outAll_colindex{'row'} && defined $outAll_colindex{'col'}) ){
		$iloop++;	
		if($iloop==100){print Dumper(\%outAll_colindex); die "Your out_all data file is not formatted properly.$!";}
		chomp ($head = <$SCREEN_DATA>);
		#if($iteration=~/1st/){print temp_DATA "$head\n";}
		$head =~ s/\s+\t|\t\s+/\t/; # trim off expra spaces before or after column labels v1.05
		@headers = split /\t/, "\L$head"; # want it to be case insensitive 
#		print "@headers\n";
		if($headers[0] && $headers[0]=~/p-value cut-off for significance:/){$v->{'pthresh'}=$headers[1];}
		%outAll_colindex=();
		@outAll_colindex{@headers} = (0..$#headers);

		# In older out-all files the column label is 'col' in newer files it is 'column'
	
		if ($outAll_colindex{'column'}) {
			$outAll_colindex{'col'} = $outAll_colindex{'column'};
			}


	}
	%outAll_colindex=();
	@outAll_colindex{@headers} = (0..$#headers); #		store header line in array
	
	# depending on how old the out-all data file is orfs and log_ratios may have different identifiers
	
	if ($outAll_colindex{'id column'}) {
		$orf_col = $outAll_colindex{'id column'};
		}
	elsif ($outAll_colindex{'orf'}) {
		$orf_col = $outAll_colindex{'orf'};
		}
	else {
		die "could not find the ORF column in your data file\n";
		}
	$log_ratio_col = $outAll_colindex{'calculated log ratio (control::exp)'} ? $outAll_colindex{'ratio'} : $outAll_colindex{'calculated log ratio (control::exp)'};


	print "header line identified:\t$headers[0]\t $headers[1]\t $headers[2]\tetc...\n";
	print "ORF column is $orf_col\n";

	
	#
	# now that we have removed and analyzed the header, iterate over the rest of an out-all file
	#
	
#	@temp = sort { $b->[$outAll_colindex{'z-score'}+1] <=> $a->[$outAll_colindex{'z-score'}+1] }
#	map { [ $_, (split /\t/) ] } grep(!/BLANK|excluded|dead/igo,<$SCREEN_DATA>); # slurp in the data, ignore blanks, deads, and excludes.  tab separated fields

	# V1.05
	@temp = sort { $b->[$outAll_colindex{'z-score'}+1] <=> $a->[$outAll_colindex{'z-score'}+1] }
	map { [ $_, (split /\t/) ] }
	#grep(!/BLANK-.*-BLANK|HIS3-BLANK|excluded-.*-excluded|dead-.*-dead/igo,<$SCREEN_DATA>);
	grep(!/BLANK\-.*\-BLANK|HIS3\-BLANK|excluded\-.*\-excluded|dead\-.*\-dead|^\n+$|^\s+$|^\t+$/igo,<$SCREEN_DATA>);
	#

	my @zscores = map { $_->[$outAll_colindex{'z-score'}+1] } @temp;
	@outAll_data = map { $_->[0] } @temp;
	@outAll_ORFs = map { $_->[$orf_col+1] } @temp;
	close $SCREEN_DATA;
	
	# outAll_data && outAll_ORFs is now sorted by z-score...

	return @outAll_ORFs;

}

################################### READ_ORF_LIST
#
#
# sub read_ORF_list {
# 	my $filename = shift; # take off parameter
# 	my (@list, @fail, $INFILE, $count);
# 	open ($INFILE, "<".$filename) || die "Could not open input file. $!";
# 	$/ = line_break_check( $INFILE ); # check what type of line breaks this file contains
# 	LINE: while (<$INFILE>) {
# 		$count++;
# 		next LINE if /^\!/; # gives it a miss with preceeding ! (meant for headers)
# 		chomp;
# 		if (/Y[A-P][LR]\d{3}[WC](-[A-D])?/) {
# 			push (@list,$_);
# 		} # end if
# 		else {
# 			push (@fail,$count);
# 		} # end else
# 	}# end while
# 	if (@fail) {
# 		print "Non ORF data exists in your input list at the following lines:\n";
# 		foreach my $x (@fail) {
# 			print "$x\n";
# 			}
# 		print "These lines were not included for interaction processing\n";
# 				} # end if fail
# 	return @list;
# }

sub read_ORF_list {
    my $filename = shift; # take off parameter
    my (@list, @fail, $INFILE, $count);    
    my $orf_exact_pattern='^[Y|y][A-Pa-p][L|R|l|r][0-9]{3}[W|w|C|c](?:-[A-Za-z])?$';

    open ($INFILE, "<".$filename) || die "Could not open input file. $!";
    $/ = line_break_check( $INFILE ); # check what type of line breaks this file contains

    LINE: while (<$INFILE>) {
        $count++;
        next LINE if /^\!/; # gives it a miss with preceeding ! (meant for headers)
        chomp;
        $_ =~ s/\s+//g; # strip out spaces
        warn "$_\n";

        # if we encounter an ORF name

        if (/$orf_exact_pattern/) {push (@list,$_);    } # end if

        else {push (@fail,$count);} # end else

    }# end while

    if (@fail) {

        print "Non ORF data exists in your input list at the following lines:\n";
        foreach my $x (@fail) {print "$x\n";    }
        print "These lines were not included for interaction processing\n";

    } # end if fail

    return @list;
}


################################### NOISE REDUCTION I : noisy query interactions A -> #
#
#

sub removeCasualEncounters {

	my $interactionsHashRef = shift;
	my $rankedORFsListRef = shift;
	my $max = shift;
	my %encounterHash = ();
	my @valueList = ();
	open ($DEBUG, ">debug.txt") || die "could not open the debug file for writing\n";
	
	foreach my $A (@{$rankedORFsListRef}) {#   A's ->
	
		#print "$A\t";
	
		my $count = scalar keys %{${$interactionsHashRef}{$A}};
		
		if ($count > $max) {
		
			$encounterHash{$A} = $count;  # put size of list into hash of ORFs in list
			delete $interactionsHashRef->{$A};
			print $DEBUG "$A\t$count\n";

			}
		
		
		} # end foreach
	
	close $DEBUG;
}

################################### NOISE REDUCTION II : noisy queries & reciprocals
#															A -> #,  # -> A
#
#

sub scrubCasualEncounters {

	my $interactionsHashRef = shift;
	my $rankedORFsListRef = shift;
	my $max = shift;
	my %encounterHash = ();
	my @valueList = ();
	open ($DEBUG, ">debug.txt") || die "could not open the debug file for writing\n";
	
	###  First run through BIOGRID hash structure is to ID and remove noisy queries.
	###  Query IDs are then used in second pass to remove interactions where the noisy
	###  ORF is detected as the 'B' ORF in an A->B pair
	
	foreach my $A (@{$rankedORFsListRef}) {#   A's ->
	
		#print "$A\t";
	
		my $count = scalar keys %{${$interactionsHashRef}{$A}};
		
		if ($count > $max) {
		
			$encounterHash{$A} = $count;  # put size of list into hash of ORFs in list
			delete $interactionsHashRef->{$A};
			print $DEBUG "$A\t$count\n";

			}
		
		
		} # end foreach
	
	###  Second level of noise reduction
	
	foreach my $A (keys %{$interactionsHashRef}) {
		
		foreach my $B (keys %{$interactionsHashRef->{$A}}) {

			if ($encounterHash{$B}) {
				
				delete $interactionsHashRef ->{$A} ->{$B};
				
				}

			} # end foreach
		
		} # end foreach
	
	close $DEBUG;
}
################################### NOISE REDUCTION III : remove noisy ORFs from ordered list
#
#
sub listNoisyQueries {

	my $interactionsHashRef = shift;
	my $max = shift;
	my %encounterHash = ();

	foreach my $A (keys %$interactionsHashRef) {

		my $count = scalar keys %{${$interactionsHashRef}{$A}};
		
		if ($count > $max) {
			
			$encounterHash{$A} = $count;
			
			}

		}

	return \%encounterHash;

}

sub trimNoisyItems {

	my $interactionsHashRef = shift;
	my $rankedORFsListRef = shift;
	my $max = shift;
	my @hashOut = ();
	
	open (my $DEBUG1, ">listDeltas.txt")||die "Can't open debug 1\n";
	open (my $DEBUG2, ">listCollapsed.txt")||die "Can't open debug 2\n";
	
	#########  Find the most promiscuous interactors

	my $encounterHashRef = &listNoisyQueries ($interactionsHashRef,$max);
	
	@hashOut = %$encounterHashRef;
	
	#########  Delete these from the rank-ordered list
	
	my $listSize = scalar @{$rankedORFsListRef};
	
	foreach (my $x=0;$x<$listSize;$x++) {
	
		if ($encounterHashRef->{$rankedORFsListRef->[$x]}) {# if the list element is in the noise list

			delete $rankedORFsListRef->[$x]; # then delete from list
			
			} # end if
			
		} # end foreach
		
	print $DEBUG1 "@hashOut\n";

	@$rankedORFsListRef = grep { defined } @$rankedORFsListRef; # collapses UNDEFS out of list

	print $DEBUG2 "@$rankedORFsListRef\n";	
	
	#########  Then spin through %interactionsHash to remove noise
	#
	#  Try a trimming routine vs. just re-calculating these from the abbreviated list
	
	# Remove noisy interactors as queries by deleting hash keys
	
	foreach my $N (keys %{$encounterHashRef}) {

		delete $interactionsHashRef->{$N};

		}
	
	# Then remove noisy guys as hits

	foreach my $A (keys %{$interactionsHashRef}) {

		foreach my $B (keys %{$interactionsHashRef->{$A}}) {

			if ($encounterHashRef->{$B}) {

				delete $interactionsHashRef ->{$A} ->{$B};

				}

			} # end foreach

		} # end foreach

	close $DEBUG1;
	close $DEBUG2;

}

################################### PRINT HELPSCREEN
#
#

sub printHelp {

print "
Mini Interaction Accumulator (minAcc) $versionID takes as input an ordered list of ORFs
derived from a functional screen performed on a library of gene disruption
strains. An interaction density profile is determined by measuring the distance
(difference in list position) between interacting ORF pairs in the list.
Interactions are taken from the BIOGRID database using the S. cerevisiae -
specific file of interactions.

CURRENT BIOGRID FILE:

BIOGRID-ORGANISM-Saccharomyces_cerevisiae-2.0.60.tab.txt

download updated files at http://www.thebiogrid.org/downloads.php

The biogrid file must be located in the same directory(folder) as the minAcc
program. The orf list file is entered on the command line. The following
switched can be entered at the command line to change program behavior.

 --randomize : randomizes the orf list as a control interaction comparison.
 --strong    : lists reciprocal interactions in addition to standard output.
 --long      : prints system and source as columns in output file.
 --millFile  : uses a screenMill \'out-all\' file as the input gene list. 
                 Default behavior is to use a simple one-colum list of ORFs 
                 with header lines preceeded by a \'!\'.
 --max       : trims promiscuous BIOGRID interactors from list before distance calculations.
                 Default is no trimming, setting --max will default to trim genes with 400
                 or more interactions from the list, but the user is prompted to enter a
                 maximum value.
 --write     : writes the rank-ordered list of ORFs as a single column text file.
 --help      : print this help screen.
 --ver       : prints version number.
";

}

####################################### Version Info
#
#			Version 1.0
#
#				Adds job control switches using long options (i.e. --)
#
#			Version 1.0.1
#
#				Adds bug fix to check for 'column' vs. 'col' as headings in screenMill out-all files
#
#			Version 1.0.2
#
#				failed bug fix
#
#			Version 1.0.3
#
#				fixes interaction timming by deleting interactors from hash elements
#				  as well as from the rank order list
#
#			Version 1.0.4  NOT IMPLEMENTED YET
#
#				Adds a switch to generate x number of randomized interaction lists for 
#				  background stats
#
#			Version 1.0.5
#
#				Bug fixes from John in read screen file sub and read orf list sub





#######################################   LEFTOVER - but possibly useful

#### I found this useful enough to leave the code fragment in...
# 
# foreach my $zed (@ranked_ORFs) {
# 
# 	print $OUT3 "$zed\t$rank{$zed}\n";
# 	
# }# end foreach
# 
# exit;




################  Get counts of ORF occurences
# my %histogram;
# 
# $histogram{$_}++ for @ranked_ORFs;
# 
# 
# foreach my $x (keys %histogram) {
# 
# 	print $OUT3 "orf $x has $histogram{$x} occurences\n";
# 	
# }
# 
# exit;

#print "orfs are:\n@ranked_ORFs\n";



################### Former Main output section #####################

# printed out the interactions with numbers of interactions

# my @new_systems_list = @all_systems_list;
# foreach (@new_systems_list) {tr/ -//d}; #remove spaces and '-' for R
# 
# my @systems_headers = join ("\t",@new_systems_list);
# 
# 
# 
# my $header = "rank\tORF_A\tORF_B\tno_interactions\tstrong_interactions\tdistance\tinverse\tsubtraction\t@systems_headers";
# 
# 
# print $OUT "$header\n";
# 
# foreach my $x (@ranked_ORFs) {#   A's ->
# 
# #	print $OUT "$rank{$x}\t$x\t";
# 
# 	
# 	my @list = keys %{${$counted_interactions_ref}{$x}};
# 
#  	foreach my $y (@list) {#      -> B's
# 
# 		my $length = $#all_systems_list;
# 		my @systems_array= map 0,0..$#all_systems_list; # zero the systems array
#  		my $inverse;
#  		my $difference = abs ($rank{$x} - $rank{$y});# position of A - position of B
#  		my $subtraction = $#ranked_ORFs - $difference + 1;
#  		if ($difference==0) {next;}
#  		else {$inverse=1/$difference;}
#  		
#  		#
#  		#	need to add granular accumulation of systems
#  		#		do a foreach w/ system names as hash key and number= traverse_sum(interaction->A->B->{system})
#  		#		this sum gets pasted into an output array by array[$system_colindex{'system'}]
#  		#		then join that sum array with tabs and append to output line
#  		#
#  		foreach my $experiment (@all_systems_list) { # fill the systems array
#  		
#  			if (${$interactions_ref}{$x}->{$y}->{$experiment}) {
#  			
# # 				print "$experiment position $system_colindex{$experiment}\n";
#  			
#  				$systems_array[$system_colindex{$experiment}] = traverse_sum(${$interactions_ref}{$x}->{$y}->{$experiment});
#  			
#  			}
#  		
#  		}
#  		
# #		print "2. systems array is then @systems_array\n";
# 
#  		my @systems_output = join ("\t",@systems_array);
#  		
#  		print $OUT "$rank{$x}\t$x\t$y\t${$counted_interactions_ref}{$x}->{$y}\t$strong_interactions{$x}->{$y}\t$difference\t$inverse\t$subtraction\t@systems_output\n";
#  	}
# 
# #	print $OUT "\n";
# }
# 
# close $OUT;

# 
# my ($OUT,$OUTNAME);
# 
# my $storage_hash = $GO_file.$screen_data_file;
# $storage_hash =~ tr/\./\_/;
# 
############################# interaction strengthening #######

# The idea here is that there are many A->B interactions that are singletons, i.e.
# one query gene resulting in many interactions.  These can be filtered out in the
# accumulated data output by removing the points with only one interaction.  However,
# if A->B exists and B->A exists each as a singleton, i.e. a bipolar interaction,
# then this should be considered stronger evidence of interaction.  Here we will
# iterate over the two interaction hashes and add one interaction to bipolars. So,
# if A->B & B->A exists, then add one to both so that A->B = 2 and B->A = 2.
# This data gets listed as a separate column in the output, so graphs can be
# generated with or without this data.

############################

# Strengthening for $counted_interactions copied and recorded into %strong_interactions
#		
#		strengthens by reciprocal even if reciprocal is from a different system or source

#print Dumper (%counted_interactions);

# store ($counted_interactions_ref, "temphash");
# 
# my %strong_interactions = %{retrieve("temphash")}; # storable module for deep cloning of structure
# 
# # print "Now STRONG interactions\n";
# # print Dumper (%strong_interactions);
# # exit;
# 
# print "Calculating reciprocal interactions...\n";
# 
# foreach my $a (keys %strong_interactions) {#  A loop
# 
# 	my @B_list = keys %{$strong_interactions{$a}};
# 
# 	foreach my $b (@B_list) {
# 	
# 		if ($strong_interactions{$b}->{$a}) {
# 		
# #			print "$strong_interactions{$a}->{$b}\n";
# 		
# 			$strong_interactions{$a}->{$b}++;
# 			
# #			print "interaction between $a and $b is reciprocal\n";
# 			
# #			print "$strong_interactions{$a}->{$b}\n";
# 
# 
# 		} #end if exists accumulation
# 	
# 	} # end foreach B
# 
# 
# } # end foreach A
# 
# print "Done with reciprocals.\n";


# Strengthening for $interactions by system & source is NOT IMPLEMENTED YET
#		
#		strengthens by reciprocal only if in the same system and source

# foreach my $system (keys %nteractions) {# system loop
# 
# 	my @source_list=keys %{${$interactions}{$system}};
# 	
# 	foreach my $source(@source_list) {
# 	
# 		my @A_list = keys %{${$interactions}{$system}->{$source}};
# 		
# 		foreach my $a (@A_list) {
# 		
# 			my @B_list = keys %{${$interactions}{$system}->{$source}->{$a}};
# 			
# 			foreach my $b (@B_list) {
# 			
# 				if (${$interactions}{$system}->{$source}->{$b}->{$a}) {
# 				
# 					${$interactions}{$system}->{$source}->{$a}->{$b}++;
# 				
# 				} #end if exists accumulation
# 			
# 			} # end foreach B
# 		
# 		} # end foreach A
# 	
# 	} # end foreach source
# 
# }# end foreach system

#
########### end of strengthening
#
#### The following is meant to be used for a table output in which the interactions
#      are listed with multiple interactions & sources for each record. The statements
#      sort the sources & systems lists from BIOGRID and makes column indices based on these
#      lists so that a table can include all the possible types of interactions for each
#      A -> B interacting pair
####
# my (@all_systems_list, @all_sources_list, %source_rank, %system_rank, %source_colindex,%system_colindex);
# 
# @all_sources_list = sort keys %{$all_sources_ref}; # nothing implemented with this yet
# @all_systems_list = sort keys %{$all_systems_ref};
# 
# %source_colindex=();
# @source_colindex{@all_sources_list} = (0..$#all_sources_list);
# 
# %system_colindex=();
# @system_colindex{@all_systems_list} = (0..$#all_systems_list);
####



#print "orfs are:\n@ranked_ORFs\n";


