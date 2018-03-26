# program that takes an output file and generates a linking metabolite file
# the inMatrix is used to find the linking metabolites. It is the metabolite 
# by disease matrix used to generate the ltc scores.
# only outputs linking metabolites for cuis that have scores greater than 0
use strict;
use warnings;
use UMLS::Interface;

#user input
my $inFile = 'output/cardiacArrestDiseases.output';
my $inMatrix = 'metaboliteDisease_1975_2015_window8_threshold1';
my $inMetabolites = 'cuis/metabolites.cuis';
my $outFile = 'output/cardiacArrestDiseases.metabolites';


#####################################
#####################################


# Input Checking
open IN, $inFile or die ("ERROR: cannot open inFile: $inFile\n");
open IN_MATRIX, $inMatrix or die ("ERROR: cannot open inMatrix: $inMatrix\n");
open IN_METABOLITES, $inMetabolites or die ("ERROR: cannot open inMetabolites: $inMetabolites\n");
open OUT, ">$outFile" or die ("ERROR: cannot open outfile: $outFile\n");


#initialize UMLS::Interface
print "   Initializing UMLS\n";
my %tHash = ();
$tHash{'t'} = 1; #default hash values are with t=1 (silence module output)
my $componentOptions = \%tHash;
# use default configuration
my $umls_interface = UMLS::Interface->new($componentOptions) or die "Error: Unable to create UMLS::Interface object.\n";


###############################
##  BEGIN CODE
###############################

#Read the inFile into a hash{$cui} = $line 
# and save the order of the output (for ties) as $cuis in @order
# each line of the in file is: score - cui - term(s)
print "   Reading Output File\n";
my %inLines = (); #save lines for output
my %scores = (); #save scores as a check
my @order = (); #save order for output (ties)
while (my $line = <IN>) {
    #grab the score and cui from the line
    $line =~ /(-?\d+) - (C\d{7}) -/;
    my $score = $1;
    my $cui = $2;

    #error checking
    if (!defined $score || !defined $cui) {
	print "ERROR reading line in outputFile: $line\n";
	exit;
    }

    #only output the scores that are > 0
    if ($score > 0) {
	#save the lin and the order
	$inLines{$cui} = $line;
	$scores{$cui} = $score;
	push @order, $cui;
    }
}
close IN;


#read in all the metabolite cuis
# cuis file contains a cui on each line (thats it)
print "   Reading MetabolitesFile\n";
my %metabolites = ();
while (my $line = <IN_METABOLITES>) {
    chomp $line;
    $metabolites{$line} = 1;
}
close IN_METABOLITES;

#construct a list of metabolite terms
my %metaboliteTerms = ();
foreach my $metabolite (keys %metabolites) {
    my $preferredTerm = $umls_interface->getAllPreferredTerm($metabolite);
    $metaboliteTerms{$metabolite} = $preferredTerm;
}

#find all the linking metabolites for each of the cuis
# matrix file contains cui\tcui\tfrequencyOfCooccurrence
# the result of this is a linking CUIs list which is of the format:
# hash{$disease} = a hash of the form $hash{$metabolite} = $score
# e.g. %linkingCuis is a hash of hashes that specifies the linking metabolites
#      and their frequency of co-occurrence for each disease
print "   Reading Matrix\n";
my %linkingCuis = ();
while (my $line = <IN_MATRIX>) {
    chomp $line;
    (my $cui1, my $cui2, my $score) = split(/\t/, $line);
    
    #find the linking CUI and disease (first or second)
    # or this is not a metabolite line, so we can ignore it
    my $metabolite;
    my $disease;
    if (exists $inLines{$cui1}) {
	if (exists $metabolites{$cui2}) {
	    $metabolite = $cui2;
	    $disease = $cui1;
	}
    }
    if (exists $inLines{$cui2}) {
	if (exists $metabolites{$cui1}) {
	    $metabolite = $cui1;
	    $disease = $cui2;
	}
    }

    #if the line did contain a linking CUI (metabolite) and disease
    if (defined $metabolite && defined $disease) {
	#add the linking CUI to the list of linking CUIs
	if (!exists $linkingCuis{$disease}) {
	    my %emptyHash = ();
	    $linkingCuis{$disease} = \%emptyHash;
	}
	
	${$linkingCuis{$disease}}{$metabolite} = $score;
    }
}
close IN_MATRIX;


#output the linking terms file
# the file is the form: 
# CUI line = numLinkingTerms - CUI - Term
# metabolites line(s) = \t $frequencyOfCo-occurrence $CUI - Term
print "   Outputting File\n";
foreach my $cui (@order) {
    #print out the CUI line
    print OUT $inLines{$cui};


    #grab the linking metabolites hash
    my $metabolitesRef = $linkingCuis{$cui};
    
    #check to make sure the number of metabolites in the metabolitesRef
    # are the same as the score (which should be the LTC)
    if ($scores{$cui} != scalar keys %{$metabolitesRef}) {
	die ("ERROR: disease score does not match the number of metabolites: $cui, $scores{$cui}, ".(scalar keys %{$metabolitesRef})."\n");
    }

    #print out each of the linking term lines (in descending by frequeny of 
    # co-occurrence)
    foreach my $metabolite (sort {${$metabolitesRef}{$b} <=> ${$metabolitesRef}{$a}} keys %{$metabolitesRef}) {
	print OUT "\t${$metabolitesRef}{$metabolite} - $metabolite - $metaboliteTerms{$metabolite}\n";
    }

}
close OUT;

print "DONE!\n";

