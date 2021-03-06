# Gets information about CUIs in a NodeInfo file including:
#
# CUI, preferred term, semantic type(s), semantic group(s)
#
# Thoutput file contains these values tab seperated, if there are multiple
# answers for a value (e.g. multiple semantic groups) then those answers are 
# comma seperated. Info about each CUI is on a single line.
#
# Example Output:
# C0000545    Fish Oil   T047,T053     DYSN     
#
#  TODO - update example above to make sure it is correct/accurate
#
# Usage: perl getCuiInfo.pl [nodeInfoFile] [cuiInfoOutFile]
# Example: perl getCuiInfo.pl RayFish_nodeInfo RayFish_cuiInfo
#
use strict;
use warnings;
use UMLS::Interface;


#user input
my $nodeInfoFile = shift;
my $outFile = shift;


##################################
#        Begin Code
##################################
########
#initialize UMLS::Interface
my %tHash = ();
$tHash{'t'} = 1; #default hash values are with t=1 (silence module output)

$tHash{'hostname'} = '192.168.24.89';
$tHash{'username'} = 'henryst';
$tHash{'password'} = 'OhFaht3eique';
$tHash{'socket'} = '/var/run/mysqld.sock';

my $componentOptions = \%tHash;
#add more component options here (if needed)
#else use default configuration
my $umls = UMLS::Interface->new($componentOptions) 
    or die "Error: Unable to create UMLS::Interface object.\n";

########
#read in all cuis from the cui file
my %cuis = ();
open IN, $nodeInfoFile or die ("ERROR: cannot open nodeInfoFile: $nodeInfoFile"); 
while (my $line = <IN>) {
    while ($line =~ /(C\d{7})/) {
	$cuis{$1} = 1;
    }
}
close IN;


########
# get all preferred terms for all cuis
my %preferredTerm = ();
foreach my $cui(keys %cuis) {
    #get the preferred term
    my $preferredTermString = $umls->getAllPreferredTerm($cui);

    #set to undef if no preferred string exists
    if ($preferredTermString eq '') {
	$preferredTermString = 'undef';
    }

    #save the preferred term
    $preferredTerm{$cui} = $preferredTermString;
}

########
# get semantic types for all cuis
my %semanticType = ();
foreach my $cui (keys %cuis) {
    #get all semantic types
    my $semanticTypeString = '';
    my $typesRef = $umls->getSt($cui);
    foreach my $type(@{$typesRef}) {
	my $abr = $umls->getStAbr($type);
	$semanticTypeString .= "$abr,";	
    }
    chop $semanticTypeString;

    #set to undef, if no type exists
    if ($semanticTypeString eq '') {
	$semanticTypeString = 'undef';
    }

    #save the semantic type
    $semanticType{$cui} = $semanticTypeString;
}

########
# get semantic groups of each semantic type
my %semanticGroup = ();
foreach my $cui (keys %semanticType) {
    
    #find the semantic group if it hasn't been already
    if(!defined $semanticGroup{$semanticType{$cui}}) {

	#get the semantic groups
	my $groups = $umls->stGetSemanticGroup($semanticType{$cui});
	my $semanticGroupString = '';
	foreach my $group (@{$groups}) {
	    $semanticGroupString .= "$group,";
	}
	chop $semanticGroupString;

	#set to under if no group exists
	if ($semanticGroupString eq '') {
	    $semanticGroupString = 'undef';
	}
    }
}

########
# Output the results
open OUT, ">$outFile" or die ("ERROR: unable to open outFile: $outFile\n");
foreach my $cui (keys %cuis) {
    print OUT "$cui\t$preferredTerm{$cui}\t$semanticType{$cui}\t$semanticGroup{$semanticType{$cui}}\n";
}
close OUT;

#Done!
print "Done!\n";
