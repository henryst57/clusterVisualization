#reads the disease file and metabolites file and makes pairs for
# the two files. This pair file can then be input into UMLS::Association

#User input
my $metaboliteFile = 'cuis/cardiacArrest.cui';
my $diseaseFile = 'cuis/dysn.cuis';
my $pairsFile = 'cardiacArrestDisease.pairs';


#test the I/O
open M_IN, $metaboliteFile or die("ERROR: cannot open metaboliteFile: $metaboliteFile\n");
open D_IN, $diseaseFile or die("ERROR: cannot open diseaseDile: $diseaseFile\n");
open OUT, ">$pairsFile" or die ("ERROR: cannot open outputFile: $pairsFile\n");


#read in all the metabolites, each line is a CUI
my @metabolites = ();
while (my $line = <M_IN>) {
    chomp $line;
    push @metabolites, $line;
}
close M_IN;

#read in all the diseases, each line is a CUI
my %diseases = ();
while (my $line = <D_IN>) {
    chomp $line;
    push @diseases, $line;
}
close D_IN;

#output the metabolite disease pairs, each line is a CUI<>CUI
foreach my $metabolite (@metabolites) {
    foreach my $disease (@diseases) {
	print OUT "$metabolite<>$disease\n";
    }
}
close OUT;

print "DONE!\n"; 

