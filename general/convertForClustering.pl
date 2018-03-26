#converts this output to clustering output (standard LBD output)
# does not convert any terms over that have 0 linking terms 
# this is indicated by a score of -1 in the output file
use strict;
use warnings;

my $inFile = 'output/cardiacArrestDiseases.output';
my $outFile = 'output/cardiacArrestDiseases.forClustering';

open IN, $inFile or die ("ERROR: unable to open inFile: $inFile\n");
open OUT, ">$outFile" or die ("ERROR: unable to open outFile: $outFile\n");

#write a header for the output
print OUT "converted metabolite data from $outFile\n\n";

#read the infile and write to the output
my $rank = 1;
while (my $line = <IN>) {
    #each line is <score> - <CUI> - <preferredTerm>
    $line =~ /-?(\d+\.?\d*) - (C\d{7}) - (.*)/;
    my $score = $1;
    my $cui = $2;
    my $preferredTerm = $3;


    #vals = score - CUI - preferred term
    if ($score > 0) {
	print OUT "$rank\t$score\t$cui\t$preferredTerm\n";
	$rank++;
    }
}
close IN;
close OUT;

print "done!\n";


