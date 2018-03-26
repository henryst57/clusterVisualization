#reduces the co-occurrence matrix to be only a matrix of the scores of interest

my $pairsFile = 'pairs/metaboliteDisease.pairs';
my $matrixInFile = '1975_2015_window8_threshold1';
my $matrixOutFile = 'metaboluteDisease_1975_2015_window8_threshold1';

#input checking
open P_IN, $pairsFile or die ("ERROR: unable to open pairs file: $pairsFile\n");
open M_IN, $matrixInFile or die ("ERROR: unable to open matrix file: $matrixInFile\n");
open OUT, ">$matrixOutFile" or die ("ERROR: unable to open out file: $matrixOutFile\n");


#read all the pairs into a hash
my %pairs = ();
while (my $line = <P_IN>) {
    chomp $line;
    $pairs{$line} = 1;
}
close P_IN;

#copy each line of interest to the new matrix file
while (my $line = <M_IN>) {
    (my $cui1, $cui2, $val) = split(/\t/,$line); 
    if (exists $pairs{"$cui1<>$cui2"}) {
	print OUT $line;
    }
}
close M_IN;
close OUT;
