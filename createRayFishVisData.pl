#script to perform discovery replication and create output for visualization

use strict;
use warnings;

##### Perform Clustering
my $vclusterDir =  'clustering/cluto-2.1.1/Linux';
my $lbdFile = '../data/discoveryRep/rayFish_ltc_threshold_targetTermList';
my $vectorFile = '../data/discoveryRep/1983_1985_window8.vectors';
my $clusterMethod = 'rb';
my $outputFile = 'TEST';
my $outputDir = '../output/discoveryRep/rayFish/';
my $parentArrayOut = 'RayFish_parentArray';
my $clustersFileOut = 'RayFish_clusters';
my $cuiInfoFileOut = 'RayFish_cuiInfo';

my $command = "perl clustering/createVisualizationOutput.pl $vclusterDir $lbdFile $vectorFile $clusterMethod $outputDir $parentArrayOut $clustersFileOut";
print "$command\n";
#`$command`;

$command = "perl general/getCuiInfo.pl $outputDir$clustersFileOut $outputDir$cuiInfoFileOut";
print "$command\n";
`$command`;


print "DONE!\n";



