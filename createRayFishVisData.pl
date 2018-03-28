#script to perform discovery replication and create output for visualization

use strict;
use warnings;

##### Perform Clustering
my $vclusterLocation =  'clustering/cluto-2.1.1/Linux/vcluster';
my $lbdFile = '../data/discoveryRep/rayFish_ltc_threshold_targetTermList';
my $vectorFile = '../data/discoveryRep/1983_1985_window8.vectors';
my $clusterMethod = 'rb';
my $assocMeasure = 'x2';
my $assocMatrix = '/home/sam/clusterVisualization/data/discoveryRep/1983_1985_window8.matrix';
my $outputFile = 'TEST';
my $outputDir = '../output/discoveryRep/rayFish/';
my $parentArrayOut = 'RayFish_parentArray';
my $clustersFileOut = 'RayFish_clusters';
my $cuiInfoFileOut = 'RayFish_cuiInfo';

my $command = "perl clustering/createVisulizationOutput.pl $vclusterLocation $lbdFile $vectorFile $clusterMethod $assocMeasure $assocMatrix $outputDir $parentArrayOut $clustersFileOut";
print STDERR "\n\n$command\n\n";
`$command`;

$command = "perl general/getCuiInfo.pl $outputDir$clustersFileOut $outputDir$cuiInfoFileOut";
print "$command\n";
`$command`;


print "DONE!\n";



