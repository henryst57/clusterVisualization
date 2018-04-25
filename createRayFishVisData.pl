#script to perform discovery replication and create output for visualization

use strict;
use warnings;

##########################################
#                User Input
##########################################
#clustering params
my $vclusterLocation =  'clustering/cluto-2.1.1/Linux/vcluster';
my $lbdFile = '../data/discoveryRep/rayFish_ltc_threshold_targetTermList';
my $vectorFile = '../data/discoveryRep/1983_1985_window8.vectors';
my $clusterMethod = 'rb';

#association params (ranking clusters)
my $assocMeasure = 'x2';
my $assocMatrix = '../data/discoveryRep/1983_1985_window8.matrix';

#ouput params
my $outputDir = '/home/henryst/clusterVis/output/rayFish_threshold3000_new/';
my $nodeInfoFileOut = 'RayFish_nodeInfo';
my $edgeInfoFileOut = 'RayFish_edgeInfo';
my $cuiInfoFileOut = 'RayFish_cuiInfo';




###########################################
#            Perform Clustering
###########################################
my $command = "perl clustering/createVisulizationOutput.pl $vclusterLocation $lbdFile $vectorFile $clusterMethod $assocMeasure $assocMatrix $outputDir $nodeInfoFileOut $edgeInfoFileOut";
print STDERR "\n\n$command\n\n";
`$command`;


###########################################
#            Get Cui Info
###########################################
#$command = "perl general/getCuiInfo.pl $outputDir$nodeInfoFileOut $outputDir$cuiInfoFileOut";
#print "$command\n";
#`$command`;


print "DONE!\n";



