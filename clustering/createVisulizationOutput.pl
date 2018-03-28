#!/usr/bin/perl
# usage: perl createClusterTree.pl [vclusterLocation] [lbdFile] [vectorFile] [clMethod] [outDir] [parentArrayOut] [clustersFileOut]
# example: TODO update example
#
# Input Parameters:
#  # input params:  
#  lbdFile - the file containing LBD output (target terms)
#  vectorFile - the file containing word2vec vectors
#
#  # clustering params:
#  vclusterLocation - the location of the CLUTO vcluster executable script
#  clMethod - the clustering method to use: rb, rbr, direct, agglo, graph, bagglo
#
#  # scoring params:
#  assocMeasure - string indicating the association measure to use (e.g. x2, ll, etc)
#  cooccurrenceMatrix - the location of the association matrix used in calculating association scores 
#
#  # output params:
#  outputDir - the directory to output clustering results to
#  parentArrayOut - the output fileName of the parent array
#  clustersFileOut - the output fileName of the clusters file
#
# Output:
#  The parent array and the clusters file out. These are formatted as TODO, describe IO
#
#
#TODO - vector file format?  well its word2vec format
#TODO - finalize doucmentation (methods are done, but how to use should be documented more)
use strict;
use warnings;

use UMLS::Association;



Main();


########################################
# Begin Code
########################################

# Main routine to perform clustering, scoring, and outputting
# Input:  none
# Output: none
sub Main {

    #get the input parameters
    my ($vclusterLocation, $lbdFile, $vectorFile, $clMethod, $assocMeasure, $cooccurrenceMatrix, $outputDir, $parentArrayOut, $clustersFileOut) = &getArgs();

    # read in data from file, since all cuis may not have a vector, the ogCuiList may be different from the cuiList
    print STDERR "   Reading Input File\n";
    my ($startingCuis, $ogCuiList, $cuiScores, $cuiTerms) = &readLBDData($lbdFile);
    my ($vectors, $cuiList) = &extractVectors($ogCuiList, $vectorFile);
    my $matrixSize = &getMatrixSize($vectors);

    # run clustering	
    print STDERR "   Clustering\n";
    my $vectorInputFile = &printVClusterInputFile($vectors, $vectorFile, $outputDir, $matrixSize);
    (my $clustersFile) = &runVCluster(
	$outputDir, $parentArrayOut, $vclusterLocation, $vectorInputFile, $clMethod, $matrixSize);

    # cluster labeling
    print STDERR "Labeling Clusters\n";
    my $clusters = &extractClusters($parentArrayOut, $cuiList);
    my $centroids = &calculateCentroids($vectors, $clusters, $cuiList, $matrixSize);
    my $clusterLabels = &labelClusters($centroids, $clusters, $vectors, $cuiTerms, $cuiList);

    # score the clusters
    print STDERR "Calculating Cluster Scores\n";
    my $nodeScores = &calculateScores_setAssociation($clusters, $startingCuis, $assocMeasure, $cooccurrenceMatrix);
    my $edgeScores = &calculateScores_sumOfDescendantScores($nodeScores, $clusters, $parentArrayOut);

    # print out the cluster info (clusterID, clusterLabel, clusterWeigt, edgeWeight, cuisInCluster)
    print STDERR "Outputting Results\n";
    open OUT, ">$clustersFileOut" or die ("ERROR: cannot open outFile: $clustersFileOut\n");
    foreach my $clusterID (keys %{$clusters}) {
	print OUT "$clusterID\t${$clusterLabels}{$clusterID}\t${$nodeScores}{$clusterID}\t${$edgeScores}{$clusterID}\t";
	my $cuisString = '';
	foreach my $cui (@{${$clusters}{$clusterID}}) {
	    $cuisString .= "$cui,";
	}
	chop $cuisString;
	print OUT $cuisString;
	print OUT "\n";
    }

    print STDERR "Done\n";
}



########################################
# Cluster Scoring
########################################
# Calculates the cluster as the sum of all leaf node scores
# Input:  $cuiScores - a hash ref of the LBD scores of each cui ($cuiScores{$cui}=$score)
#         $clusters- a hash ref of clusters and their cuis ($clusters{$id}=@cuis)
# Output: \%clusterScores - a hash ref of scores for each cluster ($clusterScores{$id}=$score)
sub calculateScores_sumOfLeafNodes {
    my ($cuiScores, $clusters) = (@_);	
    my %clusterScores;
    foreach my $clusterIndex (keys %$clusters){
	$clusterScores{$clusterIndex} = 0;
	foreach my $cui (@{$clusters->{$clusterIndex}}){
	    my $score = $cuiScores->{$cui};
	    $clusterScores{$clusterIndex} += $score;
	}
    }
    return \%clusterScores;
}
#Note: No longer used

# Calculates the cluster as the sum of all descendant node scores
# remember, a node score is seperate from the cluster score. Those
# are passed in. In this case, the node score is the set association
# of that cluster and the starting term. This is precalculated and 
# passed in. The cluster score is the sum of all descendant node scores
# ...so the cluster score will be incoming edge width, node score will
# be the size of the displayed node
# Input:  $nodeScoresRef - a hash ref of nodes and their scores ($nodes{$id}=$score)
#         $clusters- a hash ref of clusters and their cuis ($clusters{id}=@cuis)
#         $parentFile - the file containing the parent list of the hierarchy
# Output: \%clusterScores - a hash ref containing scores of each cluster ($clusterScores{$id}=$score)
sub calculateScores_sumOfDescendantNodes {
    my ($nodeScoresRef, $clusters, $parentFile) = (@_);	

    #get a list of all descendant nodes
    my $childrenListRef = &extractChildrenList($parentFile);

    #score each cluster as the sum of its descendant node scores
    my %cluster_scores;
    foreach my $clusterIndex (keys %$clusters){
	#recursively find the score
	$cluster_scores{$clusterIndex}  = &sumDescendantScores($clusterIndex, $nodeScoresRef, $childrenListRef);
    }
    return \%cluster_scores;
}

# Recursively sums all descendant scores of a node
# Input:  $nodeID - the ID of the current node (in recursion) 
#         $nodeScoresRef - a hash ref containing the scores of each node (%nodeScores{id}=score)
#         $childrenListRef - a hash ref containing the children of each node (%children{id}=\@childNodeIDs)
# Output: $score - the sum of all descendants node values of this node
sub sumDescendantScores {
    my $nodeID = shift;
    my $nodeScoresRef = shift;
    my $childrenListRef = shift;

    #sum the score of all descendants recursively
    my $score = 0;
    foreach my $childrenArray (keys ${$childrenListRef}{$nodeID}) {
	foreach my $childID (@{$childrenArray}) {
	    $score += ${$nodeScoresRef}{$childID};
	    $score += &sumDescendantScores($childID, $nodeScoresRef, $childrenListRef);
	}
    }

    return $score;
}


# calculates the score of each node as the set association between it and the starting set
# Input:  $clusters - a hash ref of clusters and their cuis (clusters{id}=@cuis)
#         $aTermsRef - an array ref of CUIs that are the A terms in LBD
#         $measure - the association measure to use (e.g. x2, ll, etc..)
# Output: \%nodeScores - a hash ref of scores for each node (%nodeScores{$id}=$score)
sub calculateScores_setAssociation {
    my $clusters = shift;
    my $aTermsRef = shift;
    my $measure = shift;
    my $cooccurrenceMatrix = shift;
    
    #initialize UMLS::Association
    my %tHash = ();
    $tHash{'t'} = 1; #default hash values are with t=1 (silence module output)
    $tHash{'vsa'} = 1; #always use overlapping B sets association, since it is implicit
    $tHash{'matrix'} = $cooccurrenceMatrix;
    $tHash{'noOrder'} = 1;
    my $componentOptions = \%tHash;
    my $association = UMLS::Association->new($componentOptions) or 
	die "Error: Unable to create UMLS::Association object.\n";
    #TODO - options{'uml'} = interface; (is that needed?)


    #create set pair lists for a to all nodes in the dataset
    # the set pair list is what is passed to UMLS::Association
    # it is two parallel arrays, the A terms, and each node set
    # UMLS::Association then returns the the set associations between
    # the A to C pair at each index
    my @cuis1 = ();
    my @cuis2 = ();
    my @order = ();
    foreach my $clusterKey (keys %{$clusters}) {
	my $cTermsRef = ${$clusters}{$clusterKey};
	push @cuis1, $aTermsRef;
	push @cuis2, $cTermsRef;
	push @order, $clusterKey;
    }

    #calculate the assocaitions
    my $scoresRef = $association->calculateAssociation_setPairList(\@cuis1,\@cuis2,$measure);

    #assign each node its calculated score
    my %nodeScores = ();
    for (my $i = 0; $i < scalar @order; $i++) {
	$nodeScores{$order[$i]} = ${$scoresRef}[$i];
    }

    #return the node scores
    return \%nodeScores;
}




########################################
# Cluster Extracting
########################################
# Converts the parent file into a children list
# the children list is a hash, where the key is a cluster id, and the
# value is an array ref to a list of cluster ids that are children
# Input:  $parentFile - the file containing the parent list of the hierarchy
# Output: \%children - a hash containing the children of each node (%children{id}=\@childNodeIDs)
sub extractChildrenList {
    my $parentFile = shift;

    #read the parent tree file
    open IN, $parentFile or die ("ERROR: cannot open parentFile: $parentFile\n");
    my %parents = ();
    my $index = 0;
    while (my $line = <IN>) {
	my ($clusterNum, $val1, $val2) = split(/ /, $line);
	$parents{$index} = $clusterNum;
	$index++;
    }
    close IN;

    #convert the parent tree into a children list
    #initialize the list
    my %children = ();
    foreach my $id (keys %parents) {
	my @emptyArray = ();
	$children{$id} = \@emptyArray;
    }
    #find the children of each id (cluster number ID)
    foreach my $id1 (keys %parents) {
	#find all children of $id1
	foreach my $id2 (keys %parents) {
	    #see if the parent of $id2 is $id1
	    if ($parents{$id2} == $id1) {
		#add $id2 to the list of children for $id1
		push @{$children{$id1}}, $id2;
	    }
	}
    }

    #return the children hash
    return \%children;
}

# Extracts the clusters and the cuis of them from files
# Input:  $parentFile - the file containing the parent list of the hierarchy
#         $cuiList - hash ref of cui indeces (cuiList{cui}=index)
# Output: \%clusters - a hash ref of clusters and their cuis (clusters{id}=@cuis)    
sub extractClusters {
    my $parentFile = shift;
    my $cuiList = shift;
    
    #get the children list from the tree file
    my $childrenListRef = &extractChildrenList($parentFile);

    #convert the children list to a list of concepts in each cluster
    my %clusters = ();
    foreach my $id (keys %{$childrenListRef}) {
	$clusters{$id} = &getChildCuiIDs($childrenListRef, $id);
    }

    #convert the CUI IDs to CUIs
    foreach my $id (keys %clusters) {
	my @cuis = ();
	foreach my $cuiID (@{$clusters{$id}}) {
	    push @cuis, ${$cuiList}{$cuiID};
	}
	$clusters{$id} = \@cuis;
    }

    #return the hash of clusters and their CUIs
    return \%clusters;
}


#Recursively gets the descendent leaf nodes (CUIs)
# if the children tree. Returns an array ref of CUI IDs
# Input:  $children - a hash ref containing the children of each node (%children{id}=\@childNodeIDs)
#         $id - the ID of the node to get all children of
# Output: \@cuis - an array ref of cuis that are leaf nodes of this hierarchy
sub getChildCuiIDs {
    my $childrenRef = shift;
    my $id = shift;

    #intialize
    my @children = @{${$childrenRef}{$id}}; 
    my @cuis = ();

    #base case, this is a leaf node (no children)
    if (scalar @children == 0) {
	push @cuis, $id;
	return \@cuis;
    }

    #else keep going down until a leaf node is hit
    foreach my $childID(@children) {
	push @cuis, @{&getChildCuiIDs($childrenRef, $childID)};
    }
    return \@cuis;
}




########################################
# Cluster Labeling
########################################
# Calculates the "centroid" (average value of all vectors) for all 
# clusters for a particular clustering solution
# Input:  $vectors - a hash ref of vectors (vectors{id}=@vector)
#         $clusters - a hash ref of clusters (clusters{id}=@cuis)
#         $cuiList - hash ref of cui indeces (cuiList{cui}=index)
# Output: \%centroids - a hash ref of cluster centroids (centroid{$id}=@vector)
sub calculateCentroids {
    my ($vectors, $clusters, $cuiList, $matrixSize) = (@_);
    my %centroids;
    my %cuiList = reverse %{$cuiList}; #makes the cui list a hash{cui}=index rather than hash{key}=index
    my $numColumns = @$matrixSize[1];

    #find the centroid of each cluster
    foreach my $clusterIndex (keys %$clusters){
	#initialize array of size num_columns to 0
	$centroids{$clusterIndex} = [(0) x $numColumns]; 

	#sum vectors of each cui in the cluster
	my $vectorsPerCluster = 0;
	foreach my $cui (@{$clusters->{$clusterIndex}}){

	    #get this cui's vector
	    my $index = $cuiList{$cui};
	    my @vector = @{$vectors->{$index}}; # get vector values from vectors hash
	    
	    #component-wise summation of the vector
	    foreach my $col (0 .. $numColumns-1){
		$centroids{$clusterIndex}[$col] += $vector[$col];
	    }
	    $vectorsPerCluster++;
	}
	#component-wise average over summation
	foreach my $col (0.. $numColumns-1){
	    $centroids{$clusterIndex}[$col] /= $vectorsPerCluster;
	}
    }
    return \%centroids;
}


# Labels cluster centroid with vector closest to the centroid, 
# determined with cosine similarity 
# Input:  $centroids - a hash ref of centroids (centroid{$id}=@vector)
#         $clusters - a hash ref of clusters (clusters{id}=@cuis)
#         $vectors - a hash ref of vectors (vectors{id}=@vector)
#         $cuiTerms - a hash ref of cuiTerms (cuiTerms{cui}=term)
#         $cuiList - hash ref of cui indeces (cuiList{cui}=index)
# Output: \%clusterNames - a hash ref of the names of each cluster (clusterNames{ID}=name)
sub labelClusters {
    my ($centroids, $clusters, $vectors, $cuiTerms, $cuiList) = (@_);
   
    my %clusterNames = ();;
    my %cuiList = reverse %$cuiList;
    foreach my $clusterIndex (keys %$centroids){
	my $centroidVectorIndex;
	my $maxCosSim = -1; # start with minimum value for cosine similarity
	my @centroid = @{$centroids->{$clusterIndex}};
	foreach my $cui (@{$clusters->{$clusterIndex}}){
	    my $vectorIndex = $cuiList{$cui};
	    my @vector = @{$vectors->{$vectorIndex}};
	    my $cosSim = &calculateCosSim(\@centroid, \@vector);
	    if ($cosSim > $maxCosSim){
		$maxCosSim = $cosSim; # update maximum cosine similarity value
		$centroidVectorIndex = $vectorIndex;	# update centroid vector index
	    }
	}
	my $cui = $cuiList->{$centroidVectorIndex};
	my $term = $cuiTerms->{$cui};
	$clusterNames{$clusterIndex} = $term;	
    }
    return \%clusterNames;
}

# Calculates cosine similarity between two vectors (aka the dot product)
# cosine similarity forumla (Vector a, Vector b):
# (a * b) / (|a|*|b|)
# Input:  $a - aVector = an array of values
#         $b - bVector = an array of values 
# Output: $cosSim - the cosine similarity between $a and $b
sub calculateCosSim {
    my ($a, $b) = (@_);
    my @a = @$a;
    my @b = @$b;
    my $numerator;
    my $a_sq_sum;
    my $b_sq_sum;

    foreach my $i (0 .. $#a){
	$numerator += $a[$i]*$b[$i];
	$a_sq_sum += $a**2;
	$b_sq_sum += $b**2;
    }
    my $denominator = sqrt($a_sq_sum)*sqrt($b_sq_sum);
    my $cosSim = $numerator/$denominator;
    return $cosSim;
}


########################################
# Read Input Files
########################################

# Reads in the LBD target term fileName
# Input:  $lbdFile - the file containing LBD output (lbd target term file)
# Output: \@startingCuis - an array ref of statrting cuis (@startingCuis[i]=$cui)
#         \%cuiList - a hash ref of cui indeces (%cuiList{$cui}=$index)
#         \%cuiScores - a hash ref of target term scores ($cuiScores{$cui}=$score)
#         \%cuiTerms - a hash ref of preferred terms ($cuiTerms{$cui}=$term)
sub readLBDData {
    my ($lbdFile) = (@_);
    my @startingCuis = ();
    my %cuiList;
    my %cuiScores;
    my %cuiTerms;
    open my $fh, '<', "$lbdFile" or die "Can't open $lbdFile: $!";
    while (my $line = <$fh>) {
	#grab starting cuis
	if ($line =~ /startCuis -> /) {
	    while ($line =~ /(C\d{7})/g) {
		push @startingCuis, $1;
	    }
	}

	#grab target cui info
	if ($line =~ /(\d+)\t(\d+.?\d*)\t(C\d{7})\t(.+)/){
	    my $index = $1 - 1;
	    my $score = $2;
	    my $cui = $3;
	    my $term = $4;
	    $cuiList{$cui} = $index;
	    $cuiScores{$cui} = $score;
	    $cuiTerms{$cui} = $term;
	}
    }
    close $fh;
    if (keys %cuiList == 0){
	print "Invalid data in ltc file: $lbdFile.\n";
	exit;
    }
    return \@startingCuis, \%cuiList, \%cuiScores, \%cuiTerms;
}
#NOTE: if we don't want to use LBD output, then we can get terms from elsewhere, and we don't actually use the individual cui scores anymore


#Reads the vector file
# Input:  $cuiList - a hash ref of cui indeces (%cuiList{$cui}=$index)
#         $vectorFile - file containing vector representations of terms
# Output: \%vectors - a hash ref of vectors (vectors{id}=@vector)
#         \%newCuiList - the updated cui list (cuis that don't have a vector 
#                        have been removed) (%newCuiList{$cui}=$index)
sub extractVectors {
    my ($cuiList, $vectorFile) = (@_);
    my %vectors; # later, sort vectors by descending rank and print
    my %newCuiList;
    my $index = 0;
    open my $fh, '<', "$vectorFile" or die "Can't open $vectorFile: $!";
    while (my $line = <$fh>){
	if ($line =~ /^(C\d{7})(.+)/){ # if line contains a cui/vector pair
	    my $cui = $1;		
	    if (exists $cuiList->{$cui}){ # if exists in target terms list (not all target terms will be in the vector list due to w2v threshold)
		my $vector = $2;
		my @vectorVals = ($vector =~ /(-?\d.\d+)/g);
		$vectors{$index} = [@vectorVals];
		$newCuiList{$index} = $cui;
		$index += 1;
	    }
	}
    }
    close $fh;
    return \%vectors, \%newCuiList;
}


# Gets the size of a matrix (stored as a hash of hashes)
# Input:  $vectors - a hash ref of vectors (vectors{id}=@vector)
# Output: $matrixSize - array ref, containing: rows,cols
sub getMatrixSize {
    my ($vectors) = (@_);
    my @matrix_size;
    $matrix_size[0] = keys %$vectors;	# number of rows
    $matrix_size[1] = @{$vectors->{0}};	# number of columns
    return \@matrix_size;
}




########################################
# VCluster Prep and Running
########################################

# Converts internal format to VCluster format so it can be run
# Input:  $vectors - a hash ref of vectors (vectors{id}=@vector)
#         $vector_file - the file containing vectors
#         $cl_soln_dir - the directory to output the clustering solution
#         $matrix_size - array ref of the dimensions of the vector file (rows, cols)
# Output: $v_ifile - the file formatted for input into CLUTO
sub printVClusterInputFile {
    # print new vector file with _v appended
    my ($vectors, $vector_file, $cl_soln_dir, $matrix_size) = (@_);

    my $num_rows = @$matrix_size[0];
    my $num_columns = @$matrix_size[1];

    $vector_file  =~ s/.*\/(.*?)/$1/;
    my $v_ifile = $vector_file . '_v';
    my $v_ifile_path = $cl_soln_dir . '/' . $v_ifile;

    if (-e $v_ifile_path && -f $v_ifile_path){
	return $v_ifile;
    }

    open my $fh, '>', "$v_ifile_path" or die "Can't open $v_ifile_path: $!"; # open file in create/write/truncate mode
    print $fh "$num_rows $num_columns\n";
    foreach my $vector_i (sort {$a <=> $b} keys %$vectors){
	print $fh join (' ', @{$vectors->{$vector_i}});
	print $fh "\n";
    }
    close $fh;

    return $v_ifile;
}

#runs Vcluster and outputs the clusters and the clusters tree
# Input:  $output_dir - the directory to output results
#         $tree_filename - the file containing the parent array (will be output to)
#         $vclusterLocation - the location of the vcluster script
#         $v_ifile - the vector file formatted for input into CLUTO
#         $cl_method - the clustering method (e.g. agglo, bagglo, rb, rbr, etc)
#         $matrix_size - array ref of the dimensions of the vector file (rows, cols)
# Output: $clusters_filename - the file name the cluserting solution was output to
sub runVCluster {
    my ($output_dir, $tree_filename, $vclusterLocation, $v_ifile, $cl_method, $matrix_size) = (@_);

    #We want a cluster tree where each leaf node is a vector
    # thereby creating a full cluster tree. SO num_clusters = num_vectors
    my $num_vectors = @$matrix_size[0];
    my $num_clusters = $num_vectors;

    #create file names
    my $file_to_cluster = $output_dir.$v_ifile;
    my $clusters_filename = $file_to_cluster.'.clusters';
    #my $tree_filename =  $file_to_cluster.'.tree';
    my $labels_filename = $file_to_cluster.'labels';

    #NOTE: the labeltree option labels with a set of features, not a single name
    #create command and perform clustering
    #my $cmd = "./$vcluster_dir/vcluster $file_to_cluster $num_clusters -clmethod=$cl_method -clustfile=$clusters_filename -showtree -cltreefile=$tree_filename";
 my $cmd = "./$vclusterLocation $file_to_cluster $num_clusters -clmethod=$cl_method -clustfile=$clusters_filename -fulltree -treefile=$tree_filename";
    print "$cmd\n";
    my $cluster_cmd = `$cmd`;

    return $clusters_filename;
}





########################################
# Get and Check Input Parameters
########################################
# gets and checks the input arguments to the program
# Input:  none
# Output: the input arguments
sub getArgs {
    #get the input arguments
    my $vclusterLocation = $ARGV[0];
    my $lbdFile = $ARGV[1];
    my $vectorFile = $ARGV[2];
    my $clMethod = $ARGV[3];
    my $assocMeasure = $ARGV[4];
    my $cooccurrenceMatrix = $ARGV[5];
    my $outputDir = $ARGV[6];
    my $parentArrayOut = $ARGV[7];
    my $clustersFileOut = $ARGV[8];
  
    #check the input args
    &checkNumArgs();
    &checkVcluster($vclusterLocation);
    &checkCLMethod($clMethod);
    &checkAssocMeasure($assocMeasure);
    &checkFileErr($cooccurrenceMatrix);
    &checkFileErr($lbdFile);
    &checkFileErr($vectorFile);
    &createDir($outputDir);
    $parentArrayOut = $outputDir.$parentArrayOut;
    &checkCreateFileErr($parentArrayOut);
    $clustersFileOut = $outputDir.$clustersFileOut;
    &checkCreateFileErr($clustersFileOut);

    #return the input arguments
    return $vclusterLocation, $lbdFile, $vectorFile, $clMethod, $assocMeasure, $cooccurrenceMatrix, $outputDir, $parentArrayOut, $clustersFileOut;
}


# Checks the number of args passed into the program
# Input:  none
# Output: none
sub checkNumArgs {
    if ($#ARGV != 8) {
	print "Incorrect number of arguments. Program requires 9 arguments to run.\n";
	print "Usage: perl DiscoveryReplication.pl [vclusterLocation] [lbdFile] [vectorFile] [clMethod] [assocMeasure] [cooccurrenceMatrix] [outDir] [parentArrayOut] [clustersFileOut]\n";
	exit;
    }
}

# Checks if the vclusterLocation parameter is correct
# Input:  $vclusterLocation - the location of the vClusterScript
# Output: none
sub checkVcluster {
    my ($vclusterLocation) = (@_);
    if (! -e "$vclusterLocation" || ! -e $vclusterLocation){
	print "Incorrect location for vcluster program: $vclusterLocation.\n";
	exit;
    }
}

# Checks if the clMethod parameter is correct
# Input:  $clMethod - the string specifying the clustering method
# Output: none
sub checkCLMethod {
    my ($clMethod) = (@_);
    my @validCLMethods = ('rb', 'rbr', 'direct', 'agglo', 'graph', 'bagglo');
    if (!(grep /$clMethod/, @validCLMethods)) {
	print "Invalid method. Supported options are 'rb', 'rbr', 'direct', 'agglo', 'graph', 'bagglo'.\n";
	print "See CLUTO manual for descriptions of each clustering method (-clmethod).\n";
	exit;
    }
}

# Checks if the associasion measure is correctly specified
# Input: $assocMeasure - string specifying the association measure
# Output: none
sub checkAssocMeasure {
    my $assocMeasure = shift;
    my @validAssocMeasures = ('freq', 'dice', 'leftFisher', 'rightFisher', 'twoTailed', 'jaccard', 'll', 'tmi', 'odds', 'pmi', 'phi', 'x2', 'ps', 'tscore');
    if (!(grep /$assocMeasure/, @validAssocMeasures)) {
	print "Invalid association measure. Supported options are ".join(', ',@validAssocMeasures)."\n";
	print "See UMLS::Association --help for information\n";
	exit;
    }

}

# Checks if input files are correct
# Input:  $file - the file to check if it exists and is a file
# Output: none
sub checkFileErr {
    my ($file) = (@_);
    if (! -e $file || ! -f $file){
	print "File not found: $file.\n";
	exit;
    }
}

# Checks if output files are correct
# Input:  $file - the file to create
# Output: none
sub checkCreateFileErr {
    my ($file) = (@_);
    open OUT, ">$file" or die("ERROR: cannot create file: $file\n");
    close OUT;
}

# Creates a directory
# Input:  $dir - the directory to create
# Output: none
sub createDir {
    my ($dir) = (@_);
    # check if dir exists, and if its a dir.
    my $made = (-e $dir || -d $dir);
    
    #create the dir if needed
    if (!$made) {
	$made = mkdir $dir;
    }
    ($made) or die ("ERROR: directory does not exists, and cannot create directory: $dir\n");
}


