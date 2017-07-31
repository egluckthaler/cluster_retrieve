#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FileHandle;
$| = 1;
sub usage {
	my $message = shift;
	if (not defined($message)) { $message = ''}
	my $usage = qq/
Usage: cluster_retrieve.pl [OPTIONS]
-g <gff_file> filename must be formatted: <omecode>.rest_of_filename
-c <coo3_file> Required if not providing a gff file (can contain multiple genomes)
-p <proteome> multifasta (can contain multiple genomes, if providing coo3)
-f <cluster_model_fasta> 
-i <cluster_model_info> 
-o <analysis_name> should include absolute path

This script searches your genome of interest for clusters containing
hits to proteins supplied in the cluster fasta file, and outputs a 
two files with: 
	1) clusters containing anchor genes
	2) clusters not containing anchor genes. 
gff files must be in the standard format adhered to by JGI, where the following collumns contain:
	1st column = contig ID
	3rd column = feature type (e.g., CDS)
	4th column = most downstream position of CDS
	5th column = most upstream position of CDS
	7th column = strand (+ or -)
	8th column = annotation info
e.g., VVO_00001 JGI CDS 1562 1588 . - 0 name "yvh1"; proteinId 110841; exonNumber 8
Headers in the multifasta file must be formatted as 
	>genomecode_proteinID, where proteinIDs have no '_'s
If supplying a coo3 file, it must be a whitespace separated file where:
	1st column = sequenceheader (ie., genomecode_proteinID, where proteinIDs have no '_'s)
	2nd column = contig name
	3rd column = most downstream position of CDS
	4th column = most upstream position
	5th column = strand (+ or -)
Criteria for hits are (can be manually changed starting at line 153 in this script):
	1) minimum 50 bitscore
	2) minimum 30% identity 
	3) alignment is 50-150% of query sequence length 
Criteria for defining clusters based on:
	1) a maximum of 6 intervening genes between hits
A coo3 file and protein multifasta, with information from multiple genomes may be provided 
to search multiple genomes at once for a cluster of interest\n\n/;
	die($usage, $message);
}

main: {

	#check if correct switches were supplied
	my %opts;
	getopt('g:c:p:f:i:o:h:', \%opts);
	Opts_check(\%opts);

	#Check if name on info file and fasta file corespond
	my ($infoname)= fileparse($opts{'i'}, ".info");
	my ($fastaname) = fileparse($opts{'f'}, ".fasta");
	usage("Error: the provided fasta file must correspond to the correct info file.. do their names match?\n") if ($infoname ne $fastaname);
	
	#capture the proteome name from proteome file
	my $proteome = basename($opts{'p'});

	#capture the directory to output files from the outputfile supplied by user
	my ($outfilename,$outpath) = fileparse($opts{'o'});
	my $blastpath = "$outpath/blast";
	
	#create a blastdb of proteome of interest
	if (! -e "$blastpath/$proteome.pdb.psd") {
		print "Now creating BLAST database of $proteome..\n";
		system("mkdir $blastpath");
		my $fail_check = system("makeblastdb -in $opts{'p'} -dbtype prot -parse_seqids -input_type fasta -out $blastpath/$proteome.pdb");
		if ($fail_check != 0) { die "Error: could not construct BLAST database.. Exiting..\n$!\n"}
	} else {
		print "BLAST database of $proteome already exists.. skipping database creation step..\n";
	}

	#query proteome blastdb for all cluster protein queries
	if (! -e "$blastpath/$outfilename.blast.parsed" || -z "$blastpath/$outfilename.blast.parsed") {
		print "Executing BLASTp to search for all best hits to queries in $fastaname..\n";
		Blastp($opts{'f'}, $blastpath, "$blastpath/$proteome.pdb", $outfilename); #return $outfilename.blast.ids, $outfilename.blast.out, $outfilename.blast.parsed in $outpath
	} else {
		print "Best hits to queries in $fastaname have already been found.. skipping BLASTp step..\n";
	}
	#Prepare coo3 file if it is not provided capture the directory to gff file
	my $coo3file;
	$coo3file = $opts{'c'} if (defined $opts{'c'});
	if (not defined $coo3file) {
		my ($gfffilename, $gffpath) = fileparse($opts{'g'});
		my $coopath = "$gffpath/coo3";
		system("mkdir $coopath") if (! -d $coopath);
		$coo3file = "$coopath/$gfffilename.coo3";
		print "Creating coo3 from supplied gff3 file\n" if (! -e $coo3file);
		gff2coo($opts{'g'}, $coo3file) if (! -e $coo3file);
	} else {
		print "A coo3 file has been provided.. skipping coo3 creation step..\n";
	}

	#search for clusters using blast.ids and user-supplied gff file that is transformed into .coo3 format
	if (! -e "$outpath/$opts{'o'}.anchor" || -z "$outpath/$opts{'o'}.anchor") {
		
		#load information about query cluster proteins into a hash
		open(my $info_in, '<', $opts{'i'}) or usage("Error: cannot read cluster info file.. Exiting..\n");
		my $inforef = Info_parse($info_in);
	
		#load all of the hits to clustered proteins in proteome of interest into a hash
		open(my $hits_in, '<', "$blastpath/$outfilename.blast.parsed") or usage("Error: cannot open $blastpath/$outfilename.blast.parsed for reading\n");
		my ($hits) = Blast_parse($hits_in);
		#search through gff of interest to find information about hits
		
		print "Now finding clusters in $proteome..\n";
		my ($hit_info) = Coo3_parse($coo3file, $hits);

		#see if protids are physically grouped into clusters, remove any that do not cluster with at least one other protein
		#tempclref is structured as clusters{"{ome}_c{cl_count}"}{prot_id} = [info]
		my ($tempclref, $singletons) = Clusterfy($hit_info);

		#check to see which clusters contain anchor gene of interest (based on MCL1000), remove them from %tempclusters and place into %trueclusters
		my ($trueclref) = Anchor_clusters($tempclref, $hits);
	
		#annotate clustered proteins with information from the queries that fetched them
		Annotate_clusters($tempclref, $inforef);
		Annotate_clusters($trueclref, $inforef);

		#print out clusters to file	
		print "Found ".scalar(keys(%{$singletons}))." unclustered hits\n";	
		open(my $noan_out, '>', "$opts{'o'}.no_anchor") or die "Error: cannot open file for output.. Exiting..\n";
		print "Found ".scalar(keys(%{$tempclref}))." clusters missing anchor\n";
		Cluster_print($tempclref, $noan_out);
		open(my $an_out, '>', "$opts{'o'}.anchor") or die "Error: cannot open file for output.. Exiting..\n";
		print "Found ".scalar(keys(%{$trueclref}))." clusters with anchor\n";
		Cluster_print($trueclref, $an_out);
	
	} else {
		print "Outputfile $opts{'o'}.anchor already exists.. Skipping outputfile creation step..\n";
	}
}


sub Opts_check {
	
	my ($opts) = @_;
	
	#check if all switches are defined
	usage() if (exists $opts->{'h'});
	usage("Error: please supply a gff3 file or a coo3 file\n") if (! defined $opts->{'g'} && ! defined $opts->{'c'});
	usage("Error: no proteome file supplied\n") if (! defined $opts->{'p'});
	usage("Error: no cluster fasta supplied\n") if (! defined $opts->{'f'});
	usage("Error: no cluster info supplied\n") if (! defined $opts->{'i'});
	usage("Error: no output file supplied\n") if (! defined $opts->{'o'});

	#check if files are accesible
	if (defined $opts->{'g'} && ! -r $opts->{'g'}) { usage("Error: cannot read gff3 file.. Exiting..\n") }
	if (defined $opts->{'c'} && ! -r $opts->{'c'}) { usage("Error: cannot read coo3 file.. Exiting..\n") }
	if (! -f $opts->{'p'}) { usage("Error: cannot locate proteome file.. Exiting..\n") }
	if (! -r $opts->{'f'}) { usage("Error: cannot read cluster fasta file.. Exiting..\n") }
	if (! -r $opts->{'i'}) { usage("Error: cannot read cluster info file.. Exiting..\n") }
}


sub Blastp {
	my ($infile, $path, $db, $outfilename) = (@_);
	my %ids;
	my $blastfile = "$path/$outfilename.blast.all";
	my $parsedblast = "$path/$outfilename.blast.parsed";
	my $idfile = "$path/$outfilename.blast.ids";

	my $fail_check = system("blastp -db $db -query $infile -evalue 1e-4 -outfmt '6 std qlen' -num_alignments 10000 -num_descriptions 10000 -num_threads 1 > $blastfile");
	if ($fail_check != 0) { die "Error: could not conduct BLASTp search.. Exiting..\n$!\n"}

	open(my $blast_in, '<', $blastfile);
	open (my $new, ">", "$parsedblast");
	$new-> autoflush(1);
	while (my $line = <$blast_in> ){
		chomp $line;
		my @fields = split ("\t", $line);
		my ($subjectid, $percentid, $bitscore, $alnlen, $qlen) = ($fields[1], $fields[2], $fields[11], $fields[3], $fields[12]); #grab alignment length, percent positive matches, and query sequence length
		if ($percentid >=30) { #change min perc id here
			if ($bitscore >= 50) { #change min bit score here
				if (($alnlen/$qlen) >= 0.5 && ($alnlen/$qlen) <= 1.5) {
					$ids{$subjectid} = 1; #this should remove any duplicates as well: not in the actual blastoutput file, but in the id file that is printed out
					print $new "$line\n";
				}
			}
		}
	}
	open (my $out1, ">", "$idfile") or die;
	$out1-> autoflush(1);
	print $out1 "$_\n" foreach (sort (keys %ids));
}

sub Clusterfy {
	my $usage = qq/
Usage: Clusterfy(<\%protein_positional_info>)
sub Clusterfy takes a hash of hash{ome}{contig}{position}{protein} = [info array]
as input (e.g., produced by Coo3_parse Ebinf::Utils and calculates distances between 
proteins on each contig. Returns a hash containing all proteins separated by less than 
6 intervening genes, grouped by cluster. 
cluster hash structure: clusters{"{ome}_c{cl_count}"}{prot_id} = [info]
loner hash structure: loners{"{ome}_ffalsecount"}{prot_id}
/;
	if (scalar @_ != 1) { die ($usage) }
	
	my ($hitref) = @_;	
	my %hit_info = %{$hitref};
	my (%clusters, %loners);
	foreach my $ome (sort keys %hit_info) {
		my ($cl_count, $falsecount) = (1, 1);
		foreach my $contig (sort keys %{$hit_info{$ome}}) {
			my @pos = sort {$a <=> $b} keys %{$hit_info{$ome}{$contig}}; #sort all protein positions from smallest to largest on each contig
			for (my $i = 0; $i < scalar @pos; $i++) { #now work through all positions on contig associated with a protein of interest
				my ($prot_id) = keys %{$hit_info{$ome}{$contig}{$pos[$i]}};
				$clusters{"${ome}_c${cl_count}"}{$prot_id} = $hit_info{$ome}{$contig}{$pos[$i]}{$prot_id}; #add protein to cluster, because it definitely is not more than 6 intervening genes away from the nearest preceeding (upstream) protein
				#now check status of the next protein downstream
				#there are 2 conditions to define the downstream boundary of a cluster: if neither are met, keep adding proteins to cluster
				if (! defined $pos[$i+1]) { #this will happen to the last element of @pos
					if (scalar keys %{$clusters{"${ome}_c${cl_count}"}} <= 1) { #clusters of size 1 excluded
						$loners{"${ome}_f$falsecount"}{$prot_id} = $hit_info{$ome}{$contig}{$pos[$i]}{$prot_id};
						$falsecount++;
						delete $clusters{"${ome}_c${cl_count}"};
					} else { #then the downstream boundary of a true cluster has been found; increment cluster count
						$cl_count++;
					}
				} elsif (($pos[$i+1] - $pos[$i]) > 7) { #is the next closest protein downstream more than 6 intervening genes away?
					if (scalar keys %{$clusters{"${ome}_c${cl_count}"}} <= 1) { #clusters of size 1 excluded
						$loners{"${ome}_f$falsecount"}{$prot_id} = $hit_info{$ome}{$contig}{$pos[$i]}{$prot_id};
						$falsecount++;
						delete $clusters{"${ome}_c${cl_count}"};
					} else { #then the downstream boundary of a true cluster has been found; increment cluster count
						$cl_count++;
					}
				}
			}
		}
	}
	return(\%clusters, \%loners);
}

sub Blast_parse {
	my $usage = qq/
	Usage: Blast_parse(<blast_FH>)
	sub Blast_parse takes as input a filehandle to a tabular BLAST 
	output file and returns a hash where {hit => query}.
	/;
	if (scalar @_ != 1) { die($usage) }
	
	my ($hits_in) = @_;
	my %hits;
	while (my $line = <$hits_in>) {
		chomp $line;
		my @column = split/\t/, $line;
		my ($hit, $query) = ($column[1], $column[0]);
		push @{$hits{$hit}}, $query;
	}
	return(\%hits);
}


sub Coo3_parse {
	my $usage = qq/
Usage: coo3_parse(<path_to_coo3>, \%proteins_of_interest)
sub coo3_parse takes as input a filehandle for a coo3_file and a hash of proteins of 
interest and returns a hash where hash{ome}{contig}{position}{protein} = [info from coo3]
ome is determined by everything up to last _
/;
	if (scalar @_ != 2) { die($usage) }
	
	my ($coo3_path, $protref) = @_;
	my %hits = %{$protref};
	
	open(my $coo3_in, '<', $coo3_path) or usage("Error: cannot read coo3 file from $coo3_path\n");
	
	my $line_counter = 1;
	my %hit_info;
	while (my $line = <$coo3_in>) {
		chomp $line;
		$line =~ s/ +/\t/g;
		my @column = split/\t/, $line;
		if (! defined $column[2]) { usage("Error: coo3 file is not whitespace separated.. Exiting..\n")}
		my ($prot_id, $contig, $position, $strand) = ($column[0], $column[1], "$column[2]-$column[3]", $column[4]);
		if (exists $hits{$prot_id}) { #if the protein is a besthit for a cluster protein
			$prot_id =~ m/^(.+?)_[^_]+$/; #grab everything up to last _
			my $ome = $1;
			push @{$hit_info{$ome}{$contig}{$line_counter}{$prot_id}}, $contig, $position, $strand, @{$hits{$prot_id}};
		}
		$line_counter++;
	}
	return(\%hit_info); 
}

sub Info_parse {
	my $usage = qq/
	Usage: Info_parse(<cluster_info_FH>)
	This sub inputs a filehandle to read a cluster info file
	that contains information about the proteins supplied in
	the cluster fasta file.
	/;
	if (scalar @_ != 1) { die($usage) }
	
	my ($info_in) = @_;
	my %info;
	while (my $line = <$info_in>) {
		chomp $line;
		my @column = split/\t/, $line;
		my ($og, $ann) = ($column[1], $column[2]);
		push @{$info{$og}}, $ann;
	}
	return(\%info);
}


sub Anchor_clusters {
	my $usage = qq/
	Usage: Anchor_clusters(<\%temporary_cluster>, <\%hits_queries>)
	sub Anchor_clusters takes as input a hash of temporary clusters where 
	temp{count}{protid} = [position, strand, query_hit]
	and a hash of BLAST hits where 
	{hit => [queries]}
	and returns a hash of clusters
	containing hits to the query anchor gene
	/;
	if (scalar @_ != 2) { die($usage) }

	my ($tempclusters, $hitref) = @_;
	my %hits = %{$hitref};
	my %trueclusters;
	my %truecount;
	foreach my $temp (sort {$a cmp $b} keys %{$tempclusters}) {
		$temp =~ m/^(.+?)_[^_]+$/;
		my $ome = $1;
		my $anchorfound = 0;
		foreach my $prot_id (sort {$a cmp $b} keys %{$tempclusters->{$temp}}) {
			foreach my $query (@{$hits{$prot_id}}) {
				if ( $query =~ m/MCL1000/) {
					$anchorfound = 1;
				}
			}
		}
		if ($anchorfound == 1) {
			$truecount{$ome}++;
			$trueclusters{"${ome}_ac$truecount{$ome}"} = $tempclusters->{$temp};
			delete $tempclusters->{$temp};
		}
	}
	return(\%trueclusters);
}

sub Annotate_clusters {
	my $usage = qq/
	Usage: Annotate_clusters(<\%clusters, <\%cluster_info)
	This sub takes in refs to cluster hashes and modifies the values
	to contain annotation information stored in \%cluster_info. Modifies
	the ref directly, does not return a new hash.
	/;
	die $usage unless scalar @_ == 2;
	my ($clusters, $info) = @_;

	foreach my $cluster (sort keys %{$clusters}) {
		foreach my $prot_id (sort keys %{$clusters->{$cluster}}) {
 			foreach my $element (@{$clusters->{$cluster}->{$prot_id}}) {
				if (exists $info->{$element}) { #push information about the query MCLgroup associated with the cluster_protein onto the cluster_protein array 
					push @{$clusters->{$cluster}->{$prot_id}}, @{$info->{$element}}; #this will just push the whole array on, so when a cluster protein hits multiple queries, this array will have info about those queries too
				}
			}	
		}
	
	}

}

sub Cluster_print {
	my $usage = qq/
	Usage: Cluster_print (<\%clusters>, <output_fh>)
	This sub prints out all of the information in \%clusters
	to the output filehandle.
	/;
	die $usage unless @_ == 2;
	
	my ($clusters, $out) = @_;
	
	foreach my $cluster (sort keys %{$clusters}) { 
		my %sortedpositions;
		foreach my $protein (sort keys %{$clusters->{$cluster}}) { 
			my $most_downstream = ${$clusters->{$cluster}->{$protein}}[1];
			$sortedpositions{$most_downstream} = "$cluster\t$protein\t".join("\t", @{$clusters->{$cluster}->{$protein}})."\n";
		}
		foreach my $position (sort keys %sortedpositions) {			
			print $out $sortedpositions{$position};
		}
	}
}

sub gff2coo {
	my $usage = qq/
gff2coo(<gff_file>, <coo_outfile)
sub gff2coo (partly courtesy of J. Slot) formats GFF files from JGI 
into coo3 files used by various clustering pipelines\n\n/;
	die($usage, "Error: incorrect number of args in gff2coo\n") if (scalar @_ != 2);
	my ($gff_file, $coo3_outfile) = @_;
	my (%coo_hash);
	open(my $gff_in, '<', $gff_file) or die($usage, "Error: could not open $gff_file for reading\n");
	my ($gfffilename) = fileparse($gff_file);
	$gfffilename =~ m/^([^\.]+)/;
	my $ome = $1;
	while (my $line = <$gff_in>) {
		next if ($line =~ m/^#/ || $line =~ m/^\s/);
		chomp $line;
		$line = $line . ";";

		#acquire info
		my @fields = split ("\t", $line);
		my ($contig, $boundary1, $boundary2, $orientation) = ($fields[0], $fields[3], $fields[4], $fields[6]);
		#the following is for gff3 files from JGI
		if ($fields[2] eq "CDS"){				
			my $accession;
			$fields[8] .= " ";
			if ($fields[8] =~ m/[Pp]roteinId (.+?)[ ;]/) {$accession = $1;$accession =~s/"//g;}#jgi
			elsif ($fields[8] =~ m/transcript_id (.+?)[ ;]/s) {$accession = $1;$accession =~s/"//g;}#broad
			elsif ($fields[8] =~ m/protein_id=(.+?)[ ;]/s) {$accession = $1;$accession =~s/"//g;}#ncbi
			else {print "no match $line\n"}
			#hash data to make .coo file
			#print "$fields[0]\n";
			my $add_protein = "$ome"."_"."$accession";
			$coo_hash{$ome}{$contig}{$add_protein}{$boundary1} = $orientation;
			$coo_hash{$ome}{$contig}{$add_protein}{$boundary2} = $orientation;
		}
	}
	open (my $coo_out, '>', $coo3_outfile) or die($usage, "Error: cannot open $coo3_outfile for writing\n");
	$coo_out->autoflush(1);
	foreach my $ome (sort keys %coo_hash) {
		my (@coo_array);
		foreach my $locus (sort {$a cmp $b} keys %{$coo_hash{$ome}}) {
			foreach my $protein (keys %{$coo_hash{$ome}{$locus}}) {
				my @sort_coords = (sort {$a <=> $b} keys %{$coo_hash{$ome}{$locus}{$protein}});
				my ($min, $max) = @sort_coords;
				my $strand = $coo_hash{$ome}{$locus}{$protein}{$min};
				my @outputline = ($protein, $locus, $min, $max, $strand);
				push (@coo_array, [@outputline]);	
			}
		}
		#process .coo info
		#$VAR10637 = 'Ecol19_NP_052627.1 Ecol19_NC_002128.1 23512 23793';
		my $i = 0;
		foreach my $protein_location ((sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @coo_array)) {
			$i++;
			#print .coo file
			#$protein, $locus, $downstream_position, $upstream_position, $strand
			print $coo_out "@$protein_location\n";
		}
	}
}


