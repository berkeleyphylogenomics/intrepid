#! /usr/bin/perl -w
use lib '/usr/lib/perl5/site_perl';
use lib '/usr/lib/perl5/5.8.3';
use lib '/usr/lib/perl5/site_perl/5.8.3';

#use lib '/sw/lib/perl5/5.8.6';
#use lib '/clusterfs/ohana/software/lib/perl/lib/perl5/site_perl/5.8.8/';
use lib '/home/sriram_s/perl/lib/perl5/site_perl';
use lib '/home/sriram_s/perl/lib/perl5/5.8.7';
use lib '/home/sriram_s/perl/lib/perl5/site_perl/5.8.7';

use strict;
use Getopt::Long;
$|=1;
use Bio::AlignIO;
use Bio::Structure::IO;
use Bio::Seq;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::TreeIO;
use Bio::Tree::NodeI;

use Class::Struct;

#
#
# Author 	- Sriram Sankararaman (sriram_s@cs.berkeley.edu)
#

#
# Canonical gap character
my $gapCharacter 	= "-";

# Character to denote subfamily conservation
my $subfamCharacter	= "Z";



struct Dmm => {
	component	=> '$',
	weights		=> '@',
	alpha		=> '$',
	aamap		=> '%',
	totalalpha 	=> '@'
};
my %parameters=();
initParameters (\%parameters);
initMaps (\%parameters);
print "verbose1=$parameters{verbose1}\n";
print "verbose2=$parameters{verbose2}\n";
print "verbose3=$parameters{verbose3}\n";
run(\%parameters);
done (\%parameters);



# The main method
#
sub run{
	my %parameters	= %{$_[0]};
	my $pic 		= $parameters{pic};
	my $verbose		= $parameters{verbose2};
	my $treeFile	= $parameters{tree_file};
	my $msaFile		= $parameters{msa_file};
	my $outdir		= $parameters{outdir};
	my $loghandle	= $parameters{loghandle};
	my @mapnames = @{$parameters{map_names}};

	my $map1 = "importance_map";

	print "Running INTREPID\n";
	globalConservation (\%parameters);
	treeTrace(\%parameters);
	for (my $i = 0; $i <@mapnames; $i++){
		print "Run $i:$mapnames[$i]\n";
		$parameters{importance_map} = $parameters{$mapnames[$i]};
		print $loghandle "Completed run $i\n";
	}
	


	my @keys = sort keys %{$parameters{indicator_maps}};
	printResult (\%parameters);
}



sub printResult{
	my %parameters		= %{$_[0]};
	my $loghandle		= $parameters{loghandle};
	my %indicatorMaps 	= %{$parameters{indicator_maps}};
	my %seqToStruct		= %{$parameters{seq_to_struct}};
	my $verbose			= $parameters{verbose2};	
	my %auxMap			= %{$parameters{aux_map}};
	my %posMap 			= %{$parameters{pos_map}};
	my $pdbId			= $parameters{structure_id};
	my $chain			= $parameters{chain};
	my $outdir			= $parameters{outdir};
	my $auxFile			= "output.aux";
	my $rankFile		= "output.rank";
	my %originalMap		= %{$parameters{original_map}};
	my @mapnames 		= @{$parameters{map_names}};
	$auxFile			= $outdir."/".$auxFile if defined $outdir;
	$rankFile			= $outdir."/".$rankFile if defined $outdir;

	open ( AUXFILE,  ">$auxFile");
	my @keys	= @mapnames;
	print AUXFILE "#Position\tResSeq\tChain\tResidue\tICode\t".join("\t",@keys)."\n";

	foreach my $originalSeqPos (sort {$a<=>$b} keys %seqToStruct){
		my @colorArray;
		foreach my $key ( @keys ){
			my %map = %{$parameters{$key}};
			my $seqPos = $originalMap{$originalSeqPos};
			my $color = 0;
			print $loghandle "Map=$key\toriginalSeqPos= $originalSeqPos\tseqPos=$seqPos\tDefined?=".defined $map{$seqPos};
			print $loghandle "\n";
			$color = $map{$seqPos} if (defined $map{$seqPos});
			if ($parameters{verbose1}){
				print $loghandle "Seq=${originalSeqPos}, Key=${key}, Color=$color\n";
			}
			push (@colorArray, $color);
		}		
		my $pdbPos = $seqToStruct{$originalSeqPos};
		print AUXFILE $originalSeqPos."|".$pdbPos."|".$auxMap{$pdbPos}."|".arrayToString("%5.3f",\@colorArray,"|")."\n";
	}
	close AUXFILE;


	open ( RANKFILE,  ">$rankFile");
	@keys	= @mapnames;
	print RANKFILE "#Position\tResSeq\tChain\tResidue\tICode\t".join("\t",@keys)."\n";

	foreach my $originalSeqPos (sort {$a<=>$b} keys %seqToStruct){
		my @colorArray;
		foreach my $key ( @keys ){
			my %map = %{$parameters{"rank_".$key}};
			my $seqPos = $originalMap{$originalSeqPos};
			my $color = 0;
			print $loghandle "Map=$key\toriginalSeqPos= $originalSeqPos\tseqPos=$seqPos\tDefined?=".defined $map{$seqPos};
			print $loghandle "\n";
			$color = $map{$seqPos} if (defined $map{$seqPos});
			if ($parameters{verbose1}){
				print $loghandle "Seq=${originalSeqPos}, Key=${key}, Color=$color\n";
			}
			push (@colorArray, $color);
		}		
		my $pdbPos = $seqToStruct{$originalSeqPos};
		print RANKFILE $originalSeqPos."|".$pdbPos."|".$auxMap{$pdbPos}."|".arrayToString("%4d",\@colorArray,"|")."\n";
	}
	close RANKFILE;

}




sub initScoreMap{
	my %parameters	= %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	my $scoreFile	= $parameters{score_matrix_file};
	my $verbose		= $parameters{verbose2};
	my @row=();
	my %scoreMap=();
	open(FILE,$scoreFile);
	while (<FILE>){
		chomp;
        # $. is the line number
		s/^\s+// if $.==1;
		@row=split (/\s+/) if $.==1;
		if ($.>1){
			@_=split (/\s+/);
			my $label=shift;
			my $min = 100;
			my $max = -100;
			foreach my $r (@row){
				$label=~s/\*|\./$gapCharacter/g;
				$r=~s/\*|\./$gapCharacter/g;
				my $value = shift;
				$scoreMap{uc($label).uc($r)}=$value;
				$scoreMap{uc($r).uc($label)}=$value;
				$max = $value if $value > $max;
				$min = $value if $value < $min;
			}
			$scoreMap{$label}="$min:$max";
		}
	}
	close FILE;
	if ($verbose){
		foreach my $key (keys %scoreMap){
			print $loghandle $key."\t".$scoreMap{$key}."\n";
		}
	}
	return %scoreMap;
}

sub slice{
	my $aln				= shift;
	my ($start, $end)	= @_;
	my $newAln				= new Bio::SimpleAlign();
	foreach my $seq ($aln->each_seq()){
		my $slice=$seq->subseq($start, $end);
		my $newSeq=Bio::LocatableSeq->new(	-seq=>$slice,
											-id=>$seq->id,
											-start=>$seq->start,
											-end=>$seq->end,
											-alphabet=>"protein"
											);
		$newAln->add_seq($newSeq);
	}
	return $newAln;
}


sub dist{
	my ($atom1,$atom2)=@_;
	my ($x1,$y1,$z1)=($atom1->x(),$atom1->y(), $atom1->z());
	my ($x2,$y2,$z2)=($atom2->x(),$atom2->y(), $atom2->z());
	return sqrt(($x1-$x2)**2+($y1-$y2)**2+($z1-$z2)**2);
}


sub cluster1d{
	my %parameters	= %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	my $pdbFile		= $parameters{pdb_file};
	my $sequenceID	= $parameters{sequence_id};
	my $msaFile		= $parameters{msa_file};
	my %posMap		= %{$parameters{pos_map}};
	my %aaMap		= %{$parameters{amino_acid_map}};
	my %structToSeq = %{$parameters{struct_to_seq}};
	my %seqToStruct	= %{$parameters{seq_to_struct}};
	my $verbose		= $parameters{verbose2};
	my %columnScoreMap 	= %{$parameters{importance_map}};

	
	my $residues = keys %columnScoreMap;
	my $window = 5;

	my @clusterVector = ();
	my %clusterMap = ();
	for (my $i = 1; $i < $residues; $i++){
		my $count = 0;
		my $avg = 0;
		for (my $j = $i-$window; $j <= $i+$window; $j++){
			next if ($j  < 1 || $j >= $residues);
			next if (!defined $columnScoreMap{$j});
			$avg += $columnScoreMap {$j};
			$count ++;
		}
		$avg = $avg/$count if ($count > 0);
		
		$clusterVector[$i] = $avg;
	}


	my @clusterz = normalize (\@clusterVector);
	for (my $i = 0; $i < @clusterz ; $i ++){
		next if (!defined $clusterz[$i]);
		$clusterMap{$i} = $clusterz[$i];
	}
	
	if ($parameters{verbose1}){
		print $loghandle "SP map\n";
		foreach my $key (keys %clusterMap){
			next if !defined $clusterMap{$key};
			print $loghandle $key."\t".$clusterMap{$key}."\n";
		}
	}
	$parameters{cluster_map}	= \%clusterMap;
	%{$_[0]}	= %parameters;
}


sub cluster{
	my %parameters	= %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	my $pdbFile		= $parameters{pdb_file};
	my $sequenceID	= $parameters{sequence_id};
	my $msaFile		= $parameters{msa_file};
	my %posMap		= %{$parameters{pos_map}};
	my %aaMap		= %{$parameters{amino_acid_map}};
	my %structToSeq = %{$parameters{struct_to_seq}};
	my %seqToStruct	= %{$parameters{seq_to_struct}};
	my $verbose		= $parameters{verbose2};
	my %columnScoreMap 	= %{$parameters{importance_map}};
	

	if (defined $parameters{importance_map}){
			print "importance_map defined\n";
		} else {
			print "importance_map undefined\n";
		}
	my $weightFunction	= $parameters{weight_function};
	my $clusterDistance;
	

	if ($weightFunction eq "step" ){
		$clusterDistance = $parameters{cluster_distance};
	} elsif ($weightFunction eq "exp" ){
		$clusterDistance = $parameters{cluster_distance};
	} else {
		die "Error: Unknown weight function\n";
	}
	
	
	
	# The sequence from the MSA
	my $sequencePrimary="";

	# The sequence from the PDB file
	my $structurePrimary="";

	my @coords=();
	# Positions in the 3d structure
	my @positions=();
	my @clusterVector = ();
	my %clusterMap = ();

	print "Computing structural proximity\n";
	print $loghandle "Computing structural proximity\n";

	if  ( $parameters{only_sequence}){
		foreach my $key (keys %columnScoreMap){
			$clusterMap{$key}  = $columnScoreMap{$key};		
		}
		$parameters{cluster_map}	= \%clusterMap;
		%{$_[0]}	= %parameters;
		return;
	}

	my $struct 		= $parameters{structure};
	my $index		= 0;
	
	for my $chain ($struct -> get_chains ){
		my $chainid = $chain->id;

		if($chainid eq $parameters{chain}){
			print $loghandle "Found chain:$chainid\n" if $verbose;
			for my $res ($struct->get_residues($chain)){
				my $resid  = $res->id;
				my ($residue,$pos)=split (/-/,$resid);
				# Checking for strange residues or HETATM records.
				next if !defined $aaMap{$residue};
				my @atoms  = $struct->get_atoms($res);
				$index++;
				foreach my $atom (@atoms){
					if ($atom->pdb_atomname =~ /CA/){
						$coords[$index]		= $atom;
						$positions[$index]	= $pos;
					}
				}
			}
		}
	}

	print $loghandle "SP score\n";
	my $counter = 0;
	for (my $i=1; $i<=$#coords; $i++){
		print $loghandle "Finding neighbors for residue :$i\n" if $verbose;
		my $countCluster=0;
		my $avgClusterScore=0;
		next if (!defined $structToSeq{$positions[$i]});
		my $sequencePosition=$structToSeq{$positions[$i]};

		$counter++;
		
		my @neighbors=();
		my @scores = ();
		for (my $j=1 ;$j<=$#coords ; $j++){
			print $loghandle "Examining residue:$j\n" if $verbose;
			if (!defined $structToSeq{$positions[$j]}) {
				my $distance = dist ($coords[$i], $coords[$j]);
				print $loghandle "Undefined j=$j,pos=".$positions[$j].",dist=".$distance."\n";
			}
			next if (!defined $coords[$i] || !defined $coords[$j]);
			next if (!defined $structToSeq{$positions[$j]});
			my $distance = dist ($coords[$i], $coords[$j]);
			my $weight = 0;
			if ($weightFunction eq "step" ){
				$weight = 1 if ($distance < $clusterDistance);
			} elsif ($weightFunction eq "exp" ){
				$weight = 2**(-$distance/$clusterDistance);
			}
			
			my $sequencePosition=$structToSeq{$positions[$j]};
			if ($sequencePosition > 0){
				$countCluster 	+= $weight;
				$avgClusterScore+= $weight*$columnScoreMap{$sequencePosition} if defined $columnScoreMap{$sequencePosition};
				push (@scores, $columnScoreMap{$sequencePosition}) if (defined $columnScoreMap{$sequencePosition} && $weight>0);
			}
			
			if ($distance < $clusterDistance){
#				if ($verbose){
					print $loghandle "j=$j, sequencePosition=".$sequencePosition."\n";
#				}
				push (@neighbors, $j);
			}
		}
		$avgClusterScore/=$countCluster if ($countCluster>0);
		print $loghandle "$i\tStructure position=$positions[$i]\tNeighbors= $countCluster\tScore=$avgClusterScore\n";
		print $loghandle "Scores=".join(",",@scores)."\n";
		print $loghandle join(",",@neighbors)."\n";
		$clusterVector[$counter] = $avgClusterScore;
	}

	my @clusterz = normalize (\@clusterVector);
#	my @clusterz = @clusterVector;
	for (my $i = 0; $i < @clusterz ; $i ++){
		next if (!defined $clusterz[$i]);
		$clusterMap{$i} = $clusterz[$i];
	}
	
	if ($parameters{verbose1}){
		print $loghandle "SP map\n";
		foreach my $key (sort {$a<=>$b} keys %clusterMap){
			next if !defined $clusterMap{$key};
			print $loghandle $key."\t".$clusterMap{$key}."\n";
		}
	}
	$parameters{cluster_map}	= \%clusterMap;
	%{$_[0]}	= %parameters;
}



#
# Compute the solvent exposure of amino acids
# 
sub solventExposure {
	my %parameters	= %{$_[0]};
	return if (!defined $parameters{dssp});
	my $buriedThreshold = -1;
	$buriedThreshold = $parameters{buried_threshold};
	my %solventMap = ();

	if ($buriedThreshold >= 0){
		%solventMap = readMap (\%parameters);
	}

	$parameters {solvent_map} = \%solventMap;
}


sub readDSSP{
	my %parameters	= %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	my $file 		= $parameters{dssp};
	my %areaMap		= %{$parameters {area_map}};
	my $chainId 	= $parameters {chain};
	my $indir  = $parameters{indir};
	my $prefix = $indir."/" if $indir;

	$file = "$prefix/$file";
	open (FILE, $file);
	my $flag =0 ;
	my %map = ();
	while (<FILE>){
		chomp;
		my @tmp  = split (// );
		next if @tmp<2;
		$flag = 1 if ($tmp[2] eq "#");
		if ($flag){
			my $line = $_;
			my $pos = join("",@tmp[6..9]);
			$pos=~s/^[\s]+//g;
			next if !($pos=~/^[0-9]+$/);
			my $chain = $tmp[11];
			next if !($chainId eq $chain);
			my $res = $tmp[13];
			next if ($res eq "X");
			next if ($res eq "");
			my $acc = join("",@tmp [35..37]);
			$acc =~s/^[\s]+//g;
			next if !($acc=~/^[0-9]+$/);
			$map{$pos} = 0;

			if ($parameters {verbose2}) {
				print "pos=$pos,res=$res:acc =$acc:".$areaMap{$res}."\n";
			}
			$map{$pos} = $acc/$areaMap{$res} if( defined $areaMap{$res} && $areaMap{$res}>0);
		}
	}

	if ($parameters {verbose2}){
		foreach my $key ( sort {$a<=>$b}keys %map){
			print $key."\t".$map{$key}."\n";
		}
	}
	return %map;
}


sub columnToPosition{
	my $sequenceID	= $_[0];
	my $aln		 	= $_[1];
	my $column		= $_[2];

	my $length	= $aln->length();
	return -1 if ($column > $length);
	
	my $sequencePrimary;
	foreach my $seq ($aln->each_seq_with_id($sequenceID)){
		 $sequencePrimary = $seq;
	}

	my $position=0;
	for(my $i=1; $i<=$length && $i<=$column; $i++){
		if ($sequencePrimary->subseq($i,$i) =~ /[A-Z]/ ){
			$position++;
		}
	}
	return $position;
}



sub dmmrank {
	my $dmm 		= $parameters{dmm_class};
	my $loghandle	= $parameters{loghandle};
	my $verbose		= $parameters{verbose1};
	my %validMap	= %{$parameters{valid_amino_acid}};

	# Using DMM
	my @dmre = ();
	my @blos = ();
	my @bg = @{$_[0]};
	my @profile = @{$_[1]};

	if ( $verbose ){ 
		print $loghandle "In dmmrank\n";
	}

	for (my $i = 0; $i < @profile; $i++){
		my %map = %{$profile[$i]};
		my %bgMap = %{$bg[$i]};

# Keep track of gap characters
		my $pseudoCount ;
		my $bgpseudoCount ;
		$pseudoCount = $map{$gapCharacter} ;
		$bgpseudoCount = $bgMap{$gapCharacter};
# We want gaps that are not present in the current division
		$bgpseudoCount -= $pseudoCount;

		my %p = ();
		my %q = ();
		my $denomp = 0;
		my $denomq = 0;
		foreach my $key (keys %validMap){
			my $n = $map{$key} ;
			my $m = $bgMap {$key} ;

# We want characters that are not present in the current division.
			$m -= $n;
			$p{$key}  = $n;
			$q{$key}  = $m;
			if ( $parameters{verbose2}) {
				print $loghandle $key."=>".$p{$key}."\t".$q{$key}."\n";
			}
			$denomp += $map{$key};
			$denomq += $bgMap{$key};
		}

		$denomq -= $denomp;

# Use pseudocount=1/20 for gap characters
		foreach my $key (keys %validMap){
			$p{$key} += $pseudoCount * 0.05;

			$q{$key} += $bgpseudoCount * 0.05;
		}


# Use DMM to estimate posterior means
		my @w = @{$dmm->weights};
		my $l = $dmm->component;
		my @alpha = @{$dmm->alpha};
#		my %serialmap = %{$dmm->aamap};

		my @p = ();
		my @q = ();

		if ($parameters{verbose2}) { 
			print $loghandle "w=@w\n";
			print $loghandle "l=$l\n";
		}

		for (my $i = 0; $i < 20; $i++){
			my $n1 = $p{$dmm->aamap($i)};
			my $n2 = $q{$dmm->aamap($i)};

			if ( $parameters{verbose2}) {
				print $loghandle "aamap($i)=".$dmm->aamap($i)."\n";
				print $loghandle "n1=$n1,n2=$n2\n";
			}
			
			my $c1 = 0;
			my $c2 = 0;
			for (my $j = 0 ; $j < $l; $j++){
				my $alpha = $alpha[$j][$i];
				my $total = $dmm->totalalpha($j);
				if ( $parameters{verbose2} ){ 
					print $loghandle "$i,$j\t$alpha,$total\n";
				}

				$c1  += $w[$j] * ( $alpha  + $n1)/($total + $denomp);
				$c2  += $w[$j] * ( $alpha  + $n2)/($total + $denomq);
			}

			push (@p, $c1);
			push (@q, $c2);
		}


		if ($i==1 && $parameters{verbose2})  {
			for (my $i = 0; $i < 20; $i++){
				print $loghandle $dmm->aamap($i)."\t";
			}
			print $loghandle "\n";
			print $loghandle arrayToString("%5.3f", \@p,",")."\n";
			print $loghandle arrayToString("%5.3f", \@q,",")."\n";
		}

# Ignore very gappy columns
		my $gaps = 0;
		my $gaps1 = 0;
		my $gaps2 = 0;
		$gaps1 = $pseudoCount/$denomp if ($denomp > 0);
		$gaps2 = $bgpseudoCount/($denomq+$denomp) if ($denomq > 0);
		$gaps = ($pseudoCount+ $bgpseudoCount)/($denomp+$denomq) if ($denomp+$denomq >0);
		my $re;
#		if ($gaps1 > 0.4 || $gaps2 > 0.4) {
#			$re = 0;
#		} else {
#			$re = relativeEntropy (\@p, \@q) ;
#		}
		$re = relativeEntropy (\@p,\@q) * (1- $gaps);
		$dmre[$i] = $re;
		$blos[$i] = getblosum (\@p,\@q);

	}
	@{$_[2]} = @dmre;
	@{$_[3]} = @blos;


	if ($parameters{verbose1}){
		for (my $i = 0 ; $i < @dmre; $i++){
			print $loghandle "$i:$dmre[$i]\t";
		}
		print $loghandle "\n";
	}
	if ( $verbose ){ 
		print $loghandle "Out of  dmmrank\n";
	}
}


sub getblosum {
	my $dmm 		= $parameters{dmm_class};
	my $verbose		= $parameters{verbose1};
	my $loghandle	= $parameters{loghandle};
	my %scoreMap	= %{$parameters{score_matrix}};
	my @p = @{$_[0]};
	my @q = @{$_[1]};

	my $pindex = -1;
	my $pmax = 0;
	my $qindex = -1;
	my $qmax = 0;

	for (my $i = 0 ;  $i < @p; $i++){
		if ($p[$i]  > $pmax ){
			$pmax = $p[$i];
			$pindex = $i;
		}
	}


	for (my $i = 0 ;  $i < @q; $i++){
		if ($q[$i]  > $qmax ){
			$qmax = $q[$i];
			$qindex = $i;
		}
	}

	my $a1  = $dmm->aamap ($pindex);
	my $a2  = $dmm->aamap ($qindex);

	my $score = 0;
	if ($a1 ne $a2 ){
		$score=  $scoreMap {$a1.$a2};
	}

	if ($verbose) { 
		print  $loghandle "$a1,$a2,$score\n";
	}

	return $score;
	
}

sub logoddsrank {
	my $loghandle	= $parameters{loghandle};
	my $verbose		= $parameters{verbose2};
	my %validMap	= %{$parameters{valid_amino_acid}};
	my %background  = %{$parameters{background_distribution}};

	my @bg = @{$_[0]};
	my @profile = @{$_[1]};
	my @q = ();
	foreach my $key (keys %validMap ) { 
		push  (@q, $background{$key});
	}


	#
	# Computing the log-odds rank
	#

	my $averagere = 0;
	my $averagejs = 0;
	my $averageCons = 0;
	my $total = 0;
	my @conservation = ();
	my @reConservation = ();
	my @jsConservation = ();
	for (my $i=0; $i<@profile; $i++){
		my $denom = 0;
		my $max = 0;
		my $maxKey;
		my %map = %{$profile[$i]};
		foreach my $key (keys %validMap){
			if ($map{$key} > $max){
				$max = $map{$key};
				$maxKey = $key;
			}
			$denom += $map{$key};
		}
		$total++ if ( defined $maxKey && defined $validMap{$maxKey});
		$conservation[$i] 	= 0;
		$conservation[$i] 	= $max/$denom if ($denom > 0 && $maxKey ne $gapCharacter);
		$conservation[$i] *= (1 - $map{$gapCharacter}/($denom+$map{$gapCharacter}));
		$averageCons 	+= $conservation[$i];


		my @p = ();
		my @r = ();
		foreach my $key (keys %validMap ) { 
			if ($denom > 0 ) {
				push (@p, $map{$key}/$denom);
				push (@r, 0.5 * $map{$key}/$denom + 0.5* $background {$key});
			} else {
				push (@p, 0);
				push (@r, 0);
			}
		}
		$reConservation[$i] = relativeEntropy (\@p, \@q);
		$jsConservation[$i] = 0.5 * relativeEntropy(\@p,\@r) + 0.5 * relativeEntropy (\@q,\@r);
		$reConservation[$i] *= (1 - $map{$gapCharacter}/($denom+$map{$gapCharacter}));
		$jsConservation[$i] *= (1 - $map{$gapCharacter}/($denom+$map{$gapCharacter}));
		$averagere += $reConservation[$i];
		$averagejs += $jsConservation[$i];
	}
	$averageCons /= $total if ($total > 0);
	$averagere /= $total if ($total > 0);
	$averagejs /= $total if ($total > 0);


	if ($parameters{verbose1}){
		print $loghandle "Log-odds score\n";
		print $loghandle "averageCons=$averageCons\n";
		print $loghandle "averagejs=$averagejs\n";
		print $loghandle arrayToString("%5.3f", \@conservation,",")."\n";
		print $loghandle arrayToString("%5.3f", \@jsConservation,",")."\n";
	}
	@{$_[2]} = @conservation;
	@{$_[3]} = @reConservation;
	@{$_[4]} = @jsConservation;
	return ($averageCons, $averagere, $averagejs);
}


sub rerank {
	my $loghandle	= $parameters{loghandle};
	my $verbose		= $parameters{verbose2};
	my %validMap	= %{$parameters{valid_amino_acid}};
	my $bitThreshold 	= $parameters{bit_threshold};

	my @bg = @{$_[0]};
	my @profile = @{$_[1]};

	if ($parameters{verbose1}) { 
		print $loghandle "In rerank\n";
	}

#
# Computing the entropy rank
#
	my $averageBits = 0;
	my $averageBits2 = 0;
	my @bits = ();
	my @re = ();
	my $total = 0;
	for (my $i=0; $i<@profile; $i++){
		$bits [$i] = log(20)/log(2);
		my %map = %{$profile[$i]};
		my %bgMap = %{$bg[$i]};
		my $pseudoCount = 0;
		my $bgpseudoCount = 0;
		$pseudoCount = $map{$gapCharacter}/20 ;
		$bgpseudoCount = $bgMap{$gapCharacter}/20;
# We want gaps that are not present in the current division
		$bgpseudoCount -= $pseudoCount;

		my @p = ();
		my @q = ();
		my $denom = 0;
		foreach my $key (keys %validMap){
			my $n = $map{$key} ;
			my $m = $bgMap {$key} ;
# New addition
			$m -= $n;
			push (@p, $n + $pseudoCount);
			push (@q, $m + $bgpseudoCount);
			$denom += $map{$key};
		}
		if ($parameters{verbose1}) { 
			print $loghandle "p,q\tPosition:$i\n";
			foreach my $key (keys %validMap){
				print $loghandle $key."\t";
			}
			print $loghandle "\n";
			print $loghandle join (",",@p)."\n";
			print $loghandle join (",",@q)."\n";
			if ($i==29) {
				print "RE:" .relativeEntropy (\@p,\@q)."\n";
			}

		}

# Ignore very gappy columns
		if ($pseudoCount > 0 && 20*$denom/$pseudoCount < $bitThreshold){
			print "Too many gaps in position $i\n";
			next;
		}
		my $entropy = entropy (\@p);
		print $loghandle "Entropy = $entropy\n" if ($parameters{verbose2});
		$bits [$i] 	= $entropy;

		my $re = relativeEntropy (\@p, \@q)  ;
		$re[$i] = $re;

		$averageBits += $bits[$i];
		$averageBits2 += $bits[$i]*$bits[$i];
		$total ++;
	}
	$averageBits /= $total if ($total > 0);
	$averageBits2 /= $total if ($total > 0);
	$averageBits2 -= $averageBits*$averageBits;
	$averageBits2 = sqrt ($averageBits2);
	$averageBits2 = 1000 if ($averageBits2 == 0);

		
	if ($parameters{verbose1}){
		print $loghandle "Entropy score\n";
		print $loghandle "averageBits=$averageBits\n";
		print $loghandle "Bits:".@bits."\n";
		print $loghandle arrayToString("%5.3f", \@bits,",")."\n";
		print $loghandle "RE:\n";
		print $loghandle arrayToString("%5.3f", \@re,",")."\n";

	}
	@{$_[2]} = @bits;
	@{$_[3]} = @re;

	return $averageBits;
}


sub treeTrace{
	my %seqToStruct		= %{$parameters{seq_to_struct}};
	my %parameters	= %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	my $sequenceID	= $parameters{sequence_id};
	my $verbose		= $parameters{verbose2};
	my %validMap	= %{$parameters{valid_amino_acid}};
	my $dmm 		= $parameters{dmm_class};
	if ($verbose){ 
		print "In treeTrace\n";
	}
		print "In treeTrace\n";
	
	my %rankImportanceMap = ();
	my %rankREMap = ();
	my %rankJSMap = ();

	my %consImportanceMap = ();
	my %consREMap = ();
	my %consJSMap = ();
	my $entropy		= $parameters{entropy};

	my @oddsImportance = ();
	my @reImportance = ();
	my @jsImportance = ();

	my @dmreMatrix 		= ();
	my @blosMatrix 		= ();
	my @blos 			= ();
	my @dmre			= ();
	my %rankdmreMap			= ();
	my %dmreMap			= ();
	
	
	my $aln = $parameters{msa};
	my $tree= $parameters{tree};
	
	my $root=$tree->get_root_node;
	my $designatedNode;
	my $subtreeScore = $parameters{subtree_score};


	# The array that contains all the sequences 
	my @seqMatrix = ();
	# Maps the sequence ids to their indices in the seqMatrix
	my %id2indexMap = ();
	# The rank of each sequence w.r.t the tree trace
	my @subtreeRank = ();

	my @desc =  $root->get_Descendents ();
	print $loghandle "descendents = ".@desc."\n";

	# Find out which node in the tree contains the reference sequence 
	foreach my $node ($root->get_Descendents()){
		if ($node->is_Leaf && (index ($node->id, $sequenceID)>=0)){
			$designatedNode = $node;
			last;
		}
	}

	my @trace = ();
	if (defined $designatedNode){
		push (@trace, $designatedNode);
		my $tmp = $designatedNode;
		while (defined $tmp->ancestor){
			$tmp = $tmp -> ancestor;
			unshift ( @trace, $tmp);
		}
	} else {
		die "Error:".$parameters{tree_file}. " does not contain $sequenceID\n";
	}
	
	# Find out which sequence object is the reference sequence
	my $seq ;
	foreach my $seqObj ($aln->each_seq){
		if ( index($seqObj->id, $sequenceID)>=0 ){
			$seq = $seqObj;
			last;
		}
	}

	# Initialize the sequence matrix
	# id2indexMap maps each node/sequence id to the sequence.
	foreach my $seq  ($aln->each_seq){
		push (@seqMatrix, $seq->seq);
		$id2indexMap{$seq->id} = $#seqMatrix;
	}
	
	$parameters {seq_matrix} = \@seqMatrix;
	$parameters {id2index_map} = \%id2indexMap;
#	$parameters {subtree_rank} = \@subtreeRank;
	@subtreeRank = @{$parameters {subtree_rank}};
	
	my $level=0;
	my $background ;

	
	# For each node in the tree, 
	# find all the sequences in the subtree.
	# Assemble into an array.
	# Compute conservation at each position and the average across positions.
	# Compute importance matrix.
	foreach my $node (@trace){
		my @descendents;
		$level++;
		if ($node->is_Leaf){
			unshift (@descendents, $node);
		}else{
			@descendents= $node->get_Descendents();
		}
		print $loghandle "level=$level\tdescendents=".@descendents."\n";
		
		foreach my $descendent (@descendents){
			if ($descendent->is_Leaf ){
				print $loghandle $descendent->id."\n";
			}
		}
		
		my @seqArray=();
		foreach my $seqObj ($aln->each_seq){
			my $flag =  0;
			foreach my $descendent (@descendents){
				if ($descendent->is_Leaf && $seqObj->id eq $descendent->id){
					unshift (@seqArray, $seqObj);

					# Get the subtree rank for each sequence
					$subtreeRank [ $id2indexMap{$descendent->id}] = $level;
					$flag = 1;
				}

			}
			print $loghandle "Added ".$seqObj->id."?".$flag."\n";
		}

		print $loghandle "level=$level\t isroot? " .(! defined $node->ancestor). "\tleaves=" .@seqArray."\n";
		for (my $i1 = 0; $i1 < @seqArray; $i1++){
			my $obj  = $seqArray[$i1];
			print $loghandle $obj->id."\n";
		}
		
		if ($verbose ){
			print "SubtreeRank map\n";
			foreach my $key (sort keys %id2indexMap){
				print $key."\t".$subtreeRank[$id2indexMap{$key}]."\n";
			}
		}

		$parameters {subtree_rank} = \@subtreeRank;
		my @profile = ();
		getProfile (\%parameters, $level, \@profile);
		if ($parameters{verbose1}){
			print $loghandle "Profile length=".@profile."\n";
		}
		print $loghandle "level=$level \t Profile length=".@profile."\tSequences=".@seqArray."\n";
		
		$background = \@profile if ($level == 1);
		
		
		# Log-odds
		my @conservation = ();
		my @jsConservation = ();
		my @reConservation = ();
		my ($averageCons, $averagere, $averagejs) = logoddsrank ($background, \@profile, \@conservation, \@reConservation, \@jsConservation);

		for (my $i=0; $i<@conservation;$i++){
			$oddsImportance[$level][$i+1]=0;
			$reImportance[$level][$i+1]=0;
			$jsImportance[$level][$i+1]=0;
			$oddsImportance[$level][$i+1]=log ($conservation[$i]/$averageCons) if ($averageCons>0 && $conservation[$i]>0);
			$reImportance[$level][$i+1]=($reConservation[$i] - $averagere) if ($averagere>0 && $reConservation[$i]>0);
			$jsImportance[$level][$i+1]=($jsConservation[$i] - $averagejs) if ($averagejs>0 && $jsConservation[$i]>0);
			if ($parameters{verbose1}){
				print $loghandle "Level=$level, i=".($i+1).":".$jsConservation[$i]."\t".$averagejs."\t".$jsImportance[$level][$i+1]."\n";
			}
		}


		# DMM
		dmmrank ($background, \@profile, \@dmre, \@blos);

		
		for (my $i=0; $i<@dmre;$i++){
			$dmreMatrix [$level][$i+1] = $dmre[$i];
		}
	}


	print $loghandle "Levels=$level\n";
	print $loghandle "Importance matrix\n";

	#
	# Determine the final importance scores for each position 
	# based on the profile matrices.
	#
	my @consImportanceVector = ();
	my @consREVector = ();
	my @consJSVector = ();
	my @reVector = ();
	my @dmreVector =  ();
	my @dmrePointer = ();
	my @repointer = ();
	my @maxLevel = ();
	my @jsmaxLevel = ();
	my $length = $aln->length;

	if ($parameters{verbose1}){
		print $loghandle "Matrices\n";
	}


	for (my $i=1, my $k=1; $i<= $length ; $i++){
		if ($parameters{verbose1}){
			print $loghandle "\n$i:\t";
		}
		# Skip gaps in the reference sequence
		next if (!defined ($validMap{$seq->subseq($i,$i)}));

		$consImportanceVector[$k]	= 0;
		$consREVector[$k]	= 0;
		$consJSVector[$k]	= 0;
		if ($verbose){
			print $loghandle $k."\t".$seq->subseq($i,$i)."\t";
		}

		
		my $tmp = 0;
		foreach (my $j=1; $j<= $level; $j++){
			if ($j == 1){
				$consImportanceVector[$k]=$oddsImportance[$j][$i];
				$consREVector[$k]=$reImportance[$j][$i];
				$consJSVector[$k]=$jsImportance[$j][$i];
				$maxLevel [$k] = $j;
				$jsmaxLevel [$k] = $j;

				$dmreVector[$k] = 0;
				$dmrePointer[$k] = 1;
			}

			
			if ($oddsImportance[$j][$i]> $consImportanceVector[$k]){
				$consImportanceVector[$k]=$oddsImportance[$j][$i];
				$maxLevel [$k] = $j;
			}
			if ($reImportance[$j][$i]> $consREVector[$k]){
				$consREVector[$k]=$reImportance[$j][$i];
			}
			if ($jsImportance[$j][$i]> $consJSVector[$k]){
				$consJSVector[$k]=$jsImportance[$j][$i];
				$jsmaxLevel [$k] = $j;
			}
	

			if ($j > 1  && $dmreVector[$k] < $dmreMatrix[$j][$i] ) {
				$dmreVector[$k] = $dmreMatrix[$j][$i];
				$dmrePointer [$k] = $j;
				
			}
		}
		print $loghandle "\n" if $verbose;
		$k++;
	}
	open (MAX,">maxlevel.txt");
	for (my $k = 1; $k <@maxLevel; $k++) {
		if (defined $seqToStruct{$k}){
			print MAX $k."|".$seqToStruct{$k}."|".$jsmaxLevel[$k]."\n";
		}
	}
	close MAX;
#	open (MAX,">jsmaxlevel.txt");
#	for (my $k = 1; $k <@jsmaxLevel; $k++) {
#		print MAX $k."|".$seqToStruct{$k}."|".$jsmaxLevel[$k]."\n";
#	}
#	close MAX;




	#
	# Normalize the scores
	#
	my @consImportancez = normalize (\@consImportanceVector);
	my @consJSz = normalize (\@consJSVector);
	my @consREz = normalize (\@consREVector);

	for (my $i = 1; $i  < @consImportanceVector; $i++){
		$consImportanceMap{$i} = $consImportancez[$i];
		$consREMap{$i} = $consREz[$i];
		$consJSMap{$i} = $consJSz[$i];
		$dmreMap {$i}  = $dmreVector[$i] ;
	}

	$parameters{cons_importance_map} = \%consImportanceMap;
	$parameters{cons_re_map} = \%consREMap;
	$parameters{cons_js_map} = \%consJSMap;
	$parameters{dmre_map} = \%dmreMap;

	# Rank based on scores
	my @rankImportance = rank (\@consImportancez);
	my @rankJS = rank (\@consJSz);
	my @rankRE = rank (\@consREz);
	my @rankDMRE = rank (\@dmreVector);
	for (my $i = 1; $i  < @consImportanceVector; $i++){
		$rankImportanceMap{$i} = $rankImportance[$i];
		$rankREMap{$i} = $rankRE[$i];
		$rankJSMap{$i} = $rankJS[$i];
		$rankdmreMap{$i}  = $rankDMRE[$i] ;
	}
	$parameters{rank_cons_importance_map} = \%rankImportanceMap;
	$parameters{rank_cons_re_map} = \%rankREMap;
	$parameters{rank_cons_js_map} = \%rankJSMap;
	$parameters{rank_dmre_map} = \%rankdmreMap;

	%{$_[0]}	= %parameters;
}



sub findMaximum {
	my @mat = @{$_[0]};
	my $level = $_[1];
	my $length = $_[2];
	my $ignoreroot = $_[3];
	my @maxlevel  = ();

	my @vector = ();
	
	for (my $i=1; $i<= $length ; $i++){
		$vector [$i] = 0;
		foreach (my $j=1; $j<= $level; $j++){
			if ($j==1) { 
				$vector[$i] = $mat[$j][$i];
				$maxlevel[$i] = 1;
			} else {
				if ($mat[$j][$i] > $vector[$i]){
					$vector[$i] = $mat[$j][$i];
					$maxlevel[$i] = $j;
				}
			}
		}
		if ($ignoreroot && $maxlevel[$i]==1) {
			$vector[$i] = 0;
		}
	}
	@{$_[4]} = @vector;
	@{$_[5]} = @maxlevel;

}


sub getProfile{
	my %parameters 	= %{$_[0]};
	my $level 		= $_[1];
	my %validMap	= %{$parameters{valid_amino_acid}};
	my $loghandle 	= $parameters{loghandle};
	my $verbose		= $parameters{verbose3};
	my $sequence 	= $parameters{sequence};

	my @seqMatrix = @{$parameters {seq_matrix}};
	my %id2indexMap = %{$parameters {id2index_map}};
	my @subtreeRank = @{$parameters {subtree_rank}};

	my @profile = ();
	if ($verbose){ 
		print "Level = $level\n";
		print "subtreeRank\n";
		for (my $i = 0; $i < @seqMatrix; $i ++){
			print $i."\t".$subtreeRank[$i]."\n";
		}
	}

	
	for (my $i = 0; $i < @seqMatrix; $i ++){
		my @seq = split (//,$seqMatrix[$i]);
		if ($i==0){
			for (my $j=0 ; $j< @seq; $j++){
				foreach my $key (keys %validMap){
					${$profile[$j]}{$key} = 0;
				}
				${$profile[$j]}{$gapCharacter} = 0;
			}
		}
		if ($verbose){
			print $loghandle "SubtreeRank = $subtreeRank[$i]\n";
		}
		print "i=$i\tsubtreerank=".$subtreeRank[$i]."\tlevel=".$level."\n";
		next if ($subtreeRank[$i] != $level);
		if ($verbose){
			print $loghandle "Length = ".@seq."\n";
			print $loghandle join ("", @seq)."\n";
		}
		for ( my $j=0 ; $j< @seq; $j++){
			${$profile[$j]}{$seq[$j]}++ ;
		}
	}
	
	if ($parameters{verbose1}){
		print $loghandle "Profile length=".@profile."\n";
	}
	if ($verbose){
		my @tmp  = split (//, $sequence->seq());
		for (my $i=0; $i<@profile; $i++){
			my %map = %{$profile[$i]};
			print $loghandle "*************\n";
			print $loghandle $i." with residue ".$tmp[$i]."\n";
			foreach my $key (keys %map){
				print $loghandle $key."\t".$map{$key}."\n";
			}
		}		
	}
	@{$_[2]} = @profile;

}


sub globalConservation{
	my %parameters	= %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	my $sequenceID	= $parameters{sequence_id};
	my $verbose		= $parameters{verbose2};
	my %validMap	= %{$parameters{valid_amino_acid}};
	my $aln 		= $parameters{msa};
	my $seq			= $parameters{sequence};
	my %background  = %{$parameters{background_distribution}};
	if (defined $verbose){
		print "In globalConservation\n";
	}
	
	# The array that contains all the sequences 
	my @seqMatrix = ();
	# Maps the sequence ids to their indices in the seqMatrix
	my %id2indexMap = ();
	# The rank of each sequence w.r.t the tree trace
	my @subtreeRank = ();

	my %globalMap = ();
	my %globalreMap     = ();
	my %globaljsMap     = ();
	my %globalRank = ();
	my %globalRankRE     = ();
	my %globalRankJS    = ();
	my @seqArray = $aln->each_seq;

	my $level = 1;
	for my $seqObj (@seqArray){
		push (@seqMatrix, $seqObj->seq);
		$id2indexMap {$level} = $#seqMatrix;
		$subtreeRank [$#seqMatrix] = $level;
	}

	$parameters {seq_matrix} = \@seqMatrix;
	$parameters {id2index_map} = \%id2indexMap;
	# Every sequence is given a rank
	# which is the rank of the subtree in which the sequence 
	# appears first.
	# This is used while computing profiles
	# at different cuts.
	$parameters {subtree_rank} = \@subtreeRank;

	my @q = ();
	foreach my $key (keys %validMap ) { 
		push  (@q, $background{$key});
	}

	
	my @profile = ();
	getProfile (\%parameters, $level, \@profile);
	my $averageCons = 0;
	my @conservation = ();
	my @reConservation = ();
	my @jsConservation = ();

	my $total  = 0;
	for (my $i=0; $i<@profile; $i++){
		my $denom = 0;
		my $max = 0;
		my $maxKey;
		my %map = %{$profile[$i]};
		
		# Simple global conservation
		my $total = 0;
		foreach my $key (keys %validMap){
			if ($map{$key} > $max){
				$max = $map{$key};
				$maxKey = $key;
			}
			$denom += $map{$key};
		}
		$total++ if ( defined $maxKey && defined $validMap{$maxKey});
		$conservation[$i] 	= 0;
		$conservation[$i] 	= $max/$denom if ($denom > 0 && $maxKey ne $gapCharacter);
		$conservation[$i] *= (1 - $map{$gapCharacter}/($denom+$map{$gapCharacter}));
		$averageCons 	+= $conservation[$i];

		my @p = ();
		my @r = ();
		foreach my $key (keys %validMap ) { 
			if ($denom >0 ) {
				push (@p, $map{$key}/$denom);
				push (@r, 0.5 * $map{$key}/$denom + 0.5* $background {$key});
			} else { 
				push (@p, 0);
				push (@r, 0);
			}
		}
		$reConservation[$i] = relativeEntropy (\@p, \@q);
		$jsConservation[$i] = 0.5 * relativeEntropy(\@p,\@r) + 0.5 * relativeEntropy (\@q,\@r);
		$reConservation[$i] *= (1 - $map{$gapCharacter}/($denom+$map{$gapCharacter}));
		$jsConservation[$i] *= (1 - $map{$gapCharacter}/($denom+$map{$gapCharacter}));
	}
	$averageCons /= $total if ($total > 0);
	my @conservationz = normalize (\@conservation);


	if ($verbose){
		print $loghandle "averageCons=$averageCons\n";
		print $loghandle arrayToString ("%5.3f", \@conservation, ",")."\n";
	}	

	
	my @rankCons = rank (\@conservationz);
	my @rankRE = rank (\@reConservation);
	my @rankJS = rank (\@jsConservation);
	for (my $i=0, my $k=1; $i< @conservationz ; $i++){
		next if (!defined ($validMap{$seq->subseq($i+1,$i+1)}));
		$globalMap{$k} 		= $conservationz[$i];
		$globalreMap {$k} 	= $reConservation [$i];
		$globaljsMap {$k}   = $jsConservation [$i];
		if ($verbose){
			print $loghandle $k."\t".$seq->subseq($i+1,$i+1)."\n";
		}

		$globalRank{$k} = $rankCons[$i];
		$globalRankRE{$k} = $rankRE[$i];
		$globalRankJS{$k} = $rankJS[$i];

		$k++;
	}
	$parameters{rank_global_map} = \%globalRank;
	$parameters{rank_global_re_map} = \%globalRankRE;
	$parameters{rank_global_js_map} = \%globalRankJS;

	$parameters{global_map} = \%globalMap;
	$parameters{global_re_map} = \%globalreMap;
	$parameters{global_js_map} = \%globaljsMap;
	%{$_[0]} = %parameters;
}

sub getMin{
#	my $min=shift @_;
	my $min;
	my $minIndex=0;
	my $index=0;
	while (@_){
		my $tmp = shift @_;
		next if !defined $tmp;
		$index++;
		if (defined $min){
			if( $tmp < $min ){
				$min = $tmp;
				$minIndex = $index;
			}
		}else{
			$min=$tmp;
			$minIndex = $index;
		}
	}
	return ($min, $minIndex);
}

sub getMax{
#	my $max=shift @_;
	my $max;
	my $maxIndex=0;
	my $index=0;
	while (@_){
		my $tmp = shift @_;
		next if !defined $tmp;
		$index++;
		if (defined $max){
			if( $tmp > $max ){
				$max = $tmp;
				$maxIndex = $index;
			}
		}else{
			$max = $tmp;
			$maxIndex = $index;
		}
	}
	return ($max, $maxIndex);

}


sub initParameters{
	my %parameter=%{$_[0]};

	$parameters{pic}				=50;
	$parameters{no_cluster}			= 0;
	$parameters{verbose1}			= 1;
	$parameters{verbose2}			= 0;
	$parameters{verbose3}			= 0;	
	$parameters{cluster_distance}	= 5;
	$parameters{cluster_cutoff}		= 0.6;
	$parameters{gap_threshold}		= 1;
	$parameters{entropy}			= 0;
	$parameters{weight_function}	= "step";
	$parameters{normalize}			= 0;

	$parameters{aa_file}		= "aa.txt";
	$parameters{dmm_file} 		= "Uniform.plib";
	$parameters{background_distribution_file} 		= "blosum62.distribution";
	$parameters{score_matrix_file}	= "blosum62.txt";
	$parameters{area_file}			= "area.txt";
	$parameters{benchmark}		= 0;
	$parameters{data_dir}		= "/home/sriram_s/bpgtools/scripts/data";
	$parameters{data_dir}		= $ENV{'DATA_DIR'} if defined $ENV{'DATA_DIR'};

	my @names 			 		= (	"cons_importance_map"	,
									"cons_js_map"		,
									"cons_re_map" 	,
									"global_map"			,
									"global_js_map"			,
									"global_re_map"			);
#	my @names 			 		= (	"cons_js_map"		);
	
	$parameters{map_names} = \@names;


	if (@ARGV < 1){
		printUsage ();
	}

	my %tmpMap = readMap ($ARGV[0]);
	foreach my $key (keys %tmpMap){
		$parameters{$key} = $tmpMap{$key};
	}
	

	$parameters{verbose2}	= 1 if $parameters{verbose3}==1;
	$parameters{verbose1}	= 1 if $parameters{verbose2}==1;
	$parameters{data_dir}	= "." if !defined $parameters{data_dir};
	my $prefix = "";
	my $indir  = $parameters{indir};
	$prefix = $indir."/" if $indir;
	
	# 1 if a PDB id is given. 0 otherwise.
	# 
	$parameters { only_sequence } = 0;

	if (	defined $parameters{msa_file}
			&& defined $parameters{tree_file} 
			&& defined $parameters{structure_id} ){
		
			my ($structId, $chainId)
				= split (/:/,$parameters{structure_id});
			$parameters{structure_id}=$structId;
			if ( !defined $chainId ){
				$parameters{chain} = "default" ;
			} else {
				$parameters{chain} = $chainId;
			}	
			if (!defined $parameters {sequence_id}){
				$parameters{sequence_id} = $parameters{structure_id};
			}	

			$parameters{msa_file} = $prefix.$parameters{msa_file};
			$parameters{tree_file}= $prefix.$parameters{tree_file};


			if ($parameters {verbose1}){
				print "MSA file:". $parameters{msa_file}."\n";
				print "Tree file:".$parameters{tree_file}."\n";
				print "Sequence ID:".$parameters{sequence_id}."\n";
				print "PDB ID:".$parameters{structure_id}."\n";
				print "Chain ID:".$parameters{chain}."\n";
				
			}
	} elsif (	defined $parameters{msa_file}
			&& defined $parameters{tree_file} 
			&& defined $parameters{sequence_id} ){
			$parameters {only_sequence} = 1;
			$parameters{msa_file} = $prefix.$parameters{msa_file};
			$parameters{tree_file}= $prefix.$parameters{tree_file};
			$parameters{structure_id}  = $parameters{sequence_id};
			$parameters{chain} = "default" ;
	} else {
		printUsage();	
	}



	$parameters{aa_file}	= $parameters{data_dir}."/".$parameters{aa_file};
	$parameters{score_matrix_file}	= $parameters{data_dir}."/".$parameters{score_matrix_file};
	$parameters{area_file} = $parameters{data_dir}."/".$parameters{area_file};
	$parameters{dmm_file}	= $parameters{data_dir}."/".$parameters{dmm_file};
	$parameters{background_distribution_file} 		= $parameters{data_dir}."/".$parameters{background_distribution_file};
	
	my %aaMap	= readMap($parameters{aa_file});
	$parameters{amino_acid_map} = \%aaMap;

	read_dmmfile ($parameters{dmm_file});
	read_distribution  ($parameters {background_distribution_file});

	
	$parameters{pdb_file}		= $prefix.$parameters{structure_id}."\.pdb";
	my $logfile					= "";

	$logfile					= $parameters{structure_id}."-".$parameters{chain}."\.log";
	$logfile					= $parameters{outdir}."/".$logfile if defined $parameters{outdir};
	$parameters{logfile}		= $logfile;
	print "Log file=$logfile\n";
	open (LOGFILE, ">$logfile") or die "Error: Unable to open log file:$logfile\n";
	$parameters{loghandle}		= *LOGFILE;

	print LOGFILE "PIC=".$parameters{pic}."\n";
	print LOGFILE "PDB file = ".$parameters{pdb_file}."\n";
	print LOGFILE "Chain = ".$parameters{chain}."\n";
	print LOGFILE "MSA file=".$parameters{msa_file}."\n";

	my $sequences	= 	Bio::AlignIO->new (	-file 	=> $parameters{msa_file}, 
											-format	=> 'fasta');
	my $aln = $sequences->next_aln();
	my $sequence;

	$aln->map_chars ('\.','-');
	$aln->map_chars ('\*','-');
 	foreach my $seq ($aln->each_seq()){
 		my $id = $seq->display_id();
 		if ( index ($id,$parameters{sequence_id}) >=0 ){
 			 $sequence = $seq;
 			 last;
 		}
	}


	if (!defined $sequence){
		die "Error: ".$parameters{sequence_id}." not present in MSA ".$parameters{msa_file}."\n";
	}

	$parameters{msa}		= $aln;
	$parameters{sequence}	= $sequence;

	$sequences	= 	Bio::AlignIO->new (	-file 	=> $parameters{msa_file}, 
											-format	=> 'fasta');
	$aln = $sequences->next_aln();
	my $originalSequence;
 	foreach my $seq ($aln->each_seq()){
 		my $id = $seq->display_id();
 		if ( index ($id,$parameters{sequence_id}) >=0 ){
 			 $originalSequence = $seq;
 			 last;
 		}
	}
	$parameters{original_sequence} = $originalSequence;

	my @seqArray;
	unshift(@seqArray, $sequence);
	print LOGFILE "Original MSA Sequence\n".$originalSequence->seq()."\n";
	print LOGFILE "MSA Sequence\n".$sequence->seq()."\n";

	my $seq =  "";
	if ($parameters {only_sequence}==0 ) {
		my $structio 	= Bio::Structure::IO -> new (-file 	=> $parameters{pdb_file}, -format=>'pdb');
		my $struct 		= $structio -> next_structure;
		$parameters{structure}	= $struct;

		my $structurePrimary = "";
		my $flag = 0;
		for my $chain ($struct -> get_chains ){
			my $chainid = $chain->id;
			if($chainid eq $parameters{chain}){
				$flag = 1;
				for my $res ($struct->get_residues($chain)){
					my $resid  = $res->id;
					my ($res,$pos)=split (/-/,$resid);
# Checking for strange residues or HETATM records.
					next if !defined $aaMap{$res};
					$structurePrimary=$structurePrimary.$aaMap{$res};
				}
			}
		}
		print LOGFILE "PDB sequence\n".$structurePrimary."\n";

		if (!$flag){
			die "Error: PDB file ".$parameters{pdb_file}." does not contain chain ".$parameters{chain}."\n";
		}


		$seq = Bio::Seq->new (-seq => $structurePrimary, -id=>'STRUCTURE:'.$parameters{structure_id}.":".$parameters{chain});
		$parameters{structure_primary}=$seq;

	} else { 
		$seq = Bio::Seq->new (-seq => $originalSequence->seq(), -id=>'STRUCTURE:'.$parameters{structure_id}.":".$parameters{chain});
		$parameters{structure_primary}=$seq;
	}

	my $treeio=new Bio::TreeIO( -format=>"newick", -file=>$parameters{tree_file});
	my $tree=$treeio->next_tree;
	$parameters{tree} = $tree;


	my $root=$tree->get_root_node;
	my $designatedNode;
	print LOGFILE "Tree\n";

	foreach my $node ($root->get_Descendents()){
		if ($node->is_Leaf && (index ($node->id, $parameters{sequence_id})>=0)){
			$designatedNode = $node;
			last;
		}
		print LOGFILE $node->id.",".$node->branch_length()."\n";
	}

	if (!defined $designatedNode){
		die " Error:".$parameters{tree_file}. " does not contain ".$parameters{sequence_id}."\n";
	}
	
	# Create an array of Bio::Seq objects containing the two sequences.
	# Align using Clustalw.
	unshift(@seqArray, $seq);
	$parameters{seq_array_ref} = \@seqArray;
	
	%{$_[0]}	= %parameters;
}


sub createMaps{
	my %parameters 	= %{$_[0]};
	my $loghandle 	= $parameters{loghandle};
	print $loghandle  "Creating data structures\n";
	my $msaFile = $parameters {msa_file};
	my $pdbFile = $parameters {pdb_file};
	my $chain 	= $parameters {chain};
	my %validMap= %{$parameters {valid_amino_acid}};
	my %aaMap	= %{$parameters {amino_acid_map}};
	my %auxMap 	= ();
	my $verbose	= $parameters{verbose2};

	my %posMap	= ();
	
	my $sequence	= $parameters{sequence};
	my $struct		= $parameters{structure};
#	print "Struct=".$struct."\n";

	# If no PDB id is given, then all the mappings are trivial
	if  ( $parameters{only_sequence} ) {
		my %seqToStruct	= ();
		my %structToSeq	= ();
		my $sequence;
		my $structure	= "";
		my $seqArrayRef	= $parameters{seq_array_ref};

		my @params 		= ('ktuple'=>2, 'matrix' =>'BLOSUM');
		my $factory 	= Bio::Tools::Run::Alignment::Clustalw->new(@params);

		my $aln = $factory->align($seqArrayRef);
		$parameters{seq_struct_align}	= $aln;
		print $loghandle  "Alignment\n";
		foreach my $seq ($aln->each_seq){
			print $loghandle ">".$seq->display_id()."\n";
		}


		foreach my $seq ($aln->each_seq){
			if ($seq->id =~/STRUCTURE/){
				$structure = $seq;
			}else {
				$sequence = $seq;
			}
		}
		my $icode = "";
		for (my $i=1; $i<=$aln->length() ; $i++){
			my $sequenceResidue = $sequence->subseq($i,$i);
			
			$posMap {$i} = $i;
			$auxMap {$i} = "$chain|$sequenceResidue|$icode";
			$seqToStruct{$i}=$i;
			$structToSeq{$i}=$i;
		}

		$parameters{pos_map}=\%posMap;
		$parameters{aux_map}= \%auxMap;
		$parameters{seq_to_struct}=\%seqToStruct;
		$parameters{struct_to_seq}=\%structToSeq;

	} else {

		my $index=0;
		for my $chain ($struct -> get_chains ){
			my $chainid = $chain->id;
			if($chainid eq $parameters{chain}){
				for my $res ($struct->get_residues($chain)){
					my $resid  = $res->id;
					my ($tmpres,$pos) = split (/-/,$resid);
# Checking for strange residues or HETATM records.
					next if !defined $aaMap{$tmpres};
					my @atoms = $struct->get_atoms($res);
					my $icode = $atoms[0]->icode();
					$icode ="" if (!defined $icode);

					$index++;
					$posMap{$index} = $pos;

					if (!defined $parameters{map}){
						$auxMap{$pos}="$chainid|$tmpres|$icode";
					}

				}
			}
		}
		if ($verbose){
			print $loghandle "posMap:index->position in PDB\n";	
			foreach my $key (sort{$a<=>$b} keys %posMap){
				print $loghandle $key."\t".$posMap{$key}."\n" if defined $posMap{$key};
			}
		}
		$parameters{pos_map}=\%posMap;
		$parameters{aux_map}= \%auxMap;

		if (defined $parameters{map}){
			readSeqStructMaps (\%parameters);
		} else {
			createSeqStructMaps (\%parameters);
		}
	}

	my $originalSequence= $parameters{original_sequence};
	my $originalPos = 0;
	my $pos = 0;
	my %originalMap=();
	for (my $i=1; $i <= $originalSequence->length(); $i++){
		my $originalSubseq = $originalSequence->subseq($i,$i);
		print $loghandle  "$i:$originalSubseq\n" if $verbose;
		my $subseq = $sequence->subseq($i,$i);
		if (defined $validMap{$originalSubseq}){
			$pos++;
			print $loghandle "pos=$pos, $originalSubseq\n" if $verbose ;
		}
#if (defined $validMap{$originalSubseq} || $originalSubseq=~/[a-z]/ || $originalSubseq eq "X"){
		if (defined $validMap{$originalSubseq} || $originalSubseq eq "-" || $originalSubseq eq "X"){
			$originalPos++;
			$originalMap{$originalPos}=$pos;
			print $loghandle "originalPos=$originalPos, $originalSubseq\n" if $verbose ;
			print $loghandle "originalMap:$originalPos->$pos\n" if $verbose;
		}
	}

	if ($verbose){
		print $loghandle "originalMap:position in MSA (in a2m)-> position in MSA (in FASTA)\n";
		foreach my $key (sort {$a<=>$b} keys %originalMap){
			print $loghandle $key."\t".$originalMap{$key}."\n";
		}
	}	
	$parameters{original_map}  = \%originalMap;

	%{$_[0]}	= %parameters;
}


sub createSeqStructMaps{
	my %parameters 	= %{$_[0]};
	my $aln 		= $parameters{seq_struct_align};
	my %validMap	= %{$parameters {valid_amino_acid}};
	my $verbose		= $parameters{verbose1};
	my $loghandle	= $parameters{loghandle};
	my %posMap 			= %{$parameters{pos_map}};
	my %seqToStruct	= ();
	my %structToSeq	= ();
	my $sequence;
	my $structure	= "";
	my $seqArrayRef	= $parameters{seq_array_ref};
	
	my @params 		= ('ktuple'=>2, 'matrix' =>'BLOSUM');
	my $factory 	= Bio::Tools::Run::Alignment::Clustalw->new(@params);

	$aln = $factory->align($seqArrayRef);
	$parameters{seq_struct_align}	= $aln;
	print $loghandle  "Alignment\n";
	foreach my $seq ($aln->each_seq){
		print $loghandle ">".$seq->display_id()."\n";
		print $loghandle $seq->seq;
	}


	foreach my $seq ($aln->each_seq){
		if ($seq->id =~/STRUCTURE/){
			$structure = $seq;
		}else {
			$sequence = $seq;
		}
	}
	my $seqPos = 0;
	my $structPos=0;
	for (my $i=1; $i<=$aln->length() ; $i++){
		my $structResidue = $structure->subseq($i,$i);
		my $sequenceResidue = $sequence->subseq($i,$i);
		if ($parameters{verbose1}){
			print $loghandle "i=$i,struct=$structResidue,seq=$sequenceResidue\n";
		}
		if ( defined $validMap{$sequenceResidue}){
			$seqPos++;
		}
		if ( defined $validMap{$structResidue}){
			$structPos++;
		}
		if ( defined $validMap{$sequenceResidue} && defined $validMap{$structResidue}){
			my $pdbPos = $posMap {$structPos};
			$seqToStruct{$seqPos}=$pdbPos;
			$structToSeq{$pdbPos}=$seqPos;
			if ($parameters{verbose1}){
				print $loghandle "Inserting seqPos=$seqPos, structPos=$structPos\n";
			}
		}
	}

	if ($verbose){
		print  $loghandle "seqToStruct\n";
		foreach my $pos (sort {$a<=>$b} keys %seqToStruct){
			print $loghandle $pos."\t".$seqToStruct{$pos}."\n" if defined $seqToStruct{$pos};
		}
		print $loghandle "StructToSeq\n";
		foreach my $pos (sort {$a<=>$b} keys %structToSeq){
			print $loghandle $pos."\t".$structToSeq{$pos}."\n" if defined $structToSeq{$pos};
		}
	}
	$parameters{seq_to_struct}=\%seqToStruct;
	$parameters{struct_to_seq}=\%structToSeq;
	%{$_[0]}	= %parameters;
	print "Exiting createSeqStructMaps\n";
}

sub readSeqStructMaps{
	my %parameters = %{$_[0]};
	my $file = $parameters {map};
	my $loghandle = $parameters{loghandle};

	my %seqToStruct;
	my %structToSeq;
	my %auxMap;

	open (FILE, $file);
	while  (<FILE>){
		chomp;
		my @tmp = split(/\|/);
		my $seqpos = $tmp[0];
		my $pdbpos = $tmp[3];
		$seqpos =~ s/\s+//g;
		$pdbpos =~ s/\s+//g;
		if (@tmp < 5){
			die "Error:Bad line in map file $file in line number ".$.."\n";
		}
		$seqToStruct{$seqpos} = $pdbpos;
		$structToSeq{$pdbpos} = $seqpos;
		$auxMap {$pdbpos} = $tmp[1]."|".$tmp[2]."|".$tmp[4];
	}
	
	close FILE;
	print $loghandle "seqToStruct\n";
	foreach my $key (sort {$a<=>$b} keys %seqToStruct){
		print $loghandle  $key."\t".$seqToStruct{$key}."\n";
	}
	$parameters{seq_to_struct}=\%seqToStruct;
	$parameters{struct_to_seq}=\%structToSeq;
	$parameters{aux_map}= \%auxMap;
	%{$_[0]}	= %parameters;
}


sub mapSequenceToScore{
	my %parameters	= %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	print $loghandle "Mapping the sequence and the structure\n";
	my $sequenceID	= $parameters{sequence_id};
	my $msaFile	 	= $parameters{msa_file};
	my $gapThreshold= $parameters{gap_threshold};
	my %validMap	= %{$parameters{valid_amino_acid}};
	my $verbose		= $parameters{verbose2};
	my %scoreMap	= %{$parameters{score_matrix}};
	
	my %conservationMap 		= ();


	# Retrieve the sequence from the MSA file
	my $sequences	= 	Bio::AlignIO->new (	-file 	=> $msaFile, 
											-format	=> 'fasta');
	my $aln=$sequences->next_aln();
	my $length=$aln->length();
	
	my $sequencePrimary = $parameters{sequence};
	if ($verbose){
		foreach my $seq ($aln->each_seq){
			print $loghandle "Seq ID=".$seq->id."\n";
		}

	}

	my $gapsAllowed = 0;
	$gapsAllowed=$gapThreshold * $aln->no_sequences;
	# Do the scoring
	my $position = 0;
	for (my $i = 1; $i <= $length ; $i++){
		my $slice = slice($aln,$i,$i);
		my @column=$slice->each_seq;
		my $gapCount = 0;
		my $score = 0;
		my $weight = 0;
		
		my $reference = $sequencePrimary->subseq($i, $i);
		if ( defined $validMap{$reference}){
			$position++;
		}else{
			next;
		}
		if ($verbose){
			print $loghandle "Position=$position, Reference=$reference\n";
		}
		my ($min, $max) = split(/:/,$scoreMap{$reference});
		for (my $j = 0; $j<@column; $j++){
			my $c1=$column[$j]->seq();
 			if (! defined $validMap{$c1}) {
 				$gapCount++;
 			} else {
 				$score+= ($scoreMap{$reference.$c1}-$min)/($max-$min);			
 				$weight ++;
 			}

			last if ($gapCount > $gapsAllowed);
		}
		$score-= 1;
		if ($gapCount <= $gapsAllowed){
			print $loghandle "Setting position $position\n" if ($verbose);
			$conservationMap{$position}=$score/$weight ;
		}else{
			print $loghandle "Avoiding position\n" if ($verbose);
		}
	}
	if ($verbose){
		print  $loghandle "Column score map\n";
		foreach my $key (sort {$a<=> $b} keys %conservationMap){
			print $loghandle $key."\t".$conservationMap{$key}."\n";
		}
	}
	$parameters{conservation_map}=\%conservationMap;
	%{$_[0]}	= %parameters;
}


sub initMaps{
	my %parameters = %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	print "Initializing data structures\n";
	print $loghandle "Initializing data structures\n";

	# The map of valid amino acids
	my %validMap=();
	my %aaMap	= %{$parameters{amino_acid_map}};
	foreach my $value (values %aaMap){
		$validMap{$value}=$value;
	}

	# Map of pairs of amino acids to scores.
	# scores here are the BLOSUM62 scores.
	my %scoreMap = initScoreMap (\%parameters);

	# Map of the surface area of amino acids
	my %areaMap = readMap ($parameters{area_file});
	$parameters{area_map} = \%areaMap;
	
	$parameters{amino_acid_map}	= \%aaMap;
	$parameters{valid_amino_acid} = \%validMap;
	$parameters{score_matrix}	= \%scoreMap;

	# Create other maps
	createMaps (\%parameters);
	mapSequenceToScore (\%parameters);
	%{$_[0]} = %parameters;
}


sub listSum{
	my $sum = 0;
	my $count = 0;
	foreach my $a (@_){
		next if (!defined $a);
		$sum+=$a;
		$count ++;
	}
	return ($sum, $count);
}

sub weightedSum{
	my %values = %{$_[0]};
	my %weights = %{$_[1]};

	my $sum=0;
	my $count=0;
	foreach my $key (keys %values){
		if (defined $weights{$key}){
			$sum+=$values{$key}*$weights{$key};
			$count+=$weights{$key};
		}
	}
	return ($sum, $count);
}



# readMap	- Loads a map of key-value pairs from a file.
#				Each key-value pair must be separated by newlines.
#				The key must be separated from the value by whitespaces.
# Param		- The file containing the key-value pairs.
# Return	- The map
sub readMap{
	my $mapFile = shift (@_);
	my %map = ();
	open (FILE,$mapFile);	
	my @lines = <FILE>;
	foreach my $line (@lines){
		chomp ($line);
		my @tmp = split (/\s+/,$line);
		$map {$tmp[0]} = $tmp[1];
	}
	close FILE;
	return %map;
}

sub printUsage{
	print "Usage:$0 <config file>\n";
    die;
}




sub done{
	my %parameters = %{$_[0]};
	my $loghandle	= $parameters{loghandle};
	close $loghandle;
}

sub arrayToString{
	my $fmt 	= shift @_;
	my @array 	= @{shift @_};
	my $delim	= shift @_;
	$delim 		= "\t" if !defined $delim;
	my @output 	= ();
	foreach my $i (@array){
		push (@output, sprintf ($fmt,$i));
	}
	return join($delim, @output);
}



sub entropy {
	my @p = @{$_[0]};
	my @q = split (//,"1"x@p);
	my $loghandle = $parameters{loghandle};
	my $re = relativeEntropy (\@p,\@q);
	if ($parameters{verbose3}){
		print $loghandle "p=".join(",",@p)."\n";
		print $loghandle "q=".join(",",@q)."\n";
		print $loghandle "relative entropy=". $re."\n" ;
	}

	return  log(@q)/log(2) - $re ;
	
}

sub relativeEntropy {
	my @p = @{$_[0]};
	my @q = @{$_[1]};
	my $loghandle = $parameters{loghandle};

	if ($parameters{verbose3}){
		print $loghandle "p=@p\n";
		print $loghandle "q=@q\n";
	}

	my $zp = 0;
	my $zq = 0;
	for (my $i = 0 ; $i<@p ; $i++){
		$zp += $p[$i];
	}
	for (my $i = 0 ; $i<@q ; $i++){
		$zq += $q[$i];
	}


	die "Error in $0-relativeEntropy: lengths of arrays p and q differ\n" if (@p != @q);

	my $re =  0;
	for (my $i = 0 ; $i<@p ; $i++){
		my $p = 0;
		my $q = 0;
		$p = $p[$i]/$zp if ($zp > 0);
		$q = $q[$i]/$zq if ($zq > 0);
		$re += $p * log ($p/$q)/log (2) if ($p > 0 && $q > 0);
	}
	return $re;
}


sub normalize {
	my $loghandle = $parameters{loghandle};
	my @r = @{$_[0]};
	return if (@r == 0);
	# The input array may have undefined entries.
	# This might be the case if the array were to be indexed from 1.
	# Defined a new array @p such that all the elements in @p are defined and operate on this array.
	my @p = ();
	foreach my $r (@r){
		next if (!defined $r);
		push (@p, $r);
	}

	# Use the count in computing averages as a safeguard againsts undefined elements.
	my ($ave, $count) = listSum (@p);
	$ave = $ave/$count if ($count >0);
	if ($parameters {verbose3}){
		print $loghandle "Normalizing ".join(",",@p)."\n";
	}

	my @q = ();
	foreach my $p (@p){
		push (@q, $p*$p);
	}
	my ($sqave, $sqcount) = listSum (@q);
	$sqave = $sqave/$sqcount;
	my $stdev = sqrt($sqave - $ave**2);

	if ($parameters{verbose3}){
		print $loghandle "Ave=".$ave."\n";
		print $loghandle "Stdev=".$stdev."\n";
	}

	# The output array retains the undefined elements in the input array
	my @z = ();
	for (my $i = 0; $i < @r; $i ++){
		next if (!defined $r[$i]);
		$z[$i] = 0;
		$z[$i] = ($r[$i] - $ave)/$stdev if ($stdev > 0);
	}
	return @z;
}


sub correlation {
	my @pz = normalize ($_[0]);
	my @qz = normalize ($_[1]);

	die "Error in $0-correlation:Inconsistent length arrays \n" if (@pz!=@qz);
	my $r = 0;
	my $n = 0;
	for (my $i = 0; $i < @pz; $i++){
		next if (!defined $pz[$i] && !defined $qz[$i]);
		die "Error in $0-correlation: Inconsistent elements in pz" if (!defined $pz[$i]);
		die "Error in $0-correlation: Inconsistent elements in qz" if (!defined $qz[$i]);
		$r += $pz[$i]*$qz[$i];
		$n++;
	}
	$r /= $n;
	return $r;
}

sub correlate {
	my $loghandle = $parameters{loghandle};
	my %p = %{$_[0]};
	my %q =%{ $_[1]};
	my @p = ();
	my @q = ();
	foreach my $key (keys %p){
#		print $loghandle  "$key\t$p{$key}\t$q{$key}\n";
		next  if (! defined $q{$key});
		push (@p, $p{$key});
		push (@q, $q{$key});
	}

	my $r = correlation (\@p, \@q);
	return $r;
}



sub read_dmmfile {
	my $file = $_[0];

	my $dmm = Dmm->new ();

	open (FILE, $file) or die "Error: unable to open file $file";
	
	my $index = 0;
	my @tmparray = ();
	while (<FILE>){
		chomp;
		if (/NumDistr/){
			$_ =~s/\s//g;
			my @tmp  = split (/=/);
#			print $tmp[1]."\n";
			$dmm->component($tmp[1]); 
		}
		if (/Number/){
			$_ =~s/\s//g;
			my @tmp  = split (/=/);
			$index  = $tmp[1]; 
#			print "index=$index\n";
		}
		
		if (/Order/){
			my @tmp = split (/\s+/);
			shift @tmp;

			for (my $i = 0; $i < @tmp; $i++){
				$dmm->aamap($i, $tmp[$i]);
			}
		}
		
		if (/Mixture/){
			$_ =~s/\s//g;
			my @tmp  = split (/=/);
			$dmm->weights ($index, $tmp[1]); 
#			print "Inserting ".$tmp[1]." at $number\n";
		}
		if (/Alpha/){
			my @tmp = split (/=/);
			$tmp[1] =~ s/^\s+//g;
			my @tmp1 = split (/\s+/,$tmp[1]);
			shift (@tmp1);
			my $total = 0;
			print "DMM=".join (",",@tmp1)."\n";
			for (my $i= 0; $i < @tmp1; $i++){
				$tmp1[$i] *= 20;
#				$tmp1[$i] *= 1;
				
				$tmparray[$index][$i] = $tmp1[$i];
				$total += $tmp1[$i];
			}
			$dmm->totalalpha ($index, $total);
		}
	}
	$dmm->alpha (\@tmparray);
	close FILE;

	$parameters{dmm_class} = $dmm;
}


sub read_distribution  {
	my $f = $_[0];
	print "F=$f\n";
	open (F,$f);
	my @mapkeys = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V');

	my @tmp  =();
	while (<F> ) {
		chomp;
		next if (/^#/);
		@tmp = split (/\s+/);
		last;
	}

	close F;

	my %map = ();
	for (my $i = 0; $i < @tmp; $i++) { 
		$map{$mapkeys[$i]} = $tmp[$i];
	}
	$parameters{background_distribution} = \%map;
	for (my $i = 0; $i < @tmp; $i++) { 
		print $mapkeys[$i]." = ". $tmp[$i]."\n";
	}
}


sub rank {
	my @p = @{$_[0]};
	return if (@p == 0);
	# The input array may have undefined entries.
	# This might be the case if the array were to be indexed from 1.
	# Defined a new array @p such that all the elements in @p are defined and operate on this array.
	my @r = ();
	foreach my $p (@p){
		next if (!defined $p);
		push (@r, $p);
	}

	my @ind = ();
	my @interrank = ();
	for (my $i = 0 ;  $i < @r; $i++){
		$ind[$i] = $i;
		$interrank[$i] = 0;
	}

	for (my $i = 0 ;  $i < @r; $i++){
		for (my $j = $i + 1; $j < @r; $j++){
			if ( $r[$i] < $r[$j]){
				my $tmp = $r[$i];
				$r[$i] = $r[$j];
				$r[$j] = $tmp;

				$tmp  = $ind[$i];
				$ind[$i] = $ind[$j];
				$ind[$j] = $tmp;
			}
		}
	}

	my $tmprank = 1;
	my $tmpscore = $r[0];
	my $tiedrank = 0 ;
	for (my $i = 0 ;  $i < @r; $i++){
		if ($r[$i] < $tmpscore ) { 
			$tmprank += $tiedrank;
			$tmpscore = $r[$i];
			$interrank[$ind[$i]] = $tmprank;
			$tiedrank = 1 ;
		} else {
			$tiedrank ++;
			$interrank[$ind[$i]] = $tmprank;
		}
	}

	# Output array contains the undefined elements in the input
	my @rank =  ();
	my $j = 0 ;
	for (my $i = 0 ; $i < @p ;$i++){
		next if (!defined $p[$i]);
		$rank[$i] = $interrank[$j];
		$j++;
	}


	return @rank;
}
