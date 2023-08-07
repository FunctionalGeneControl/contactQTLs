#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;

my %restriction;
read_restriction_fragments();

find_potential_probes();

open (OUT,'>','probe_positions.txt') or die $!;
open (FAILED,'>','failed_fragments.txt') or die $!;

foreach my $chr (sort keys %restriction) {

  my @fragments = @{$restriction{$chr}};

#  print Dumper(\@keepers);

  my $good_for_probes = 0;
  my $good_back_probes = 0;

  foreach my $fragment (@fragments) {

			my $found_probe = 0;

			if (exists $fragment->{front_probe_offset}) {
					print OUT join("\t",($chr,$fragment->{start}+$fragment->{front_probe_offset},$fragment->{start}+$fragment->{front_probe_offset}+119,'+',$fragment->{front_probe_sequence},$fragment->{start},$fragment->{end})),"\n";
					++$good_for_probes;
					$found_probe = 1;
			}
			if (exists $fragment->{back_probe_offset}) {
					print OUT join("\t",($chr,$fragment->{start}+$fragment->{back_probe_offset},$fragment->{start}+$fragment->{back_probe_offset}+119,'-',$fragment->{back_probe_sequence},$fragment->{start},$fragment->{end})),"\n";
					++$good_back_probes;
					$found_probe = 1;
			}

			unless ($found_probe) {
					print FAILED join("\t",($chr,$fragment->{start},$fragment->{end})),"\n";
			}

  }

  warn "For Chr $chr there were ".(scalar @fragments)." fragments ($good_for_probes for, $good_back_probes back)\n";
}


sub find_potential_probes {

		foreach my $chr (keys %restriction) {

				warn "Finding suitable probes for chr $chr\n";

				foreach my $fragment (@{$restriction{$chr}}) {

						my $front_pos;
						my $front_seq;
						my $back_pos;
						my $back_seq;

						my $seq = $fragment->{sequence};
						my $masked_seq = $fragment->{masked};

						# Do a stupid check to make sure that we're looking at the same sequence
						#for my $i (0..(length($seq)-1)) {
						#		my $m = substr($masked_seq,$i,1);
						#		if ($m ne 'N') {
						#				if ($m ne substr($seq,$i,1)) {
						#						die "Not equivalent\n$seq\n$masked_seq\n";
						#				}
						#		}
						#}


						for my $index (0..length($seq)-121) {
	
								my $subseq = substr($seq,$index,120);
								my $masked = substr($masked_seq,$index,120);
								next if ($masked =~ /N{3}/);

								my $gc = ($subseq=~tr/CGcg//);
								my $at = ($subseq=~tr/ATat//);

								if ($gc+$at == 0) {
										die "No sequence from\n$subseq\nmasked is\n$masked full is\n$seq\nfull masked is\n$masked_seq\nchr=".$fragment->{chr}." start=".($fragment->{start}+$index)." end=".($fragment->{start}+$index+120)."\n";;
								}
								
								my $percent_gc = ($gc/($gc+$at))*100;

								next if ($percent_gc < 25);
								next if ($percent_gc > 65);
								
								unless (defined $front_pos) {
										$front_pos = $index;
										$front_seq = $subseq;
										$index +=60;
								}
								else {
										$back_pos = $index;
										$back_seq = $subseq;
										$back_seq = reverse($back_seq);
										$back_seq =~ tr/GATCgatc/CTAGctag/;
								}
						}

						# The sonnicated fragments are size restricted to <500bp, so if our
						# probe doesn't sit within this region we don't want to keep it.

						if (defined $front_pos and $front_pos <= 380) {
								$fragment->{front_probe_offset} = $front_pos;
								$fragment->{front_probe_sequence} = $front_seq;
						}

						if (defined $back_pos and $back_pos >= ($fragment->{end}-$fragment->{start})-500) {
								$fragment->{back_probe_offset} = $back_pos;
								$fragment->{back_probe_sequence} = $back_seq;
						}
				}
		}
}


sub read_masked_sequence {

		my ($chr) = @_;

		return read_fasta_file("/rds/general/user/hrayjone/projects/lms-spivakov-analysis/live/HiCUP_settings/GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.${chr}.fa");
}

sub read_restriction_fragments {

  my @files = </rds/general/user/hrayjone/projects/lms-spivakov-analysis/live/HiCUP_settings/GRCh38/Homo_sapiens.GRCh38.dna.chromosome.*.fa>;

  open (FRAG,'>','all_restriction_fragments.txt') or die $!;

  foreach my $file (@files) {

    my $chr = (split(/\./,$file))[-2];

		#next unless ($chr eq '18');

    my $masked = read_masked_sequence ($chr);

    warn "Getting fragments from chr $chr\n";

    my $sequence = read_fasta_file($file);

    my $last_pos = 0;

    while ($sequence =~ /GATC/ig) {
				my $pos = pos($sequence)-4;
				if ($last_pos) {

						my $masked_seq = substr($masked,$last_pos-1,($pos-$last_pos));
						my $seq = substr($sequence,$last_pos-1,($pos-$last_pos));

						print FRAG join("\t",($chr,$last_pos,$pos)),"\n";
						push @{$restriction{$chr}},{chr=>$chr,start=>$last_pos,end=>$pos,sequence=>$seq,masked=>$masked_seq};
				}
				$last_pos = $pos+1;
    }

  }

  close FRAG or die $!;

}

sub read_fasta_file {
		my ($file) = @_;
	
		open (FASTA,$file) or die "$file: $!";

		$_ = <FASTA>;

		my $seq;
		while (<FASTA>) {
				chomp;
				$seq .= $_;
		}

		close FASTA;

		return $seq;
}

