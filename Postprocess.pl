#!/usr/bin/perl

use strict;
use warnings;

open(my $in, '<', 'output.csv')
    or die "Cannot open input.txt: $!";

open(my $out, '>', 'EnergyVals.csv')
    or die "Cannot open output.txt: $!";

while (<$in>) {
      if (index($_, "EnergyOutput ") != -1) {
      my $line = $_;
      $line =~ s/EnergyOutput //;
      print $out $line;
    }
  }

  close($in);
  close($out);
