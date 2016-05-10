#!/usr/bin/perl

use strict;
use warnings;

open(my $in, '<', 'output.csv')
    or die "Cannot open input.txt: $!";

open(my $out, '>', 'EnergyVals.csv')
    or die "Cannot open output.txt: $!";

while (<$in>) {
      print $out $_ unless /CLAW/;
  }

  close($in);
  close($out);
