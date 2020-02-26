#!/usr/bin/perl
#use strict;
my @file=glob("23E3*-F_*.seq");
foreach my $file(@file){
open F,$file;
print ">$file\n";
while(<F>){
chomp;
print $_;
}
print "\n";
}

#
#
