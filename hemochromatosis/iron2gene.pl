#!/usr/bin/perl
my @file=glob("*iron.pathway.txt");
foreach my $file(@file){
open F,$file;
while(<F>){
next if /^\s+/;
my($id,$gene,$tmp)=split/[\s+|;]/,$_;
print "$gene\n";
}
}
