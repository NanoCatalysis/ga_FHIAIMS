#!/bin/bash
for i in 9 10 11 12 13 14 15
do
 echo "$i" > Au$i/Au$i.xyz
 grep "Total energy of the DFT / Hartree-Fock s.c.f. calculation" Au10/au10.out | head -1 | awk '{print "# eV" $12}'>>Au$i/Au$i.xyz	
 awk '{if(NR >= 6) print $5"   "$2"   "$3"   "$4 }' Au$i/geometry.in.next_step >>Au$i/Au$i.xyz
done
