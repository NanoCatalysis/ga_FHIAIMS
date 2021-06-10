#!/bin/bash
for i in 10 11 12 13 14 15
do
 cd Au$i && bsub < qsub_fhi.sh && cd ..
done
