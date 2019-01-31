#!/bin/bash

mkdir -p slurmoutput

jout1=$(sbatch runfirst.h $1 $2)

jid1="${jout1//[!0-9]/}"

echo $jid1

jout2=$(sbatch --dependency=afterok:$jid1 runsecond.h $1 $2)
