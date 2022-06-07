#!/bin/bash

for munum in {0..5..1}
do
    export MUNUM=$munum
    sbatch job_WD.slurm
done