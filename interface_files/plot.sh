#!/bin/bash

#SBATCH --job-name=PMMHD_kelvin_helmholz
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --account=krasny0

module load python/3.10.4 numpy/1.22.3 matplotlib/3.5.3 ffmpeg/6.0.0
# python plot_particles.py -ps --time 0
# python plot_particles.py -ps --time 1
for t in {0..20}; do
  python plot_particles.py -ps --time "$t"
done