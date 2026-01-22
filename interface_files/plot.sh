#!/bin/bash

#SBATCH --job-name=PMMHD_kelvin_helmholz
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=krasny0

module load python/3.10.4 numpy/1.22.3 matplotlib/3.5.3 ffmpeg/6.0.0
python plot_particles.py -ps --time 0
# python plot_particles.py -ps --time 1
# python plot_particles.py -ps --time 2
# python plot_particles.py -ps --time 3
# python plot_particles.py -ps --time 4
python plot_particles.py -ps --time 5
# python plot_particles.py -ps --time 6
# python plot_particles.py -ps --time 7
# python plot_particles.py -ps --time 8
# python plot_particles.py -ps --time 9
python plot_particles.py -ps --time 10
# python plot_particles.py -ps --time 11
# python plot_particles.py -ps --time 12
# python plot_particles.py -ps --time 13
# python plot_particles.py -ps --time 14
python plot_particles.py -ps --time 15
# python plot_particles.py -ps --time 16
# python plot_particles.py -ps --time 17
# python plot_particles.py -ps --time 18
# python plot_particles.py -ps --time 19
python plot_particles.py -ps --time 20