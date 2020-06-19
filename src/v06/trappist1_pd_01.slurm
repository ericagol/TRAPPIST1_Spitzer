#!/bin/bash
## Run computation of run Photodynamic model for TRAPPIST-1 planets with limbprior, contamination, incprior toggled on and off (8 options x 3 runs)
#SBATCH --job-name=T1_pd_run
## Allocation Definition
## On mox and ikt, the account and partition options should be the same.
#SBATCH --account=astro
#SBATCH --partition=astro
## Resources
## Nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
## Walltime (2 days). Do not specify a walltime substantially more than your job needs.
#SBATCH --time=4-00:00:00
## Memory per node. It is important to specify the memory since the default memory is very small.
## For mox, --mem may be more than 300G depending on the memory of your nodes.
## For ikt, --mem may be 58G or more depending on the memory of your nodes.
## See above section on "Specifying memory" for choices for --mem.
#SBATCH --mem=50G
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/astro/agol/
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=agol@uw.edu
## export all your environment variables to the batch job session
#SBATCH --export=ALL


export JULIA_DEPOT_PATH="/gscratch/astro/agol/julia-1.0.5/bin/"
export PATH=/gscratch/astro/agol/julia-1.0.5/bin/:$PATH
export JULIA_PKGDIR="/usr/lusers/agol/.julia"

julia run_pd_background_01.jl 001 &
julia run_pd_background_01.jl 002 &
julia run_pd_background_01.jl 003 &

julia run_pd_background_02.jl 001 &
julia run_pd_background_02.jl 002 &
julia run_pd_background_02.jl 003 &

julia run_pd_background_03.jl 001 &
julia run_pd_background_03.jl 002 &
julia run_pd_background_03.jl 003 &

julia run_pd_background_04.jl 001 &
julia run_pd_background_04.jl 002 &
julia run_pd_background_04.jl 003 &

julia run_pd_background_05.jl 001 &
julia run_pd_background_05.jl 002 &
julia run_pd_background_05.jl 003 &

julia run_pd_background_06.jl 001 &
julia run_pd_background_06.jl 002 &
julia run_pd_background_06.jl 003 &

julia run_pd_background_07.jl 001 &
julia run_pd_background_07.jl 002 &
julia run_pd_background_07.jl 003 &

julia run_pd_background_08.jl 001 &
julia run_pd_background_08.jl 002 &
julia run_pd_background_08.jl 003 &

wait
