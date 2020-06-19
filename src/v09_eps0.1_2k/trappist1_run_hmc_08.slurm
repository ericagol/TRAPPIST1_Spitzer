#!/bin/bash
## Run computation of run HMC for TRAPPIST-1 planets
#SBATCH --job-name=trappist1_hmc_04
## Allocation Definition
## On mox and ikt, the account and partition options should be the same.
#SBATCH --account=astro
#SBATCH --partition=astro
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (10 days). Do not specify a walltime substantially more than your job needs.
#SBATCH --time=10-00:00:00
## Memory per node. It is important to specify the memory since the default memory is very small.
## For mox, --mem may be more than 300G depending on the memory of your nodes.
## For ikt, --mem may be 58G or more depending on the memory of your nodes.
## See above section on "Specifying memory" for choices for --mem.
#SBATCH --mem=120G
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

julia trappist1_run_hmc.jl 291 &> trappist1_hmc_291.txt &
julia trappist1_run_hmc.jl 292 &> trappist1_hmc_292.txt &
julia trappist1_run_hmc.jl 293 &> trappist1_hmc_293.txt &
julia trappist1_run_hmc.jl 294 &> trappist1_hmc_294.txt &
julia trappist1_run_hmc.jl 295 &> trappist1_hmc_295.txt &
julia trappist1_run_hmc.jl 296 &> trappist1_hmc_296.txt &
julia trappist1_run_hmc.jl 297 &> trappist1_hmc_297.txt &
julia trappist1_run_hmc.jl 298 &> trappist1_hmc_298.txt &
julia trappist1_run_hmc.jl 299 &> trappist1_hmc_299.txt &
julia trappist1_run_hmc.jl 300 &> trappist1_hmc_300.txt &

julia trappist1_run_hmc.jl 301 &> trappist1_hmc_301.txt &
julia trappist1_run_hmc.jl 302 &> trappist1_hmc_302.txt &
julia trappist1_run_hmc.jl 303 &> trappist1_hmc_303.txt &
julia trappist1_run_hmc.jl 304 &> trappist1_hmc_304.txt &
julia trappist1_run_hmc.jl 305 &> trappist1_hmc_305.txt &
julia trappist1_run_hmc.jl 306 &> trappist1_hmc_306.txt &
julia trappist1_run_hmc.jl 307 &> trappist1_hmc_307.txt &
julia trappist1_run_hmc.jl 308 &> trappist1_hmc_308.txt &
julia trappist1_run_hmc.jl 309 &> trappist1_hmc_309.txt &
julia trappist1_run_hmc.jl 310 &> trappist1_hmc_310.txt &

julia trappist1_run_hmc.jl 311 &> trappist1_hmc_311.txt &
julia trappist1_run_hmc.jl 312 &> trappist1_hmc_312.txt &
julia trappist1_run_hmc.jl 313 &> trappist1_hmc_313.txt &
julia trappist1_run_hmc.jl 314 &> trappist1_hmc_314.txt &
julia trappist1_run_hmc.jl 315 &> trappist1_hmc_315.txt &
julia trappist1_run_hmc.jl 316 &> trappist1_hmc_316.txt &
julia trappist1_run_hmc.jl 317 &> trappist1_hmc_317.txt &
julia trappist1_run_hmc.jl 318 &> trappist1_hmc_318.txt &

wait
