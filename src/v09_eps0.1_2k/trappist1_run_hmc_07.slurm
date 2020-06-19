#!/bin/bash
## Run computation of run HMC for TRAPPIST-1 planets
#SBATCH --job-name=trappist1_hmc_03
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
## For mox, --mem may be more than 100G depending on the memory of your nodes.
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

julia trappist1_run_hmc.jl 261 &> trappist1_hmc_261.txt &
julia trappist1_run_hmc.jl 262 &> trappist1_hmc_262.txt &
julia trappist1_run_hmc.jl 263 &> trappist1_hmc_263.txt &
julia trappist1_run_hmc.jl 264 &> trappist1_hmc_264.txt &
julia trappist1_run_hmc.jl 265 &> trappist1_hmc_265.txt &
julia trappist1_run_hmc.jl 266 &> trappist1_hmc_266.txt &
julia trappist1_run_hmc.jl 267 &> trappist1_hmc_267.txt &
julia trappist1_run_hmc.jl 268 &> trappist1_hmc_268.txt &
julia trappist1_run_hmc.jl 269 &> trappist1_hmc_269.txt &
julia trappist1_run_hmc.jl 270 &> trappist1_hmc_270.txt &

julia trappist1_run_hmc.jl 271 &> trappist1_hmc_271.txt &
julia trappist1_run_hmc.jl 272 &> trappist1_hmc_272.txt &
julia trappist1_run_hmc.jl 273 &> trappist1_hmc_273.txt &
julia trappist1_run_hmc.jl 274 &> trappist1_hmc_274.txt &
julia trappist1_run_hmc.jl 275 &> trappist1_hmc_275.txt &
julia trappist1_run_hmc.jl 276 &> trappist1_hmc_276.txt &
julia trappist1_run_hmc.jl 277 &> trappist1_hmc_277.txt &
julia trappist1_run_hmc.jl 278 &> trappist1_hmc_278.txt &
julia trappist1_run_hmc.jl 279 &> trappist1_hmc_279.txt &
julia trappist1_run_hmc.jl 280 &> trappist1_hmc_280.txt &

julia trappist1_run_hmc.jl 281 &> trappist1_hmc_281.txt &
julia trappist1_run_hmc.jl 282 &> trappist1_hmc_282.txt &
julia trappist1_run_hmc.jl 283 &> trappist1_hmc_283.txt &
julia trappist1_run_hmc.jl 284 &> trappist1_hmc_284.txt &
julia trappist1_run_hmc.jl 285 &> trappist1_hmc_285.txt &
julia trappist1_run_hmc.jl 286 &> trappist1_hmc_286.txt &
julia trappist1_run_hmc.jl 287 &> trappist1_hmc_287.txt &
julia trappist1_run_hmc.jl 288 &> trappist1_hmc_288.txt &

wait
