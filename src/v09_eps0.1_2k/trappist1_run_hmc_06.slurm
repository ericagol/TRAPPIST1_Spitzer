#!/bin/bash
## Run computation of run HMC for TRAPPIST-1 planets
#SBATCH --job-name=trappist1_hmc_02
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

julia trappist1_run_hmc.jl 231 &> trappist1_hmc_231.txt &
julia trappist1_run_hmc.jl 232 &> trappist1_hmc_232.txt &
julia trappist1_run_hmc.jl 233 &> trappist1_hmc_233.txt &
julia trappist1_run_hmc.jl 234 &> trappist1_hmc_234.txt &
julia trappist1_run_hmc.jl 235 &> trappist1_hmc_235.txt &
julia trappist1_run_hmc.jl 236 &> trappist1_hmc_236.txt &
julia trappist1_run_hmc.jl 237 &> trappist1_hmc_237.txt &
julia trappist1_run_hmc.jl 238 &> trappist1_hmc_238.txt &
julia trappist1_run_hmc.jl 239 &> trappist1_hmc_239.txt &
julia trappist1_run_hmc.jl 240 &> trappist1_hmc_240.txt &

julia trappist1_run_hmc.jl 241 &> trappist1_hmc_241.txt &
julia trappist1_run_hmc.jl 242 &> trappist1_hmc_242.txt &
julia trappist1_run_hmc.jl 243 &> trappist1_hmc_243.txt &
julia trappist1_run_hmc.jl 244 &> trappist1_hmc_244.txt &
julia trappist1_run_hmc.jl 245 &> trappist1_hmc_245.txt &
julia trappist1_run_hmc.jl 246 &> trappist1_hmc_246.txt &
julia trappist1_run_hmc.jl 247 &> trappist1_hmc_247.txt &
julia trappist1_run_hmc.jl 248 &> trappist1_hmc_248.txt &
julia trappist1_run_hmc.jl 249 &> trappist1_hmc_249.txt &
julia trappist1_run_hmc.jl 250 &> trappist1_hmc_250.txt &

julia trappist1_run_hmc.jl 251 &> trappist1_hmc_251.txt &
julia trappist1_run_hmc.jl 252 &> trappist1_hmc_252.txt &
julia trappist1_run_hmc.jl 253 &> trappist1_hmc_253.txt &
julia trappist1_run_hmc.jl 254 &> trappist1_hmc_254.txt &
julia trappist1_run_hmc.jl 255 &> trappist1_hmc_255.txt &
julia trappist1_run_hmc.jl 256 &> trappist1_hmc_256.txt &
julia trappist1_run_hmc.jl 257 &> trappist1_hmc_257.txt &
julia trappist1_run_hmc.jl 258 &> trappist1_hmc_258.txt &

wait
