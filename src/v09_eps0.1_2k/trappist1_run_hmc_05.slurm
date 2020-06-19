#!/bin/bash
## Run computation of run HMC for TRAPPIST-1 planets
#SBATCH --job-name=trappist1_hmc_01
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

julia trappist1_run_hmc.jl 201 &> trappist1_hmc_201.txt &
julia trappist1_run_hmc.jl 202 &> trappist1_hmc_202.txt &
julia trappist1_run_hmc.jl 203 &> trappist1_hmc_203.txt &
julia trappist1_run_hmc.jl 204 &> trappist1_hmc_204.txt &
julia trappist1_run_hmc.jl 205 &> trappist1_hmc_205.txt &
julia trappist1_run_hmc.jl 206 &> trappist1_hmc_206.txt &
julia trappist1_run_hmc.jl 207 &> trappist1_hmc_207.txt &
julia trappist1_run_hmc.jl 208 &> trappist1_hmc_208.txt &
julia trappist1_run_hmc.jl 209 &> trappist1_hmc_209.txt &
julia trappist1_run_hmc.jl 210 &> trappist1_hmc_210.txt &

julia trappist1_run_hmc.jl 211 &> trappist1_hmc_211.txt &
julia trappist1_run_hmc.jl 212 &> trappist1_hmc_212.txt &
julia trappist1_run_hmc.jl 213 &> trappist1_hmc_213.txt &
julia trappist1_run_hmc.jl 214 &> trappist1_hmc_214.txt &
julia trappist1_run_hmc.jl 215 &> trappist1_hmc_215.txt &
julia trappist1_run_hmc.jl 216 &> trappist1_hmc_216.txt &
julia trappist1_run_hmc.jl 217 &> trappist1_hmc_217.txt &
julia trappist1_run_hmc.jl 218 &> trappist1_hmc_218.txt &
julia trappist1_run_hmc.jl 219 &> trappist1_hmc_219.txt &
julia trappist1_run_hmc.jl 220 &> trappist1_hmc_220.txt &

julia trappist1_run_hmc.jl 221 &> trappist1_hmc_221.txt &
julia trappist1_run_hmc.jl 222 &> trappist1_hmc_222.txt &
julia trappist1_run_hmc.jl 223 &> trappist1_hmc_223.txt &
julia trappist1_run_hmc.jl 224 &> trappist1_hmc_224.txt &
julia trappist1_run_hmc.jl 225 &> trappist1_hmc_225.txt &
julia trappist1_run_hmc.jl 226 &> trappist1_hmc_226.txt &
julia trappist1_run_hmc.jl 227 &> trappist1_hmc_227.txt &
julia trappist1_run_hmc.jl 228 &> trappist1_hmc_228.txt &

wait
