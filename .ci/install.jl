ENV["PYTHON"]=""
using Pkg
Pkg.add("Conda")
using Conda
Conda.update()
Conda.add("matplotlib")
Conda.add_channel("conda-forge")
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("PyPlot")
Pkg.add("SpecialFunctions")
Pkg.add("ForwardDiff")
Pkg.add("DiffResults")
Pkg.add("Optim")
Pkg.add("Documenter")
Pkg.add("GSL")
Pkg.add("DelimitedFiles")
Pkg.add("IterativeSolvers")
Pkg.add("JLD2")