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
if VERSION <= v"1.0.1"
  Pkg.add(name="Optim",version="0.22.0")
else
  Pkg.add("Optim")
end
Pkg.add("Documenter")
Pkg.add("GSL")
Pkg.add("DelimitedFiles")
Pkg.add("IterativeSolvers")
Pkg.add("JLD2")
Pkg.add("MCMCDiagnostics")
Pkg.clone("https://github.com/rodluger/Limbdark.jl.git")
