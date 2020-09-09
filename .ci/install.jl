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
Pkg.add("ForwardDiff")  # 6/30/2020
Pkg.add("DiffResults")  #
Pkg.add(PackageSpec(name="Optim",version="0.22.0"))
Pkg.free(PackageSpec("Optim")
Pkg.pin(PackageSpec(name="Optim",version="0.22.0"))
Pkg.add("Documenter")
Pkg.add("GSL")
Pkg.add("DelimitedFiles")
Pkg.add("IterativeSolvers") # 4/6/2020
Pkg.add("JLD2")
Pkg.add("MCMCDiagnostics")
if VERSION <= v"1.0.1"
  Pkg.clone("https://github.com/rodluger/Limbdark.jl.git")
else
  Pkg.add(PackageSpec(url="https://github.com/rodluger/Limbdark.jl.git"))
end
