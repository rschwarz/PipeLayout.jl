FROM leethargo/scip-julia

# needed to compile GLPK from source
RUN apt-get install --no-install-recommends -qq libgmp-dev libglpk-dev

# prepare julia dependencies (fixed version 0.6; not SCIP!)
RUN /test/julia-0.6/bin/julia -e "Pkg.add(\"Colors\"); Pkg.add(\"Combinatorics\"); Pkg.add(\"DataStructures\"); Pkg.add(\"GLPKMathProgInterface\"); Pkg.add(\"JSON\"); Pkg.add(\"JuMP\"); Pkg.add(\"LightGraphs\"); Pkg.add(\"RecipesBase\"); Pkg.add(\"SCS\"); Pkg.add(\"StaticArrays\"); Pkg.add(\"Unitful\")"
