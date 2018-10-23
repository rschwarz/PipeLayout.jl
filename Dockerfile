FROM leethargo/scip-julia

# needed to compile GLPK from source
RUN apt-get install --no-install-recommends -qq libgmp-dev libglpk-dev

# prepare julia dependencies (fixed version 0.7; not SCIP!)
RUN /test/julia-0.7/bin/julia -e "using Pkg; Pkg.add(\"GLPKMathProgInterface\"); Pkg.add(\"SCS\")"
