FROM leethargo/scip-julia

# needed to compile GLPK from source
RUN apt-get install --no-install-recommends -qq libgmp-dev libglpk-dev

# prepare julia dependencies (fixed version 1.0; not SCIP!)
RUN /test/julia-1.0/bin/julia -e "using Pkg; Pkg.add(\"GLPKMathProgInterface\"); Pkg.add(\"SCS\")"
