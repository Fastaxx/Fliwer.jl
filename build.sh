#!/bin/bash
set -x
julia --project --color=yes \
       --check-bounds=yes \
    -e 'using Pkg;
       Pkg.develop(PackageSpec(path="/Users/Louis/Desktop/Doc/CartesianGeometry.jl")); 
       Pkg.instantiate()'
julia --project --color=yes --check-bounds=yes \
       -e 'using Pkg; Pkg.instantiate(); Pkg.build();'
julia test/runtests.jl