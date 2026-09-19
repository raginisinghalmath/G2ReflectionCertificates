#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","G2CapdCertificates.jl"))
using .G2CapdCertificates

try
    build_capd_driver(force=true)
    println("CAPD SETUP PASS")
catch error
    println(stderr,"CAPD SETUP FAILED: ",sprint(showerror,error))
    rethrow()
end
