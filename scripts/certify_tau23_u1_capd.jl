#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","G2CapdCertificates.jl"))
using .G2CapdCertificates

try
    certify_capd(:tau23_u1;
                 start_order=8,
                 start_delta=BigInt(1)//BigInt(1024),
                 radius_scale=1,
                 capd_order=24,
                 centre_tolerance="1e-32",
                 box_tolerance="1e-24")
catch error
    println(stderr,"NOT CERTIFIED: ",sprint(showerror,error))
    rethrow()
end
