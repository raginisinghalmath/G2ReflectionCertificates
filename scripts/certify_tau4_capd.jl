#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","G2CapdCertificates.jl"))
using .G2CapdCertificates

try
    certify_capd(:tau4;
                 start_order=24,
                 start_delta=BigInt(1)//BigInt(8),
                 radius_scale=BigInt(1)//BigInt(10_000),
                 capd_order=24,
                 centre_tolerance="1e-32",
                 box_tolerance="1e-24")
catch error
    println(stderr,"NOT CERTIFIED: ",sprint(showerror,error))
    rethrow()
end
