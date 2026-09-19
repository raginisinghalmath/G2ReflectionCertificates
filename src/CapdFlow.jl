module CapdFlow

using IntervalArithmetic
using LinearAlgebra
using SHA
using ..ValidatedArithmetic
using ..G2Equations

export build_capd_driver, capd_field, capd_flow

const DRIVER_SOURCE = normpath(joinpath(@__DIR__, "..", "deps",
                                         "g2_capd_flow.cpp"))
const BUILD_DIRECTORY = normpath(joinpath(@__DIR__, "..", "deps", "build"))
const DRIVER_BINARY = joinpath(BUILD_DIRECTORY,
                               Sys.iswindows() ? "g2_capd_flow.exe" :
                                                 "g2_capd_flow")
const DRIVER_STAMP = DRIVER_BINARY * ".sha256"

function rhs_strings(q)
    q1,q2,q3,q4,q5,q6,q7,q8=q
    d1="($q2+$q7)"
    d2="($q3+$q6)"
    d3="(-($q4+$q5))"
    p1="(-sqrt(($d1*$d2)/(lam*$d3)))"
    p2="sqrt(($d1*$d3)/(lam*$d2))"
    p3="sqrt(($d2*$d3)/(lam*$d1))"
    # On the certified chamber d1,d2,d3>0, the original coefficient
    # lam/(2*p1*p2*p3) equals this expression, with its negative sign fixed.
    ell="(-lam/(2*sqrt(($d1*$d2*$d3)/(lam*lam*lam))))"
    a1="($q1*$q8)"; a2="($q2*$q7)"
    a3="($q3*$q6)"; a4="($q4*$q5)"
    total="($a1+$a2+$a3+$a4)"
    [
      "$ell*($q1*(2*$a1-$total)+2*$q2*$q3*$q5)",
      "$p3-$ell*($q2*(2*$a2-$total)+2*$q1*$q4*$q6)",
      "$p2-$ell*($q3*(2*$a3-$total)+2*$q1*$q4*$q7)",
      "-$p1+$ell*($q4*(2*$a4-$total)+2*$q2*$q3*$q8)",
      "$p1-$ell*($q5*(2*$a4-$total)+2*$q1*$q6*$q7)",
      "-$p2+$ell*($q6*(2*$a3-$total)+2*$q2*$q5*$q8)",
      "-$p3+$ell*($q7*(2*$a2-$total)+2*$q3*$q5*$q8)",
      "-$ell*($q8*(2*$a1-$total)+2*$q4*$q6*$q7)"
    ]
end

"""CAPD parser expression for `(q,T)'=((T-delta)F(q),0)`."""
function capd_field(kind)
    if kind==:tau4
        names=["q$i" for i in 1:8]
        logical=copy(names)
        selected=1:8
    elseif kind==:tau23_u1
        names=["q1","q2","q4","q5","q6","q8"]
        logical=["q1","q2","q2","q4","q5","q6","q6","q8"]
        selected=[1,2,4,5,6,8]
    else
        error("Unknown case")
    end
    rhs=rhs_strings(logical)
    scaled=["(tt-delta)*($(rhs[i]))" for i in selected]
    "time:s;par:lam,delta;var:"*join(vcat(names,["tt"]),",")*
        ";fun:"*join(vcat(scaled,["0"]),",")*";"
end

function flag_output(command)
    String.(split(strip(read(command,String))))
end

function build_attempt(command)
    mktemp() do _,errors
        process=run(pipeline(ignorestatus(command),stdout=devnull,stderr=errors))
        flush(errors); seekstart(errors)
        success(process),read(errors,String)
    end
end

"""Compile the audited multiprecision CAPD driver, following Alberto's setup."""
function build_capd_driver(; force=false)
    isfile(DRIVER_SOURCE) || error("Missing CAPD driver source: $DRIVER_SOURCE")
    digest=bytes2hex(sha256(read(DRIVER_SOURCE)))
    if !force && isfile(DRIVER_BINARY) && isfile(DRIVER_STAMP) &&
       strip(read(DRIVER_STAMP,String))==digest
        return DRIVER_BINARY
    end
    mkpath(BUILD_DIRECTORY)
    compiler=get(ENV,"CXX",Sys.isapple() ? "/usr/bin/clang++" : "g++")
    prefix=get(ENV,"CAPD_PREFIX",joinpath(homedir(),"capd"))
    configs=String[]
    if haskey(ENV,"CAPD_CONFIG") && isfile(ENV["CAPD_CONFIG"])
        push!(configs,ENV["CAPD_CONFIG"])
    end
    local_config=joinpath(prefix,"bin","capd-config")
    isfile(local_config) && push!(configs,local_config)
    path_config=Sys.which("capd-config")
    path_config===nothing || push!(configs,path_config)
    unique!(configs)

    attempts=Tuple{Vector{String},Vector{String}}[]
    for config in configs
        try
            push!(attempts,(flag_output(`$config --cflags`),
                            flag_output(`$config --libs`)))
        catch
        end
    end
    fallback_cflags=["-I"*joinpath(prefix,"include")]
    if Sys.isapple() && Sys.ARCH in (:aarch64,:arm64)
        append!(fallback_cflags,["-D__USE_NATIVE__","-D__HAVE_MPFR__"])
    end
    push!(attempts,(fallback_cflags,["-L"*joinpath(prefix,"lib"),"-lcapd"]))

    gmp_prefix=get(ENV,"GMP_PREFIX",
                   Sys.isapple() && isdir("/opt/homebrew") ? "/opt/homebrew" : "")
    extra_compile=isempty(gmp_prefix) ? String[] : ["-I"*joinpath(gmp_prefix,"include")]
    extra_link=isempty(gmp_prefix) ? String[] : ["-L"*joinpath(gmp_prefix,"lib")]
    messages=String[]
    for (cflags,lflags) in attempts
        args=vcat([compiler,"-O2","-std=c++17",DRIVER_SOURCE],cflags,
                  extra_compile,["-o",DRIVER_BINARY],lflags,extra_link,
                  ["-lmpfr","-lgmp"])
        ok,err=build_attempt(Cmd(args))
        if ok
            write(DRIVER_STAMP,digest*"\n")
            println("CAPD driver compiled: ",DRIVER_BINARY)
            return DRIVER_BINARY
        end
        lines=filter(x->!isempty(x),split(strip(err),'\n'))
        push!(messages,isempty(lines) ? "unknown compiler error" : first(lines))
    end
    error("CAPD driver compilation failed:\n"*join(messages,"\n")*
          "\nSet CAPD_CONFIG or CAPD_PREFIX to the multiprecision CAPD installation.")
end

decimal(x::BigFloat)=string(x)
decimal(x::Real)=string(BigFloat(x))

function parse_interval(lower,upper)
    checked(iv(parse(BigFloat,lower),parse(BigFloat,upper)))
end

function parse_output(output,n)
    state=Vector{IV}(undef,n)
    derivative=Matrix{IV}(undef,n,n)
    state_seen=falses(n); derivative_seen=falses(n,n)
    margins=fill(BigFloat(NaN),3)
    steps=-1; ok=false
    for line in split(output,'\n')
        fields=split(strip(line)); isempty(fields) && continue
        if fields[1]=="STATE"
            i=parse(Int,fields[2])+1
            state[i]=parse_interval(fields[3],fields[4]); state_seen[i]=true
        elseif fields[1]=="MAT"
            i=parse(Int,fields[2])+1; j=parse(Int,fields[3])+1
            derivative[i,j]=parse_interval(fields[4],fields[5])
            derivative_seen[i,j]=true
        elseif fields[1]=="MARGIN"
            margins[parse(Int,fields[2])+1]=parse(BigFloat,fields[3])
        elseif fields[1]=="STEPS"
            steps=parse(Int,fields[2])
        elseif fields[1]=="OK"
            ok=true
        end
    end
    require(ok,"CAPD driver did not report success")
    require(all(state_seen) && all(derivative_seen),"Incomplete CAPD output")
    require(all(isfinite,margins) && all(x->x>0,margins),
            "CAPD did not report positive chamber margins")
    require(steps>0,"CAPD reported no validated steps")
    checked(state); checked(derivative)
    (state=state,derivative=derivative,margins=margins,steps=steps)
end

"""
Rigorous CAPD C1 integration of the rescaled endpoint problem
`dq/ds=(T-delta)F(q), T'=0`, `0<=s<=1`.
"""
function capd_flow(start,terminal_time,lam,kind;
                   order=24,tolerance="1e-26",bits=256)
    d,_=dimensions(kind); n=d+1
    T=checked(terminal_time)
    require(inf(T)>BigFloat(start.delta),"Terminal time is not after delta")
    z=checked(vcat(start.q,[T]))
    require(length(z)==n,"Wrong CAPD state dimension")
    binary=build_capd_driver()
    output=mktemp() do path,io
        write(io,capd_field(kind)); close(io)
        args=String[binary,kind==:tau4 ? "4" : "23",path,string(order),
                    string(tolerance),string(bits),decimal(inf(lam)),
                    decimal(sup(lam)),decimal(start.delta)]
        for x in z
            push!(args,decimal(inf(x)),decimal(sup(x)))
        end
        read(pipeline(Cmd(args),stderr=stderr),String)
    end
    result=parse_output(output,n)
    chamber(result.state[1:d],kind)
    println("CAPD PASS: steps=",result.steps,
            " chamber margins >= ",result.margins)
    (q=result.state[1:d],state=result.state,M=result.derivative,
     chamber_margins=result.margins,steps=result.steps)
end

end
