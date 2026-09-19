module G2Equations

using IntervalArithmetic
using ..ValidatedArithmetic

export vectorfield, regularized, initial_pair, q_from_y, lift, select
export dimensions, matching, chamber, parameters, lambda_value

lambda_value()=iv(6)/sqrt(iv(5))
dimensions(kind) = kind==:tau4 ? (8,3) : kind==:tau23_u1 ? (6,2) : error("Unknown case")
lift(y,kind) = kind==:tau4 ? y : [y[1],y[2],y[2],y[3],y[4],y[5],y[5],y[6]]
select(q,kind) = kind==:tau4 ? q : q[[1,2,4,5,6,8]]
function parameters(kind; center=false, radius_scale=1)
    cs,rs = if kind==:tau4
        (["0.042329982765393195","-0.06469940434457604",
          "0.006759263945977987","0.7656914313714904"],fill("1e-6",4))
    elseif kind==:tau23_u1
        (["0.025880416406250467","-0.07764124921873236","0.8578017401911906"],
         ["1e-7","3e-7","5e-7"])
    else
        error("Unknown case")
    end
    scale=checked(iv(radius_scale))
    require(inf(scale)>0,"Parameter radius scale must be positive")
    c=checked(iv.(cs)); r=checked(iv.(rs).*scale)
    (center=c,radii=r,box=center ? c : c+r.*iv(-1,1))
end

function p_from_q(q,lam)
    d1=q[2]+q[7]; d2=q[3]+q[6]; d3=q[4]+q[5]
    (-sqrt(-d1*d2/(lam*d3)),sqrt(-d1*d3/(lam*d2)),sqrt(-d2*d3/(lam*d1)))
end
function vectorfield(u,lam,kind)
    q1,q2,q3,q4,q5,q6,q7,q8=lift(u,kind)
    p1,p2,p3=p_from_q(lift(u,kind),lam)
    a1=q1*q8; a2=q2*q7; a3=q3*q6; a4=q4*q5
    s=a1+a2+a3+a4; d=iv(2)*p1*p2*p3; l=lam/d
    # Cancel removable q_i denominators algebraically BEFORE interval evaluation.
    f=[l*(q1*(iv(2)*a1-s)+iv(2)*q2*q3*q5),
       p3-l*(q2*(iv(2)*a2-s)+iv(2)*q1*q4*q6),
       p2-l*(q3*(iv(2)*a3-s)+iv(2)*q1*q4*q7),
       -p1+l*(q4*(iv(2)*a4-s)+iv(2)*q2*q3*q8),
       p1-l*(q5*(iv(2)*a4-s)+iv(2)*q1*q6*q7),
       -p2+l*(q6*(iv(2)*a3-s)+iv(2)*q2*q5*q8),
       -p3+l*(q7*(iv(2)*a2-s)+iv(2)*q3*q5*q8),
       -l*(q8*(iv(2)*a1-s)+iv(2)*q4*q6*q7)]
    select(f,kind)
end

function q_from_y(u,a,c,t,kind)
    y=lift(u,kind); t2=t*t
    select([a+t2*y[1],t*y[2],t*y[3],c+t2*y[4],
           -iv(3)*a+t2*y[5],t*y[6],t*y[7],-iv(3)*c+t2*y[8]],kind)
end

"""Analytic G in t*y'=G(t,y), with ALL t=0 cancellations performed."""
function regularized(u,a,c,t,lam,kind)
    y1,y2,y3,y4,y5,y6,y7,y8=lift(u,kind); t2=t*t
    q1=a+t2*y1; q4=c+t2*y4; q5=-iv(3)*a+t2*y5; q8=-iv(3)*c+t2*y8
    u1=y2+y7; u2=y3+y6; u3=c-iv(3)*a+t2*(y4+y5)
    P1=sqrt(-u1*u2/(lam*u3)); P2=sqrt(-u1*u3/(lam*u2)); P3=sqrt(-u2*u3/(lam*u1))
    l=lam/(-iv(2)*P1*P2*P3)
    v=a*y8-iv(3)*c*y1-c*y5+iv(3)*a*y4+t2*(y1*y8-y4*y5)
    s2=y2*y7; s3=y3*y6; w=s2+s3; h=-q1*q8-q4*q5
    g=[l*(q1*(v-w)+iv(2)*y2*y3*q5)-iv(2)*y1,
       P3-l*(y2*(h+t2*(s2-s3))+iv(2)*q1*q4*y6)-y2,
       P2-l*(y3*(h+t2*(s3-s2))+iv(2)*q1*q4*y7)-y3,
       P1+l*(q4*(-v-w)+iv(2)*y2*y3*q8)-iv(2)*y4,
       -P1-l*(q5*(-v-w)+iv(2)*q1*y6*y7)-iv(2)*y5,
       -P2+l*(y6*(h+t2*(s3-s2))+iv(2)*y2*q5*q8)-y6,
       -P3+l*(y7*(h+t2*(s2-s3))+iv(2)*y3*q5*q8)-y7,
       -l*(q8*(v-w)+iv(2)*q4*y6*y7)-iv(2)*y8]
    select(g,kind)
end

"""Exact y_0,y_1 formulas from the (3,-1), mu=-1 local existence lemma."""
function initial_pair(p,lam,kind)
    a,c=p[1:2]; nu=kind==:tau4 ? p[3] : zero(a)
    A=sqrt(iv(3)*a-c); B=sqrt(-a*c); C=sqrt(-lam*a*c*(iv(3)*a-c))
    D=sqrt(-lam*a*c/(iv(3)*a-c)); den=iv(8)*(iv(3)*a-c)
    b1=-iv(3)*lam^2*a*(iv(5)*a-c)/den-C/(iv(2)*c)
    b2=-lam*B+A/sqrt(lam)
    b4=-lam^2*c*(iv(9)*a-iv(5)*c)/den-C/(iv(2)*a)+D
    b5=iv(9)*lam^2*a*(iv(5)*a-c)/den+C/(iv(2)*c)-D
    b6=iv(3)*lam*B-A/sqrt(lam)
    b8=iv(3)*lam^2*c*(iv(9)*a-iv(5)*c)/den+C/(iv(2)*a)
    r=(A+iv(3)*lam*sqrt(lam)*B)/(A-lam*sqrt(lam)*B)
    z=zero(a)
    select([b1,b2,b2,b4,b5,b6,b6,b8],kind),select([z,-nu,nu,z,z,-r*nu,r*nu,z],kind)
end
function matching(q,kind)
    kind==:tau4 ? [q[1]-q[8],q[2]+q[4],q[3]-q[6],q[5]+q[7]] :
                 [q[1]-q[6],q[2]-q[5],q[3]-q[4]]
end
function chamber(u,kind)
    q=lift(u,kind); checked(q)
    ds=[q[2]+q[7],q[3]+q[6],-(q[4]+q[5])]
    for i in eachindex(ds)
        require(inf(ds[i])>0,
                "Real G2 chamber component $i not certified on this tube")
    end
    inf.(ds)
end

end
