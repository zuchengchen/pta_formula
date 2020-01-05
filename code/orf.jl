using Distributed
cpuNum = length(Sys.cpu_info()) 
# number of cpu cores

addprocs(cpuNum - 2)
# using all cpu cores

@everywhere using HCubature
using DelimitedFiles:readdlm

@everywhere polars = ["TT", "ST", "VL", "SL"]

@everywhere function e_factor(xi, L1, L2, f, theta, phi)
    phase1 = 2im*f*L1*pi*(1.0 + cos(theta))
    phase2 = -2im*f*L2*pi*(1.0 + cos(theta)*cos(xi) + cos(phi)*sin(theta)*sin(xi))
    temp1 = -1.0 + exp(phase1)
    temp2 = -1.0 + exp(phase2)
    real(temp1*temp2)
end

@everywhere function Gamma(polar, xi, L1, L2, f; psrTerm=true)
    if polar == "TT"
        Gamma0 = GammaTT0
    elseif polar == "ST"
        Gamma0 = GammaST0
    elseif polar == "VL"
        Gamma0 = GammaVL0
    elseif polar == "SL"
        Gamma0 = GammaSL0
    end
    if psrTerm
        integral = x -> Gamma0(xi, x[1], x[2])*e_factor(xi, L1, L2, f, x[1], x[2])
    else
        integral = x -> Gamma0(xi, x[1], x[2])
    end
    hcubature(integral, [0.0, 0.0], [pi, 2pi], rtol=1e-4, maxevals=10^8, initdiv=2)[1]
end

function Gammas(polar, xi, L1, L2, fs; psrTerm=true)
    pmap(f -> Gamma(polar, xi, L1, L2, f, psrTerm=psrTerm), fs)
end

@everywhere function GammaTT0(xi, theta, phi)
    return 3.0/8.0/pi * (sin(theta/2.0)^2.0 * sin(theta) * (
        (cos(xi)*sin(theta) - cos(theta)*cos(phi)*sin(xi))^2.0 
        - sin(xi)^2.0 * sin(phi)^2.0)
     )/(1.0 + cos(theta)*cos(xi) + cos(phi)*sin(theta)*sin(xi))
end

function GammaTTA(xi, L1, L2, f)
    k = (1.0-cos(xi))/2.0
    return 3.0*(1.0/3.0 + k*(log(k)-1.0/6.0))
end

@everywhere function GammaST0(xi, theta, phi)
    return -3.0/8.0/pi * (sin(theta/2.0))^2.0 * sin(theta) * (
        -1 + cos(theta)*cos(xi) + cos(phi)*sin(theta)*sin(xi)
    )
end

# Analytical version of GammaST
function GammaSTA(xi, L1, L2, f)
    return 1/4 * (3 + cos(xi))
end

@everywhere function GammaVL0(xi, theta, phi)
    return -3.0/8.0/pi * cos(theta) * sin(theta) * (
        sin(2.0*theta)*(-cos(xi)^2.0 + cos(phi)^2.0 * sin(xi)^2.0)
        + cos(2.0*theta) * cos(phi) * sin(2.0*xi) 
     )/(1.0 + cos(theta))/(cos(xi)*cos(theta)/sin(theta) + 1/sin(theta) + cos(phi)*sin(xi))
end

@everywhere function GammaSL0(xi, theta, phi)
    return 3.0/16.0/pi * cos(theta)^2.0 * sin(theta) * (
        cos(theta)*cos(xi) + cos(phi)*sin(theta)*sin(xi)
    )^2.0 /(1.0 + cos(theta))/(1 + cos(theta)*cos(xi) + cos(phi)*sin(theta)*sin(xi))
end

M2Data = readdlm("backup/prs_pair_info.txt");

function saveOnePair(polar, xi, L1, L2)
    nmode = 50 
    fmin = 2.774455633345974e-09
    fmax = nmode*fmin
    fs = range(fmin, fmax, length=3nmode)
    file = string("backup/", polar, "_", xi, "_", L1, "_", L2, ".txt")
    
    yr = 365.25 * 24.0 * 3600.0
    kpc_to_sec = 3.26156e3 * yr
    L1 *= kpc_to_sec
    L2 *= kpc_to_sec
    Gs = Gammas(polar, xi, L1, L2, fs)
    open(file, "w") do io
        for i in 1:length(fs)
            write(io, string(fs[i]) * "  " * string(Gs[i]) * "\n")
        end
    end
end


function saveAllPair(polar)
    for i in 1:size(M2Data)[1]
        xi, L1, L2 = M2Data[i,:]
        @time saveOnePair(polar, xi, L1, L2)
    end
end

polarToCalculate = ["VL", "SL"]
for polar in polarToCalculate
    saveAllPair(polar)        
end