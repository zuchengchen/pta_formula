include("constant.jl")

"""
    designMatrixRed(
        toas;
        nMode=30, 
        TSpan=Nothing, 
        logFreq=false
    )

Construct an ``N_{toas} × 2N_{mode}`` Fourier design matrix from Eq.(11) of Lentati et al, 2013 
(https://arxiv.org/pdf/1210.3578.pdf)

# Arguments
- `toas::Array`: vector of time of arrial series in seconds.
- `nMode::Integer=30`: number of fourier coefficients to use.
- `TSpan:Float64=Nothing`: option to some other Tspan.
- `logFreq::Bool=false`: use log frequency spacing.

# Examples
```jldoctest
julia> toas = [1, 2]
julia> F, FFreqs = designMatrixRed(toas, nMode=2);
julia> F
2×4 Array{Float64,2}:
 -2.44929e-16  1.0  -4.89859e-16  1.0
 -4.89859e-16  1.0  -9.79717e-16  1.0
julia> FFreqs
4-element Array{Float64,1}:
 1.0
 1.0
 2.0
 2.0
```

# See also
[`designMatrixEnvelope`](@ref), [`designMatrixChromatic`](@ref), [`designMatrixDM`](@ref), [`designMatrixEphemeris`](@ref)
"""
function designMatrixRed(
        toas;
        nMode=30, 
        TSpan=Nothing, 
        logFreq=false
    )
    
    # get the span time
    T = (TSpan==Nothing) ? (maximum(toas) - minimum(toas)) : TSpan
    
    
    
    fMin = 1/T # minimum (starting) frequency 
    fMax = nMode/T # maximum (ending) frequency 
     
    # define sampling frequencies
    # frequency modes
    if logFreq
        fModes = 10 .^ range(log10(fMin), log10(fMax), length=nMode)
    else
        fModes = range(fMin, fMax, length=nMode)
    end
    
    # frequency modes of design matrix F
    FFreqs = repeat(fModes, inner=2)
    
    # F is a N_toa × 2N_mode design matrix
    F = zeros(length(toas), 2nMode)
    
    # set values by column
    for i in 1:nMode        
        F[:, 2i-1] .= sin.(2π * toas * fModes[i])
        F[:, 2i] .= cos.(2π * toas * fModes[i])
    end
    
    return F, FFreqs
end


"""
    designMatrixEnvelope(
        toas;
        log10Amp=-7.0,
        log10Q=log10(300.0),
        t0=53000.0day,
        nMode=30, 
        TSpan=Nothing, 
        logFreq=false
    )

Construct an ``N_{toas} × 2N_{mode}`` Fourier design matrix with gaussian envelope. 
Current normalization expresses DM signal as a deviation [seconds] at fref [MHz]

# Arguments
- `toas::Array`: vector of time of arrial series in seconds.
- `log10Amp::Float64=-7.0`: log10 of the amplitude [s].
- `log10Q::Float64=log10(300.0)`: log10 of standard deviation of gaussian envelope [days].
- `t0::Float64=53000.0day`: mean of gaussian envelope [s].
- `nMode::Integer=30`: number of fourier coefficients to use.
- `TSpan:Float64=Nothing`: option to some other Tspan.
- `logFreq::Bool=false`: use log frequency spacing.

# Examples
```jldoctest
julia> toas = [1, 1.015]*53000.0day
julia> F, FFreqs = designMatrixEnvelope(toas, nMode=2);
julia> F
2×4 Array{Float64,2}:
 -8.66025e-8  -5.0e-8      8.66025e-8  -5.0e-8    
 -2.58591e-9  -1.49298e-9  2.58591e-9  -1.49298e-9
julia> FFreqs
4-element Array{Float64,1}:
 1.455858374097388e-8
 1.455858374097388e-8
 2.911716748194776e-8
 2.911716748194776e-8
```

# See also
[`designMatrixRed`](@ref), [`designMatrixChromatic`](@ref), [`designMatrixDM`](@ref), [`designMatrixEphemeris`](@ref)
"""
function designMatrixEnvelope(
        toas;
        log10Amp=-7.0,
        log10Q=log10(300.0),
        t0=53000.0day,
        nMode=30, 
        TSpan=Nothing, 
        logFreq=false
    )
    
    # get base fourier design matrix and frequencies
    F, FFreqs = designMatrixRed(toas, nMode=nMode, TSpan=TSpan, logFreq=logFreq)  

    # compute gaussian envelope
    A = 10.0^log10Amp
    Q = 10.0^log10Q * day
    env = A * exp.(-(toas .- t0) .^ 2.0 / 2.0 / Q^2.0)
    
    for i in 1:size(F)[2]
        F[:, i] .= F[:, i] .* env
    end
    
    return F, FFreqs
end


"""
    designMatrixChromatic(
        toas,
        freqs;
        nMode=30, 
        TSpan=Nothing, 
        fRef=1400.0,
        index=4.0,
        logFreq=false
    )

Construct an ``N_{toas} × 2N_{mode}`` scattering-variation Fourier design matrix. 
Current normalization expresses scattering signal as a deviation [seconds] at fRef [MHz]

# Arguments
- `toas::Array`: vector of time of arrial series in seconds.
- `freqs:Array`: radio frequencies of observations [MHz].
- `nMode::Integer=30`: number of fourier coefficients to use.
- `TSpan:Float64=Nothing`: option to some other Tspan.
- `fRef:Float64=1400.0`: reference frequency [MHz].
- `index::Float64=4.0`: index of chromatic effects.
- `logFreq::Bool=false`: use log frequency spacing.

# Examples
```jldoctest
julia> toas = [1, 2]
julia> freqs = [3, 4]
julia> F, FFreqs = designMatrixChromatic(toas, freqs, nMode=2);
julia> F
2×4 Array{Float64,2}:
 -1.16163e-5  4.74272e10  -2.32326e-5  4.74272e10
 -7.35094e-6  1.50063e10  -1.47019e-5  1.50063e10
julia> FFreqs
4-element Array{Float64,1}:
 1.0
 1.0
 2.0
 2.0
```

# See also
[`designMatrixRed`](@ref), [`designMatrixEnvelope`](@ref), [`designMatrixDM`](@ref), [`designMatrixEphemeris`](@ref)
"""
function designMatrixChromatic(
        toas,
        freqs;
        nMode=30, 
        TSpan=Nothing, 
        fRef=1400.0,
        index=4.0,
        logFreq=false
    )
    
    # get base fourier design matrix and frequencies
    F, FFreqs = designMatrixRed(toas, nMode=nMode, TSpan=TSpan, logFreq=logFreq)

    # compute the DM-variation vectors
    Dm = (fRef ./ freqs) .^ index
    
    for i in 1:size(F)[2]
        F[:, i] .= F[:, i] .* Dm
    end

    return F, FFreqs
end


"""
    designMatrixDM(
        toas,
        freqs;
        nMode=30, 
        TSpan=Nothing, 
        fRef=1400.0,
        logFreq=false
    )

Construct an ``N_{toas} × 2N_{mode}`` dispersion measure (DM)-variation Fourier design matrix. 
Current normalization expresses DM signal as a deviation [seconds] at fRef [MHz].
It is a special case of designMatrixChromatic with index=2.

# Arguments
- `toas::Array`: vector of time of arrial series in seconds.
- `freqs:Array`: radio frequencies of observations [MHz].
- `nMode::Integer=30`: number of fourier coefficients to use.
- `TSpan:Float64=Nothing`: option to some other Tspan.
- `fRef:Float64=1400.0`: reference frequency [MHz].
- `logFreq::Bool=false`: use log frequency spacing.

# Examples
```jldoctest
julia> toas = [1, 2]
julia> freqs = [3, 4]
julia> F, FFreqs = designMatrixDM(toas, freqs, nMode=2)
julia> F
2×4 Array{Float64,2}:
 -5.33402e-11       2.17778e5  -1.0668e-10        2.17778e5
 -6.00077e-11  122500.0        -1.20015e-10  122500.0 
julia> FFreqs
4-element Array{Float64,1}:
 1.0
 1.0
 2.0
 2.0
```

# See also
[`designMatrixRed`](@ref), [`designMatrixEnvelope`](@ref), [`designMatrixChromatic`](@ref), [`designMatrixEphemeris`](@ref)
"""
function designMatrixDM(
        toas,
        freqs;
        nMode=30, 
        TSpan=Nothing, 
        fRef=1400.0,
        logFreq=false
    )

    return designMatrixChromatic(
        toas, 
        freqs,
        nMode=nMode, 
        TSpan=TSpan, 
        fRef=fRef,
        index=2.0,
        logFreq=logFreq
    )
end


"""
    designMatrixEphemeris(
        toas,
        pos;
        nMode=30, 
        TSpan=Nothing, 
        logFreq=false
    )

Construct an ``N_{toas} × 6N_{mode}`` ephemeris perturbation Fourier design matrix. 
The matrix contains ``6N_{mode}`` columns, ordered as by frequency first, Cartesian coordinate second:
```math
    sin(f0)[x], sin(f0)[y], sin(f0)[z], cos(f0)[x], cos(f0)[y], cos(f0)[z],\\

    sin(f1)[x], sin(f1)[y], sin(f1)[z], cos(f1)[x], cos(f1)[y], cos(f1)[z],
    ...
```

The corresponding frequency vector repeats every entry six times.
This design matrix should be used with monopole_orf overlap function 
and with a powerlaw that specifies components=6.

# Arguments
- `toas::Array`: vector of time of arrial series in seconds.
- `pos:Array`: radio frequencies of observations [MHz].
- `nMode::Integer=30`: number of fourier coefficients to use.
- `TSpan:Float64=Nothing`: option to some other Tspan.
- `logFreq::Bool=false`: use log frequency spacing.

# Examples
```jldoctest
julia> toas = [1, 2]
julia> pos = [3, 4, 5]
julia> F, FFreqs = designMatrixEphemeris(toas, pos, nMode=2);
julia> F
2×12 Array{Float64,2}:
 -7.34788e-16  -9.79717e-16  -1.22465e-15  …  -2.44929e-15  3.0  4.0  5.0
 -1.46958e-15  -1.95943e-15  -2.44929e-15     -4.89859e-15  3.0  4.0  5.0
julia> FFreqs
12-element Array{Float64,1}:
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 2.0
 2.0
 2.0
 2.0
 2.0
 2.0
```

# See also
[`designMatrixRed`](@ref), [`designMatrixEnvelope`](@ref), [`designMatrixChromatic`](@ref), [`designMatrixDM`](@ref)
"""
function designMatrixEphemeris(
        toas,
        pos;
        nMode=30, 
        TSpan=Nothing, 
        logFreq=false
    )
    
    # get base fourier design matrix and frequencies
    F0, F0Freqs = designMatrixRed(toas, nMode=nMode, TSpan=TSpan, logFreq=logFreq)

    # frequency modes of design matrix F
    FFreqs = repeat(F0Freqs, inner=3)
    
    # F is a N_toa × 6N_mode design matrix
    F = zeros(length(toas), 6nMode)
    
    # set values by column
    for i in 1:nMode 
        F[:, 6i-5] .= F0[:, 2*i-1]*pos[1]
        F[:, 6i-4] .= F0[:, 2*i-1]*pos[2]
        F[:, 6i-3] .= F0[:, 2*i-1]*pos[3]
        F[:, 6i-2] .= F0[:, 2*i]*pos[1]
        F[:, 6i-1] .= F0[:, 2*i]*pos[2]
        F[:, 6i] .= F0[:, 2*i]*pos[3]
    end

    return F, FFreqs
end