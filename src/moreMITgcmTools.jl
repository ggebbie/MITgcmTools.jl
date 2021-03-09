#module moreMITgcmTools
#using MeshArrays,MITgcmTools, LaTeXStrings
#export extract_timeseries,matmul,position_label,nancount_gcmarray

"""
    extract_timeseries(froot,years,γ,xval,yval,fval)
# Arguments
- `froot::String`: filename root
- `years::StepRange`: iterator for multiple files
- `γ::gcmarray`: GCM grid from MeshArrays.jl
- `xval::Integer`: x grid index
- `yval::Integer`: y grid index
- `fval::Integer`: face index
# Output    
- `tseries`: timeseries
- `nseries`: nseries = number of timeseries elements in each year
# Examples
```jldoctest
julia> rand(3)
0.5 0.4 0.3
```
"""
function extract_timeseries(froot,years,γ,xval,yval,fval)
    tseries = []
    nseries = []
    nyr = length(years)
    for tt = 1:nyr
        fname = froot*string(years[tt])
        println(fname)
        field = read_bin(fname,Float32,γ)
        series = extract_timeseries(field,xval,yval,fval)
        push!(nseries,length(series))
        append!(tseries,series[:])
    end
    return tseries,nseries
end

################################################################################################
"""
    extract_timeseries(flux,xval,yval,fval)
# Arguments
- `flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: input array of arrays
- `years::StepRange`: iterator for multiple files
- `xval::Integer`: x grid index
- `yval::Integer`: y grid index
- `fval::Integer`: face index
# Output    
- `series`: timeseries
"""
function extract_timeseries(flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},xval,yval,fval)
    nt = size(flux,2)
    series = Float32[]
    
    for tval ∈ 1:nt
        tmp = flux[fval,tval]
        push!(series,tmp[xval,yval])
    end
    return series
end

###############################################################################################
"""
    extract_timeseries(flux,xval,yval,fval)
# Arguments
- `flux::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}`: input array of arrays
- `years::StepRange`: iterator for multiple files
- `xval::Integer`: x grid index
- `yval::Integer`: y grid index
- `fval::Integer`: face index
# Output    
- `series`: timeseries
"""
function extract_timeseries(flux::MeshArrays.gcmarray{Float64,2,Array{Float64,2}},xval,yval,fval)
    nt = size(flux,2)
    series = Float64[]
    
    for tval ∈ 1:nt
        tmp = flux[fval,tval]
        push!(series,tmp[xval,yval])
    end
    return series
end

#################################################################################################
# tell julia how to do matrix multiplication
function matmul(M::Array{Float64,2},flux::MeshArrays.gcmarray{Float64,2,Array{Float64,2}},γ) ::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}

    nM = size(M,1)  # matrix size M
    nN = size(flux,2) # matrix size N
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    product = 0.0 .* MeshArray(γ,Float64,nM)
    
    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            for nn = 1:nN # inner product over time
                product[qq,mm] += flux[qq,nn] * M[mm,nn]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
            end
        end
    end
    return product
end

function matmul(M::Array{Float32,2},flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},γ::gcmgrid) ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}

    nM = size(M,1)  # matrix size M
    nN = size(flux,2) # matrix size N
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0f0 .* MeshArray(γ,Float32,nM)
    println("sub-optimal initialization")
    product = MeshArray(γ,Float32,nM) # some nans here
    tmp1=zeros(Float32,Tuple(γ.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=γ.read(tmp1,MeshArray(γ,Float32))
    end

    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            for nn = 1:nN # inner product over time
                product[qq,mm] += flux[qq,nn] * M[mm,nn]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
            end
        end
    end
    return product
end

function matmul(M::Array{Float32,1},flux::MeshArrays.gcmarray{Float32,1,Array{Float32,2}},γ::gcmgrid) ::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}

    nM = size(M,1)  # matrix size M
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0f0 .* MeshArray(γ,Float32,nM)
    println("sub-optimal initialization")
    product = MeshArray(γ,Float32,nM) # some nans here
    tmp1=zeros(Float32,Tuple(γ.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=γ.read(tmp1,MeshArray(γ,Float32))
    end

    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            product[qq,mm] += flux[qq] * M[mm]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
        end
    end
    return product
end

function matmul(M::Array{Float64,1},flux::MeshArrays.gcmarray{Float64,1,Array{Float64,2}},γ) ::MeshArrays.gcmarray{Float64,2,Array{Float64,2}}

    nM = size(M,1)  # matrix size M
    nQ  = size(flux,1) # repeat matmul Q times
    
    # initialize product
    #product = 0.0 .* MeshArray(γ,Float64,nM)
    println("sub-optimal initialization")
    product = MeshArray(γ,Float64,nM) # some nans here
    tmp1=zeros(Float64,Tuple(γ.ioSize))
    for tt= 1:nM
        # initialize a sub-optimal way
        product[:,tt]=γ.read(tmp1,MeshArray(γ,Float64))
    end
    
    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            product[qq,mm] += flux[qq] * M[mm]
                #product[ff,pp] +=  transpose(F[pp,tt]'*flux[ff,tt]')
        end
    end
    return product
end

function position_label(lon,lat)
    # produce label for a title by rounding to nearest whole integer.
    if lat >= 0
        latlbl = string(round(Integer,lat))* L"{\degree}N"
    else
        latlbl = string(round(Integer,-lat))* L"{\degree}S"
    end

    if lon >= 0
        lonlbl = string(round(Integer,lon))* L"{\degree}E"
    else
        lonlbl = string(round(Integer,-lon))* L"{\degree}W"
    end

    lbl = latlbl * " " * lonlbl
    return lbl
end

function nancount_gcmarray(field)
    # function nancount_gcmarray(field)
    s1,s2 = size(field)

    nancount = zeros(s1,s2)
    for i1 = 1:s1
        for i2 = 1:s2
            nancount[i1,i2] = sum(isnan,field[i1,i2])
        end
    end
    return nancount
end

end
