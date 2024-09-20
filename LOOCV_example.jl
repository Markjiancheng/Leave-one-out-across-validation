using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Distributions")
Pkg.add("Plots")
Pkg.add("LinearAlgebra")
Pkg.add("Statistics")
Pkg.add("Random")
Pkg.add("DelimitedFiles")

using Distributions
using LinearAlgebra
using Random
using CSV
using DataFrames
using SparseArrays

# run basic functions for LOOCV
# for fixed effect
function getX(factor,df;cov=false)
    data = df[:,Symbol(factor)]
    n = size(data,1)
    if cov==false
        str = data
        val = 1.0
    else 
        str = fill(factor,n)
        val = data
    end

    dict,names   = mkDict(str)
    ii    = 1:n                    # row numbers 
    jj   = [dict[i] for i in str]  # column numbers
    X    = sparse(ii,jj,val)    
end  

# Function for LOOCV
function loiEBVEff(i,Vi,C,y)
    r     = 1/Vi[i,i]
    q     = -Vi[i,:] *r
    yi    = copy(y)
    yi[i] = 0.0
    Ci    = C[i,:]
    Ci[i] = 0.0
    (Ci'*Vi - (Ci'q) * q'/r)*yi
end




# read data file
phenotypes = CSV.read("Leave-one-out-across-validation/data/example.dat", DataFrame, types=Dict(:ID => String, :Died => String), delim = ' ',header=true, missingstrings=["NA"])

# read geno file, missing geno not allowed
M= CSV.read("Leave-one-out-across-validation/data/example_geno.txt",DataFrame,types=Dict(:1 => String),
           delim = ' ',header=false, missingstrings=["NA"] )
           
# change format to matrix
M2= Matrix(M[!,2:size(M, 2)])

# remove missing individuals
ZFull = Matrix{Float64}(I,size(phenotypes, 1),size(phenotypes, 1))
res = phenotypes[:,:Nur2ADG]
phenotypes=dropmissing(phenotypes, :Nur2ADG)

sel = .!(ismissing.(res))
Z = ZFull[sel,:]
     
# Add genetic variance, sum2pq
varG = 0.295395E-02
sum_2pq=165899.511107207
varAlpha = varG/sum_2pq
@time G = Z*M2*M2'*Z'*varAlpha


# treat fixed effect as random with large variance, because residual var/ batch var is too small, nothing was added to diagonal.
X1 = getX("Batch",phenotypes)
X2 = getX("Died",phenotypes)
X3 = getX("EntryAge",phenotypes,cov=true)

var_Batch = 1
var_Died  = 1
var_Age=0.103
X=var_Batch*X1*X1' + var_Died*X2*X2' + var_Age*X3*X3'
Matrix(X)

# random effect of pen and sow
Z1 = getX("NurPenBatch",phenotypes)
Matrix(Z1*Z1')
Z2 = getX("SowID",phenotypes)
Matrix(Z2*Z2')

var_pen= 0.107150E-02
var_sow= 0.927471E-03
U=var_pen*Z1*Z1' + var_sow*Z2*Z2'

# residual variance
varRes = 0.115504E-01

# Total variance
V = G + U + X + I*varRes


# input for LOOCV run
Vi = inv(V)
C = G
n = size(V,1)
y  = phenotypes[!,:Nur2ADG]


# loocv run
@time EBV=[loiEBVEff(i,Vi,C,y) for i=1:n]

# check correlation
cor(EBV, y)


















