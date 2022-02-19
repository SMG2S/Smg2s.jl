using SMG2S, Test

#Define size of matrix to be generated
size = 100

#maximum number of continous `1` in the diagonal of generated nilpotent matrix
nbOne = 4
#offsets of lower diagonals between which the matrix will be initially filled
diag_l = -20
diag_u = -5

#a function to generate the user-provided spectrum, `idx` is the indexing of eigenvalues
function f(idx::Integer, size = size)
    return 10 * cos(((idx-1) / 2) * 2 * π / size) + 5
end

#create a empty vector to store the generated eigenvalues
#Attention that this vector should always in complex scalar type
spec = zeros(ComplexF64, size)

#Setting up the spectrum with function `f`

Spectrum!(spec, f, size)

#Generating a non-symmetric sparse matrix and store it into `genMat`
genMat = nonsym(nbOne, size, diag_l, diag_u, spec);
