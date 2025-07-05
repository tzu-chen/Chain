module Utils

using ITensors

# Translation operator for periodic chains
function translation_op(sites::Vector{Index}; inverse::Bool=false)
    N = length(sites)
    ops = Vector{ITensor}(undef, N-1)
    for j in 1:N-1
        ops[j] = BondGate(SiteSet(sites), j, j+1)
    end
    if inverse
        ops = reverse(ops)
    end
    A = Vector{ITensor}(undef, N-1)
    B = Vector{ITensor}(undef, N-1)
    for j in 1:N-1
        A[j], B[j] = factor(ops[j], (sites[j], prime(sites[j])))
    end
    mpo = MPO(N)
    mpo[1] = A[1]
    for j in 2:N-1
        mpo[j] = B[j-1] * A[j]
    end
    mpo[N] = B[N-1]
    return mpo
end

function delta3_itensor(s1::Index, s2::Index, s3::Index)
    t = ITensor(s1, s2, s3)
    for j in 1:dim(s1)
        t[s1(j), s2(j), s3(j)] = 1
    end
    return t
end

end
