using YaoBase, SparseArrays, BitBasis, YaoArrayRegister, SymEngine
export @ket_str, @bra_str

YaoArrayRegister._warn_type(raw::AbstractArray{Basic}) = nothing

const SymReg{B, MT} = ArrayReg{B, Basic, MT}
const AdjointSymReg{B, MT} = AdjointArrayReg{B, Basic, MT}
const SymOrAdjointSymReg{B, MT} = Union{SymReg{B, MT}, AdjointSymReg{B, MT}}

function parse_str(s::String)
    v = 0; k = 1
    for each in reverse(filter(x->x!='_', s))
        if each == '1'
            v += 1 << (k - 1)
            k += 1
        elseif each == '0'
            k += 1
        elseif each == '_'
            continue
        else
            error("expect 0 or 1, got $each at $k-th bit")
        end
    end
    return v, k-1
end

function ket_m(s)
    v, N = parse_str(s)
    st = spzeros(Basic, 1 << N, 1)
    st[v+1] = 1
    return ArrayReg{1}(st)
end

function bra_m(s)
    adjoint(ket_m(s))
end

"""
    @ket_str

Create a ket register. See also [`@bra_str`](@ref).

# Example

a symbolic quantum state can be created simply by

```jldoctest
julia> ket"110" + 2ket"111"
|110⟩ + 2|111⟩
```

qubits can be partially actived by [`focus!`](@ref)

```jldoctest
```

"""
macro ket_str(s)
    ket_m(s)
end

macro bra_str(s)
    bra_m(s)
end

function SymEngine.expand(x::SymReg{B}) where B
    ArrayReg{B}(expand.(x.state))
end

function Base.show(io::IO, r::SymReg{1})
    dropzeros!(r.state)
    rows = rowvals(r.state)
    nnz = nonzeros(r.state)
    if size(r.state, 2) == 1 # all actived
        nzr = nzrange(r.state, 1)
        for i in nzr
            k = rows[i]
            v = nnz[i]
            if !isone(v)
                print(io, v)
            end
            print(io, "|", string(k-1, base=2, pad=nactive(r)), "⟩")

            if i != last(nzr)
                print(io, " + ")
            end
        end
    else
        m, n = size(r.state)
        nnz_count = 0
        nnnz = length(nonzeros(r.state))

        for j in 1:n
            nzr = nzrange(r.state, j)
            for i in nzr
                row = rows[i]
                val = nnz[i]

                if !isone(val)
                    print(io, val)
                end
                print(io, "|")
                printstyled(io, string(j-1, base=2, pad=nremain(r)), color=:light_black)
                print(io, string(row-1, base=2, pad=nactive(r)), "⟩")
                nnz_count += 1
                if nnz_count != nnnz
                    print(io, " + ")
                end    
            end
        end
    end
end

function Base.show(io::IO, r::AdjointSymReg{1})
    dropzeros!(parent(state(r)))
    rows = rowvals(parent(state(r)))
    nnz = nonzeros(parent(state(r)))
    nzr = nzrange(parent(state(r)), 1)
    if size(parent(state(r)), 2) == 1 # all actived
        for i in nzr
            k = rows[i]
            v = adjoint(nnz[i])
            if !isone(v)
                print(io, v)
            end

            print(io, "⟨", string(k-1, base=2, pad=nactive(r)), "|")

            if i != last(nzr)
                print(io, " + ")
            end
        end
    else
        m, n = size(state(r))
        nnz_count = 0
        nnnz = length(nonzeros(r.state))

        for j in 1:n
            nzr = nzrange(r.state, j)
            for i in nzr
                row = rows[i]
                val = nnz[i]

                if !isone(v)
                    print(io, v)
                end

                print(io, "⟨")
                print(io, string(row-1, base=2, pad=nactive(r)))
                printstyled(io, string(j-1, base=2, pad=nremain(r)), color=:light_black)
                print(io, "|")

                if i != last(nzr)
                    print(io, " + ")
                end
            end
        end
    end
end

Base.:(*)(x::SymReg{B, MT}, y::SymReg{B, MT}) where {B, MT} = SymReg{B, MT}(kron(state(x), state(y)))
Base.:(^)(x::SymReg{B, MT}, n::Int) where {B, MT} = SymReg{B, MT}(kron(state(x) for _ in 1:n))

Base.:(*)(x::AdjointSymReg{B, MT}, y::AdjointSymReg{B, MT}) where {B, MT} = adjoint(parent(x) * parent(y))
Base.:(^)(x::AdjointSymReg{B, MT}, n::Int) where {B, MT} = adjoint(parent(x)^n)