struct MReal <: Real
    ex
end

struct MComplex <: Number
    ex
end

Base.:(x::MReal, y::MReal) = MReal(weval(W"Times"(x.ex, y.ex)))
