module BilinearCheck

using Random

import Random: rand

export randComplex,inner,isLeftLinear,isRightLinear,isHomogenousBi,
       isHomogenousSus,nonDegenFirst,nonDegenSecond,nonDegen,
       isSymmetric,isHermitian,isBilinear,isSusquilinear

function randComplex(max)
    r = rand(Int64,3,1) .% max
    i = rand(Int64,3,1) .% max
    return r .+ i.*im
end

inner(v,w) = (w'*v)[1]

function isLeftLinear(β,a,b,c)
    β(a+b,c) == β(a,c) + β(b,c)
end

function isRightLinear(β,a,b,c)
    β(c,a+b) == β(c,a) + β(c,b)
end

function isHomogenousBi(β,a,b,λ)
    (β(λ*a,b) == λ*β(a,b)) && (β(a,λ*b) == λ*β(a,b))
end

function isHomogenousSus(β,a,b,λ)
    (β(λ*a,b) == λ*β(a,b)) && (β(a,λ*b) == conj(λ)*β(a,b))
end

function nonDegenFirst(β,a,b)
    a = zeros(size(a))
    norm(β(a,b)) == 0
end

function nonDegenSecond(β,a,b)
    b = zeros(size(b))
    norm(β(a,b)) == 0
end

function nonDegen(β,a,b)
    return nonDegenFirst(β,a,b) && nonDegenSecond(β,a,b)
end

function isSymmetric(β,a,b)
    β(a,b) == β(b,a)
end

function isHermitian(β,a,b)
    β(a,b) == conj(β(b,a))
end

function isBilinear(β,a,b,c,λ)
    return isLeftLinear(β,a,b,c) && isRightLinear(β,a,b,c) && isHomogenousBi(β,a,b,λ)
end

function isSusquilinear(β,a,b,c,λ)
    return isLeftLinear(β,a,b,c) && isRightLinear(β,a,b,c) && isHomogenousSus(β,a,b,λ)
end

end #module BilinearCheck
