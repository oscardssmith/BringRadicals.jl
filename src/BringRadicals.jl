module BringRadicals
export bringrad

bringrad(a::Integer) = bringrad(float(a))
function bringrad(a::T) where T<Real
    x = abs(a)
    if x < .22
        at4 = a^4
        return a*evalpoly(at4, T(1, -1, 5, -35, 285))
    end
    if x<2
        guess = x*evalpoly(x, T(1.0098200485727955, -0.1885556510194153, 1.2609058087320157,
                              -3.7000253271311796, 3.4819520619583963, -1.109350971521824))
    else
        guess = @fastmath exp2(log2(x)*inv(T(5))))-inv(x+10)
    end
    tol = sqrt(eps(guess))
    guess = copysign(guess, -a)
    while true
        gt4 = (guess^2)^2 # avoid expense of floating point pow
        newguess = guess - muladd(guess, gt4, guess+a)/muladd(5, gt4, 1)
        abs(newguess - guess) <= tol && (guess = newguess; break)
        guess = newguess
    end
    gt4 = guess^4
    return guess - fma(guess, gt4, guess+a)/fma(5, gt4, 1)
end
end
