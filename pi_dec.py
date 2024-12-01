import math

from time import sleep
from time import time

from time import process_time

_full_double = 1024. * 1024. * 1024. * 1024. * 1024. * 8.  # 2^53
_m = 0
_invm = 0

def initialize_modulo(m):
    global _m, _invm
    _m = m
    _invm = 1. / float(m)

def mul_mod(a, b):
    q = int(_invm * float(a) * float(b))
    return a * b - q * _m

def sum_mul_mod(a, b, c, d):
    q = int(_invm * (float(a) * float(b) + float(c) * float(d)))
    return a * b + c * d - q * _m

def my_time():
    return process_time()

def easy_round(x):
    y = x + _full_double
    y -= _full_double
    return y

def extended_gcd(a):
    A0, A1 = 1, 0
    r0, r1 = a, _m

    while r1 > 0:
        q = r0 // r1

        A0, A1 = A1, A0 - q * A1
        r0, r1 = r1, r0 - q * r1

    return r0, A0

def inv_mod(a):
    a = a % _m
    if a < 0:
        a += _m
    gcd, A = extended_gcd(a)
    if gcd != 1:
        print("Error, gcd should be 1")
    return A

def pow_mod(a, b):
    r = 1
    aa = a
    while True:
        if b & 1:
            r = mul_mod(r, aa)
        b >>= 1
        if b == 0:
            break
        aa = mul_mod(aa, aa)
    return r

def sum_binomial_mod(n, k):
    if k > n / 2:
        s = pow_mod(2, n) - sum_binomial_mod(n, n - k - 1)
        if s < 0:
            s += _m
        return s

    prime_factors = []
    mm = _m
    for p in range(3, int(math.sqrt(mm)) + 1, 2):
        if mm % p == 0:
            mm //= p
            if p <= k:
                prime_factors.append(p)
            while mm % p == 0:
                mm //= p

    if mm > 1 and mm <= k:
        prime_factors.append(mm)

    binomial_power = [1] * len(prime_factors)
    next_denom = [pf for pf in prime_factors]
    next_num = [pf * (n // pf) for pf in prime_factors]

    binomial_num0 = 1
    binomial_denom = 1
    sum_num = 1
    binomial_secondary = 1

    for j in range(1, k + 1):
        num = n - j + 1
        denom = j
        binomial_secondary_update = 0

        for i, p in enumerate(prime_factors):
            if next_num[i] == n - j + 1:
                binomial_secondary_update = 1
                next_num[i] -= p
                binomial_power[i] *= p
                num //= p
                while num % p == 0:
                    binomial_power[i] *= p
                    num //= p
            if next_denom[i] == j:
                binomial_secondary_update = 1
                next_denom[i] += p
                binomial_power[i] //= p
                denom //= p
                while denom % p == 0:
                    binomial_power[i] //= p
                    denom //= p

        if binomial_secondary_update:
            binomial_secondary = 1
            for bp in binomial_power:
                binomial_secondary = mul_mod(binomial_secondary, bp)

        binomial_num0 = mul_mod(binomial_num0, num)
        binomial_denom = mul_mod(binomial_denom, denom)

        if binomial_secondary != 1:
            sum_num = sum_mul_mod(sum_num, denom, binomial_num0, binomial_secondary)
        else:
            sum_num = mul_mod(sum_num, denom) + binomial_num0

    sum_num = mul_mod(sum_num, inv_mod(binomial_denom))
    return sum_num

def digits_of_fraction(n, a, b):
    initialize_modulo (b)
    pow = pow_mod(10, n)
    c = mul_mod(pow, a)
    return float(c) / float(b)

def digits_of_series(n, m):
    x = 0.
    for k in range(0, m, 2):
        x += digits_of_fraction(n, 4, 2 * k + 1) - digits_of_fraction(n, 4, 2 * k + 3)
        # print(f"ER: {easy_round(x):.0f}")
        x = x - easy_round(x)
    return x

def digits_of_pi(n):
    logn = math.log(float(n))
    M = 2 * int(3. * n / logn / logn / logn)  # M is even
    N = 1 + int((n + 15.) * math.log(10.) / (1. + math.log(2. * M)))  # n >= N
    N += N % 2  # N should be even
    mmax = M * N + N
    # print(f"Parameters : M={M}, N={N}, M*N+M={mmax:.0f}")
    st = my_time()
    x = digits_of_series(n, mmax)
    # print(f"Series time : {my_time() - st:.2f}")
    for k in range(N):
        m = 2 * M * N + 2 * k + 1
        initialize_modulo(m)
        s = sum_binomial_mod(N, k)
        s = mul_mod(s, pow_mod(5, N))
        s = mul_mod(s, pow_mod(10, n - N))  # n-N is always positive
        s = mul_mod(s, 4)
        x += (2 * (k % 2) - 1) * float(s) / float(m)  # 2*(k%2)-1 = (-1)^(k-1)
        x = x - math.floor(x)
    return x

# def main():
#     n = int(input("Pidec, direct computation of decimal digits of pi at a given position n.\n"
#                   "(http://numbers.computation.free.fr/Constants/constants.html for more details)\n"
#                   "Enter n : "))
#     if n < 50:
#         print("Error, n should be bigger than 50. Please retry")
#         return
#     st = my_time()
#     x = digits_of_pi(n)
#     pow = 1.e9
#     y = x * pow
#     while pow > 10 and (y - math.floor(y) < 0.05 or y - math.floor(y) > 0.95):
#         pow /= 10.
#         y = x * pow
#     print(f"Digits of pi after n-th decimal digit : {math.floor(y):.0f}")
#     print(f"Total time: {my_time() - st:.2f}")

# 


def PDC(n):
    assert n > 50, "N should be bigger than 50"

    x = digits_of_pi(n)
    pow = 1.e9
    y = x * pow
    while pow > 10 and (y - math.floor(y) < 0.05 or y - math.floor(y) > 0.95):  # 0.83 Ultra
        pow /= 10.
        y = x * pow

    # print(f"РАСХ: {math.floor(y)} {y} {pow=} ∂{abs(y - math.floor(y))}")
    l = len(str(math.floor(y)))
    y = str(math.floor(y))

    return y

# PDC(218)
    # (" " * (n - 50)) + (
    #     "\033[31m0\033[0m" + y if l !=9 else y
    # ) + f"    {l} {n=}\n "


first_50 = "141592653589793238462643383279502884197169399375105820974944592307816406286208998"
trace =                                                       "82097494"

# Compute n'th digit in base 10
PAUSE_MIN_TIME = 0.05

n = -1
print("π = 3.", end="", flush=True)
shift = -1

try:
    while True:
        now = time()

        n += 1
        nth = PDC(n) if n > 50 else first_50[n]

        
        nth_crop = nth[:shift]


        if n > 50:
            p = trace.find(nth_crop, -12 + shift)

            if p != -1:
                # print("\n", trace, nth, nth_crop)
                # exit()

                trace = trace[:p] + nth
                shift = -1
            else:
                shift -= 1

                
        delta = time() - now
        if delta < PAUSE_MIN_TIME:
            sleep(PAUSE_MIN_TIME - delta)

        if n > 50:
            print(trace[n - 51], end='', flush=True)
        else:
            print(first_50[n], end='', flush=True)

except KeyboardInterrupt:
    # Delete ^C character and print ellipsis
    print("\b\b", end='', flush=True)
    print("...")

except IndexError as e:
    print(e)

    print(f"{trace=} {n=} {nth=} {shift=}")

# 230781640
# +307816406 1
#  ++78164062 1
#    781640628 1
#    +816406286 1
#     +164062862 1
#      +640628620 1
#       +40628620 0 (2)
#        ++628620 0 (3)
#          628620899 3

def j(prev: str, new: str):
    for s in range(len(new)):
        if (p := prev.find(new[:len(new) - s], -10)) != -1:
            prev = prev[:p] + new

            return prev

# This is solution
n = 51
a = j(PDC(n), PDC(n := n + 1)); print(a)

for x in range(1000):
    a = j(a, PDC(n := n + 1))

print(a)

