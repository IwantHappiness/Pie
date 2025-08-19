import math

from time import sleep
from time import time

modulo = 0
inverse_modulo = 0


def initialize_modulo(m):
    global modulo, inverse_modulo

    modulo = m
    inverse_modulo = float(1) / float(m)


def mul_mod(a, b):
    q = int(inverse_modulo * float(a) * float(b))
    return a*b - q*modulo


def sum_mul_mod(a, b, c, d):
    q = int(inverse_modulo * (float(a)*float(b) + float(c)*float(d)))
    return a*b + c*d - q*modulo


def easy_round(x):
    full_double = float(2**53)

    y = x + full_double
    y -= full_double

    return y


def extended_gcd(a):
    A0, A1 = 1, 0
    r0, r1 = a, modulo

    while r1 > 0:
        q = r0 // r1

        A0, A1 = A1, A0 - q * A1
        r0, r1 = r1, r0 - q * r1

    return r0, A0


def inv_mod(a):
    a = a % modulo

    if a < 0:
        a += modulo

    gcd, A = extended_gcd(a)

    if gcd != 1:
        print("Error, gcd should be 1")
        exit()

    return A


def pow_mod(a, b):
    r = 1

    while True:
        if b & 1:
            r = mul_mod(r, a)

        b >>= 1
        if b == 0:
            break

        a = mul_mod(a, a)

    return r


def sum_binomial_mod(n, k):
    if k > n / 2:
        s = pow_mod(2, n) - sum_binomial_mod(n, n - k - 1)
        if s < 0:
            s += modulo

        return s

    prime_factors = []
    mm = modulo
    for p in range(3, int(math.sqrt(mm)) + 1, 2):
        if mm % p == 0:
            mm //= p

            if p <= k:
                prime_factors.append(p)

            while mm % p == 0:
                mm //= p

    if mm > 1 and mm <= k:
        prime_factors.append(mm)

    binomial_powers = [1] * len(prime_factors)
    next_denom = [pf for pf in prime_factors]
    next_num = [pf * (n // pf) for pf in prime_factors]

    binomial_num0, binomial_denom = 1, 1
    sum_num = 1
    binomial_secondary = 1

    for j in range(1, k + 1):
        # New binomial: b(n, j) = b(n, j - 1) * (n - j + 1) / j
        binomial_secondary_update = 0
        num = n - j + 1
        denom = j

        for i, p in enumerate(prime_factors):
            # Test if p is a prime factor of num0
            if next_num[i] == n - j + 1:
                binomial_secondary_update = 1
                next_num[i] -= p
                binomial_powers[i] *= p
                num //= p

                while num % p == 0:
                    binomial_powers[i] *= p
                    num //= p

            # Test if p is a prime factor of denom0
            if next_denom[i] == j:
                binomial_secondary_update = 1
                next_denom[i] += p
                binomial_powers[i] //= p
                denom //= p

                while denom % p == 0:
                    binomial_powers[i] //= p
                    denom //= p

        if binomial_secondary_update:
            binomial_secondary = 1
            for binomial_power in binomial_powers:
                binomial_secondary = mul_mod(binomial_secondary, binomial_power)

        binomial_num0 = mul_mod(binomial_num0, num)
        binomial_denom = mul_mod(binomial_denom, denom)

        if binomial_secondary != 1:
            sum_num = sum_mul_mod(sum_num, denom, binomial_num0, binomial_secondary)
        else:
            sum_num = mul_mod(sum_num, denom) + binomial_num0

    sum_num = mul_mod(sum_num, inv_mod(binomial_denom))
    return sum_num


def digits_of_fraction(n, a, b):
    initialize_modulo(b)
    power = pow_mod(10, n)
    c = mul_mod(power, a)

    return float(c) / float(b)


def digits_of_series(n, m):
    x = float()

    for k in range(0, m, 2):
        x += (
            digits_of_fraction(n, 4, 2 * k + 1) -
            digits_of_fraction(n, 4, 2 * k + 3)
        )

        x = x - easy_round(x)

    return x


def digits_of_pi(n):
    logn = math.log(float(n))
    M = 2 * int(3. * n / logn / logn / logn)  # M is even
    N = 1 + int(  # n >= N
        (n + 15.) * math.log(10.) / (1. + math.log(2. * M))
    )

    N += N % 2  # N should be even
    m_max = M * N + N

    x = digits_of_series(n, m_max)

    for k in range(N):
        m = 2 * M * N + 2 * k + 1
        initialize_modulo(m)

        s = sum_binomial_mod(N, k)
        s = mul_mod(s, pow_mod(5, N))
        s = mul_mod(s, pow_mod(10, n - N))  # n - N is always positive
        s = mul_mod(s, 4)

        x += (2 * (k % 2) - 1) * float(s) / float(m)  # 2*(k % 2) - 1 = (-1)^(k - 1)
        x = x - math.floor(x)

    return x


def PDC(n):
    assert n > 50, "Error, n should be bigger than 50"

    x = digits_of_pi(n)
    pow = 1.e9
    y = x * pow

    while pow > 10 and (y - math.floor(y) < 0.05 or y - math.floor(y) > 0.95):
        pow /= float(10)
        y = x * pow

    y = str(math.floor(y))
    return y


# Compute digits in base 10
if __name__ == "__main__":
    PAUSE_MIN_TIME = 0.05

    # PDC doesn't works before 51 digit
    n = 51
    prefix = "π = 3.141592653589793238462643383279502884197169399375105"

    # Start printing
    try:
        for key in prefix:
            print(key, end='', flush=True)
            sleep(PAUSE_MIN_TIME)

        trace = PDC(n)

        def trace_join(prev: str, new: str):
            for s in range(len(new)):
                if (pos := prev.find(new[:len(new) - s], -10)) != -1:
                    prev = prev[:pos] + new

                    return prev

        while True:
            timestamp = time()

            delta = time() - timestamp
            trace = trace_join(trace, PDC(n))

            if delta < PAUSE_MIN_TIME:
                sleep(PAUSE_MIN_TIME - delta)

            print(trace[n - 51], end='', flush=True)
            n += 1
            
    except KeyboardInterrupt:
        # Delete ^C character and print ellipsis
        print("\b\b", end='', flush=True)
        print("...")
