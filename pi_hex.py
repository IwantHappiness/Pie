from math import floor
from time import time
from time import sleep


def mod(x, n):
    x %= n
    return  x + n if x < 0 else x


def modPow(b, e, m):
    b = b % m;
    y = 1

    while (e > 0):
        if (e & 1): y = (y * b) % m
        b = (b * b) % m
        e >>= 1

    return y


def S(j, n):
    # Left sum
    left = 0
    for k in range(0, n + 1):
        r = 8 * k + j
        left = mod((left + modPow(16, n - k, r) / r), 1)

    # Right sum
    right = 0
    k = n + 1
    while True:
        rnew = right + 16 ** (n - k) / (8 * k + j)
        if right == rnew: break
        right = rnew

        k += 1

    return left + right


def BBP(d, n):
    d -= 1

    return hex(
        floor(
            16 ** n * mod(4 * S(1, d) - 2 * S(4, d) - S(5, d) - S(6, d), 1)
        )
    )[2:]  # Remove 0x prefix


# Compute digits in base 16
if __name__ == "__main__":
    PAUSE_MIN_TIME = 0.05

    n = 1
    prefix = "π = 3."

    # Start printing
    try:
        for key in prefix:
            print(key, end='', flush=True)
            sleep(PAUSE_MIN_TIME)

        while True:
            timestamp = time()
            trace = BBP(n, 1)

            delta = time() - timestamp
            if delta < PAUSE_MIN_TIME:
                sleep(PAUSE_MIN_TIME - delta)

            print(trace, end='', flush=True)
            n += 1

    except KeyboardInterrupt:
        # Delete ^C character and print ellipsis
        print("\b\b", end='', flush=True)
        print("...")
