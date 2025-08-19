#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <time.h>

// Formatted original pidec code from here: http://numbers.computation.free.fr/Constants/Algorithms/pidec.cpp

// longint should be a 64-bit long integer
// it could be a double floating point type, provided
// adaptations of mulMod and sumMulMod functions are done
typedef __int64_t longint;

longint modulo;
double inverse_modulo;

void initializeModulo(longint m) {
    modulo = m;
    inverse_modulo = 1. / (double) m;
}

double myTime() {
    return ((double) clock()) / CLOCKS_PER_SEC;
}

// Compute a*b % modulo
inline longint mulMod(longint a, longint b) {
    // Classical trick to bypass the 64-bit limitation, when a*b does not fit into the longint type
    // Works whenever a * b / modulo is less than 2^52 (double type maximal precision)
    longint q = (longint) (inverse_modulo * (double) a * (double) b);

    return a*b - q*modulo;
}

// Compute a*b + c*d % modulo
inline longint sumMulMod(longint a, longint b, longint c, longint d) {
    longint q = (longint) (inverse_modulo * ((double) a * (double) b + (double) c * (double) d));
    return a*b + c*d - q*modulo;
}

double fullDouble = 1024. * 1024. * 1024. * 1024. * 1024. * 8.;  // 2^53
inline double easyRound(double x) {
    double y = x + fullDouble;

    y -= fullDouble;
    return y;
}

// Return g, A such that g = gcd(a, modulo) and a * A = g % modulo
longint extendedGcd(longint a, longint &A) {
    longint A0 = 1, A1 = 0;
    longint r0 = a, r1 = modulo;

    while (r1 > 0.) {
        longint q = r0 / r1;

        longint tmp = A0 - q * A1;
        A0 = A1;
        A1 = tmp;

        tmp = r0 - q * r1;
        r0 = r1;
        r1 = tmp;
    }

    A = A0;
    return r0;
}

longint invMod(longint a) {
    longint A;

    a = a % modulo;
    if (a < 0)
        a += modulo;

    longint gcd = extendedGcd(a, A);

    if (gcd != 1) {
        printf("pb, gcd should be 1\n");
        exit(1);
    }

    return A;
}

longint powMod(longint a, long b) {
    longint r = 1;

    while (true) {
        if (b & 1)
            r = mulMod(r, a);

        b >>= 1;
        if (b == 0)
            break;

        a = mulMod(a, a);
    }

    return r;
}

// Compute sum_{j = 0}^k binomial(n, j) mod m
longint sumBinomialMod(long n, long k) {
    // Optimisation: when k>n/2 we use the relation
    // sum_{j = 0}^k binomial(n, j) = 2^n - sum_{j = 0}^{n - k - 1} binomial(n, j)

    // Note: additionnal optimization, not afforded here, could be done when k is near n / 2
    // using the identity sum_{j = 0}^{n / 2} = 2^(n - 1) + 0.5 binomial(n, n / 2)
    // A global saving of 20% or 25% could be obtained
    if (k > n / 2) {
        longint s = powMod(2, n) - sumBinomialMod(n, n - k - 1);
        if (s < 0)
            s += modulo;

        return s;
    }

    // Compute prime factors of modulo which are smaller than k
    const long nbMaxFactors = 20;  // No more than 20 different prime factors for numbers < 2^64
    long primeFactor[nbMaxFactors];
    long nbPrimeFactors = 0;
    longint mm = modulo;

    // Modulo is odd, thus has only odd prime factors
    for (longint p = 3; p * p <= mm; p += 2) {
        if (mm % p == 0) {
            mm = mm / p;

            if (p <= k)  // Only prime factors <= k are needed
                primeFactor[nbPrimeFactors++] = p;

            while (mm % p == 0)
                mm = mm / p;  // Remove all powers of p in mm
        }
    }

    // Last factor: if mm is not 1, mm is necessarily prime
    if (mm > 1 && mm <= k) {
        primeFactor[nbPrimeFactors++] = mm;
    }

    long binomialPowers[nbMaxFactors];  // Powers of primeFactor[i] in binomial
    long nextDenom[nbMaxFactors], nextNum[nbMaxFactors];  // Next multiples of primeFactor[i]

    for (long i = 0; i < nbPrimeFactors; i++) {
        binomialPowers[i] = 1;
        nextDenom[i] = primeFactor[i];
        nextNum[i] = primeFactor[i] * (n / primeFactor[i]);
    }

    longint binomialNum0 = 1, BinomialDenom = 1;
    longint sumNum = 1;
    longint binomialSecondary = 1;

    for (long j = 1; j <= k; j++) {
        // New binomial: b(n, j) = b(n, j - 1) * (n - j + 1) / j
        int binomialSecondaryUpdate = 0;
        longint num = n - j + 1;
        longint denom = j;

        for (long i = 0; i < nbPrimeFactors; i++) {
            long p = primeFactor[i];

            // Test if p is a prime factor of num0
            if (nextNum[i] == n - j + 1) {
                binomialSecondaryUpdate = 1;
                nextNum[i] -= p;
                binomialPowers[i] *= p;
                num /= p;

                while (num % p == 0) {
                    binomialPowers[i] *= p;
                    num /= p;
                }
            }

            // Test if p is a prime factor of denom0
            if (nextDenom[i] == j) {
                binomialSecondaryUpdate = 1;
                nextDenom[i] += p;
                binomialPowers[i] /= p;
                denom /= p;

                while (denom % p == 0) {
                    binomialPowers[i] /= p;
                    denom /= p;
                }
            }
        }

        if (binomialSecondaryUpdate) {
            binomialSecondary = binomialPowers[0];
            for (long i = 1; i < nbPrimeFactors; i++)
                binomialSecondary = mulMod(binomialSecondary, binomialPowers[i]);
        }

        binomialNum0 = mulMod(binomialNum0, num);
        BinomialDenom = mulMod(BinomialDenom, denom);

        if (binomialSecondary != 1) {
            sumNum = sumMulMod(sumNum, denom, binomialNum0, binomialSecondary);
        }
        else {
            sumNum = mulMod(sumNum, denom) + binomialNum0;
        }
    }

    sumNum = mulMod(sumNum, invMod(BinomialDenom));
    return sumNum;
}

// Return fractionnal part of 10^n*(a/b)
double digitsOfFraction(long n, longint a, longint b) {
    initializeModulo(b);
    longint pow = powMod(10, n);
    longint c = mulMod(pow, a);

    return (double) c / (double) b;
}

// Return fractionnal part of 10^n*S where S = 4 * sum_{k = 0}^{m - 1} (-1)^k/(2*k + 1); m is even
double digitsOfSeries(long n, longint m) {
    double x = 0.;

    for (longint k = 0; k < m; k += 2) {
        x += (
            digitsOfFraction(n, 4, 2 * k + 1) -
            digitsOfFraction(n, 4, 2 * k + 3)
        );
        x = x - easyRound(x);
    }

    return x;
}

double digitsOfPi(long n) {
    double logn = log((double) n);
    long M = 2 * (long) (3. * n / logn / logn / logn);  // M is even
    long N = 1 + (long) ((n + 15.) * log(10.) / (1. + log(2. * M)));  // n >= N

    N += N % 2;  // N should be even
    longint mMax = (longint) M * (longint) N + (longint) N;

    printf("Parameters: M = %ld, N = %ld, M*N + M = %.0lf\n", M, N, (double) mMax);

    double timestamp = myTime();
    double x = digitsOfSeries(n, mMax);
    printf("Series time: %.2lf\n", myTime() - timestamp);

    for (long k = 0.; k < N; k++) {
        longint m = (longint) 2 * (longint) M * (longint) N + (longint) 2 * (longint) k + 1;
        initializeModulo(m);

        longint s = sumBinomialMod(N, k);
        s = mulMod(s, powMod(5, N));
        s = mulMod(s, powMod(10, n - N)); // n-N is always positive
        s = mulMod(s, 4);

        x += (2 * (k % 2) - 1) * (double) s / (double) m;  // 2*(k%2) - 1 = (-1)^(k - 1)
        x = x - floor(x);
    }

    return x;
}

int main() {
    printf("Pidec, direct computation of decimal digits of pi at a given position n\n");
    printf("(http://numbers.computation.free.fr/Constants/constants.html for more details)\n");

    long n;
    printf("Enter n: "); scanf("%ld", &n);

    if (n < 50) {
        printf("Error, n should be bigger than 50. Please retry\n");
        exit(1);
    }

    double timestamp = myTime();
    double x = digitsOfPi(n);
    double pow = 1.e9;
    double y = x * pow;

    // To be able to know exactly the digits of pi at position n,
    // the value (pow * x) should be not too close to an integer
    while (pow > 10 && (y - floor(y) < 0.05 || y - floor(y) > 0.95)) {
        pow /= 10.;
        y = x * pow;
    }

    printf("Digits of pi after n-th decimal digit: %.0lf\n", floor(y));
    printf("Total time: %.2lf\n", myTime() - timestamp);
}
