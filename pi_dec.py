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

n = 51
a = j(PDC(n), PDC(n := n + 1)); print(a)

for x in range(1000):
    a = j(a, PDC(n := n + 1))


print(a)


314159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692346034861045432664821339360726024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414695194151160943305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194912983367336244065664308602139494639522473719070217986094370277053921717629317675238467481846766940513200056812714526356082778577134275778960917363717872146844090122495343014654958537105079227968925892354201995611212902196086403441815981362977477130996051870721134999999837297804995105973173281609631859502445945534690830264252230825334468503526193118817101000313783875288658753320838142061717766914730359825349042875546873115956286388235378759375195778185778053217122680661300192787661119590921642019893809525720106548586327886593615338182796823030195203530185296899577362259941389124972177528347913151557485724245415069595082953311686172785588907509838175463746493931925506040092770167113900984882401285836160356370766010471018194295559619894676783744944825537977472684710404753464620804668425906949129331367702898915210475216205696602405803815019351125338243003558764024749647326391419927260426992279678235478163600934172164121992458631503028618297455570674983850549458858692699569092721079750930295532116534498720275596023648066549911988183479775356636980742654252786255181841757467289097777279380008164706001614524919217321721477235014144197356854816136115735255213347574184946843852332390739414333454776241686251898356948556209921922218427255025425688767179049460165346680498862723279178608578438382796797668145410095388378636095068006422512520511739298489608412848862694560424196528502221066118630674427862203919494504712371378696095636437191728746776465757396241389086583264599581339047802759009946576407895126946839835259570982582262052248940772671947826848260147699090264013639443745530506820349625245174939965143142980919065925093722169646151570985838741059788595977297549893016175392846813826868386894277415599185592524595395943104997252468084598727364469584865383673622262609912460805124388439045124413654976278079771569143599770012961608944169486855584840635342207222582848864815845602850601684273945226746767889525213852254995466672782398645659611635488623057745649803559363456817432411251507606947945109659609402522887971089314566913686722874894056010150330861792868092087476091782493858900971490967598526136554978189312978482168299894872265880485756401427047755513237964145152374623436454285844479526586782105114135473573952311342716610213596953623144295248493718711014576540359027993440374200731057853906219838744780847848968332144571386875194350643021845319104848100537061468067491927819119793995206141966342875444064374512371819217999839101591956181467514269123974894090718649423196156794520809514655022523160388193014209376213785595663893778708303906979207734672218256259966150142150306803844773454920260541466592520149744285073251866600213243408819071048633173464965145390579626856100550810665879699816357473638405257145910289706414011097120628043903975951567715770042033786993600723055876317635942187312514712053292819182618612586732157919841484882916447060957527069572209175671167229109816909152801735067127485832228718352093539657251210835791513698820914442100675103346711031412671113699086585163983150197016515116851714376576183515565088490998985998238734552833163550764791853589322618548963213293308985706420467525907091548141654985946163718027098199430992448895757128289059232332609729971208443357326548938239119325974636673058360414281388303203824903758985243744170291327656180937734440307074692112019130203303801976211011004492932151608424448596376698389522868478312355265821314495768572624334418930396