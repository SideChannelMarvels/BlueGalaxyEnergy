
Q.<y> = GF(2)[]
P.<x> = GF(2^8, name='x', modulus=y^8+y^4+y^3+y+1)
assert (x^6+x^4+x^2+x+1) * (x^7+x+1) == x^7+x^6+1

# p is a polynomial of P
def poly_to_int(p):
    return sum([1 << power for power, value in enumerate(p.polynomial().coefficients(sparse=False)) if value])

# s is an integer
def int_to_poly(s):
    L = Integer(s, 16)
    k = 0
    p = P(0)
    while L != 0:
        p = p + GF(2)(L & 1) * x^k
        k += 1
        L >>= 1
    return p

P2  = [P(0)] + [P(1)]
P4  = [P(0)] + [int_to_poly(0xbc)^i for i in range(1,  4)]
P16 = [P(0)] + [int_to_poly(0x0d)^i for i in range(1, 16)]

# Assert that all the elements of B are not in a subfield of P
def assert_not_in_any_subfield(B):
    for b in B:
        assert (b not in [poly_to_int(k) for k in P2])
        assert (b not in [poly_to_int(k) for k in P4])
        assert (b not in [poly_to_int(k) for k in P16])

# h is a integer
def get_matrice_multiplication(h):
    M = []
    ph = int_to_poly(h)

    for i in range(7, -1, -1):
        p = int_to_poly(1 << i)
        r = poly_to_int(p * ph)
        M.append([GF(2)((r >> k) & 1) for k in range(7, -1, -1)])
    return matrix(M)


def convertResult(B):
    res = []
    for b, coeff in B:
        p = get_matrice_multiplication(b).characteristic_polynomial().coefficients(sparse=False)
        v = sum([1 << power for power, value in enumerate(p) if value])
        res.append((v, b, coeff))
    return res

def generateCombinaison(isEncrypt=True):

    if isEncrypt:
        pos = ((2, 3, 1, 1), (1, 2, 3, 1), (1, 1, 2, 3), (3, 1, 1, 2))
    else:
        pos = ((14, 11, 13, 9), (9, 14, 11, 13), (13, 9, 14, 11), (11, 13, 9, 14))

    comb = []
    for i0 in range(4):
        for i1 in range(4):
            if i0 == i1:
                continue
            d = sorted([(pos[0][i0], pos[0][i1]), (pos[1][i0], pos[1][i1]),
                        (pos[2][i0], pos[2][i1]), (pos[3][i0], pos[3][i1])])
            if d not in comb:
                comb.append(d)
    comb2 = []
    for c in comb:
        for b1 in range(3):
            b0 = 3
            b2 = ((b1 + 1) % 3)
            b3 = ((b1 + 2) % 3)
            r = sorted([sorted([c[b0], c[b1]]), sorted([c[b2], c[b3]])])

            if r not in comb2:
                comb2.append(r)
    return comb2

def generateCoeffList(isEncrypt):

    existingCars = {}
    for comb in generateCombinaison(isEncrypt):
        currentCoeff = []
        for c in comb:
            ((a00, a01), (a10, a11)) = c
            c0 = sorted([a00, a11])
            c1 = sorted([a10, a01])
            a0 = int_to_poly(a00) * int_to_poly(a11)
            a1 = int_to_poly(a01) * int_to_poly(a10)
            currentCoeff.append(( poly_to_int(a0 / a1), (c0, c1)))
            currentCoeff.append(( poly_to_int(a1 / a0), (c1, c0)))
        assert_not_in_any_subfield([x for x, _ in currentCoeff])
        currentCoeff = sorted(convertResult(currentCoeff))

        req = tuple([c[0] for c in currentCoeff])
        res = tuple([c[1] for c in currentCoeff])
        coeff = tuple([c[2] for c in currentCoeff])
        if req in existingCars:
            assert existingCars[req] == (res, coeff), f"{req}, {res}, {coeff}"
        else:
            existingCars[req] = (res, coeff)
            print(f"{req}, {res}, {coeff}")

generateCoeffList(True)
generateCoeffList(False)
