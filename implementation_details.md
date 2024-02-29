# BGE implementation

This file intends to explain some of the choices that diverge from the
original [[BGE](#BGE)] paper and explain some of our design choices.

We try to do our best to write the note in a unified way, by trying to
keep the notation of [[BGE](#BGE)]. The basis vectors of the matrix are
row-wise oriented, so multiplying a matrix `A` by a vector `v` will be denoted
`v * A`.

## Setup

Suppose that we have access to a white-box, where each round `r` and
each column of a state `c` can be viewed as

```txt
                    ---------------------------------------------
  8 bits x[r, 3, c] | P[r, 3, c] | T[r, 3, c] |    | Q[r, 3, c] | y[r, 3, c] 8 bits
                    |--------------------------    -------------|
  8 bits x[r, 2, c] | P[r, 2, c] | T[r, 2, c] |    | Q[r, 2, c] | y[r, 2, c] 8 bits
                    |-------------------------- MC -------------|
  8 bits x[r, 1, c] | P[r, 1, c] | T[r, 1, c] |    | Q[r, 1, c] | y[r, 1, c] 8 bits
                    |--------------------------    -------------|
  8 bits x[r, 0, c] | P[r, 0, c] | T[r, 0, c] |    | Q[r, 0, c] | y[r, 0, c] 8 bits
                    ---------------------------------------------
```

where `T[round, i, c](v)` is `SubBytes(v) + RoundKey[idx, c]`

where `P^(-1)[r + 1, i, c] = Q[r, i, c]`. Since we will work always on
the same round, we will remove the index of the round in the notation if
it is not necessary. We also work on the same column `c`, so we will
also remove this component.

## Recovering non-linear parts, Section 3.1

**Note: in this section, `(S, o)` denotes the set `S` with the operation
`o`, which is the composition of functions of `S`.**

Referring to the diagram above, the output `y[0]` is a function of
`(x[0], x[1], x[2], x[3])`. Let `c1, c2` and `c3` be three constants,
the function `y[0]` can be written as `y[0](x[0], c1, c2, c3) = Q[0](a *
T[0](P[0](x[0]) + β[c1, c2, c3]))`. Let `c1'` be another constant, we
can write `y[0](x[0], c1', c2, c3)^(-1) = P[0]^(-1)(T[0]^(-1)(a^(-1) *
(Q[0]^(-1)(x) + β[c1', c2, c3])))`.  We can therefore write that
`y[0](y[0]^(-1)(x[0], c1', c2, c3), c1, c2, c3) =
Q[0](Q[0]^(-1)[0](x[0]) + β)`, where `β = β[c1, c2, c3] + β[c1', c2,
c3]`.

Let `S` be the set of functions of the form `Q(Q^(-1)(x) + β)`. This set
contains 256 functions, from `GF(2^8)` to `GF(2^8)`. The goals of the
first step of the [[BGE](#BGE)] attack is to recover `Q` up to an affine
mapping `A`, that is finding `Qtilde(x) = Q(A(x))`. In order to recover
`Qtilde`, the first step is to recover a group isomorphism `ψ` from `(S,
o)` to `(GF(2)^8, +)`. By following [[BGE](#BGE)], we recover `Qtilde`
by applying on each function `f` of `S` the equality: `f(0) = Qtilde(ψ(f))
`.

### Recover `ψ`

The function `ψ` maps `(S, o)` to `(GF(2)^8, +)`. Let `φ` be the map
from `(S, o)` to `(GF(2)^8, +)` where the binary representation of `β`
is viewed as a vector of `GF(2)^8` in the canonical basis of the vector
space `GF(2)^8`. However, since `β`is unknown, we may not find `φ`, but
another maps `ψ` that maps functions of `S` into vector of `GF(2)^8` in
an arbitrary basis. A change from a basis to another is a linear map,
hence the fact that `Qtilde` is equal to `Q` up to a linear map.

An algorithm to find `ψ` is the following:

0. Let `i` be equal to `1`.
1. Pick a random function `f` from `S` and remove `f` from `S`.
2. If `(f, ·)` is not in the list `R`, include in `R` the tuple `(f,
   ei)`, where `ei` is the `i`-th basis vector (i.e. `ei = (0, 0, …,0,
   1, 0, 0, …, 0)` where the `1` is at the `i`-th place). Set `ψ(f) =
   e`. For all the tuples `(g, η)` in `R`, set `ψ(f(g)) = ψ(f) + ψ(g) =
   e + η`.
3. If `S` is not empty, increment `i` and go to step `1`.

**Note**: this presentation is a bit different from the one of
[[BGE](#BGE)], since the constructed basis is built by repeatedly
multiplying by `02` the element `ei` initialized to `01`, which gives
the sequence `(2, 4, 8, 16, 32, 64, 128, 27)`.

#### Classical way

**Note**: you can find the initial description in [[BGE](#BGE)]. What we
present here mimic the algorithm by incorporating our modification with
respect to the basis vector.

Let us describe first how we think that object must be managed and
created.

##### Objects

Let `f` be a function of `S`, described by its 256 values in `GF(2^8)`
(f.val).  A unique identifier of a function `f` is its evaluation in `0`
(f.id = f.val[0]) for example, as in[[Tol](#Tol)].

The set `S` can be considered as a stack where the functions are stored
(only the pop() function will be used). Since at a point of the
algorithm, we need to refer to the functions stored in `S`, we consider
`S` as a hash map which maps function identifier to the function `f` in
the chosen basis and `s` as a list of the 256 values in [0, 256[, not
necessarily sorted.

The set `R` can be considered as a hash map that maps the function
identifier to an 8-bit integer that represents `β` in a basis of
`GF(2)^8`.

Finally, the function `ψ`, which maps `(S, o)` to `(GF(2)^8, +)`, can be
represented by a hash map that maps function identifier to decomposition
of `β` in the chosen basis.

Note that, since the hash function of the hash map is the function that
maps a function `f` of `S` to its function identifier, hash map are
simple array. We denote by `R.nb` the number of elements inserted into
`R`, by `R.contains(id)` the request to know if something was defined at
the index `id` of `R`, by `R.keys()` the list of all the keys (i.e.
indexes) defined in `R` (remember that the keys is the function
identifier) and by `R.find(id)` the functions that returns the element
at the index `id`.

##### Algorithm

We now try to write the algorithm described in [[BGE](#BGE)] in
pseudo-code, with our modification.

```txt
R[0] <- 0
ψ[0] <- 0
e <- 1
i = 0
while R.nb < 2^8:
  f <- S.find(s[i])
  i <- i + 1
  if not R.contains(f.id):
    ψ[f.id] <- e
    for g in R.keys():
      fogid <- f.val[g.id] // remember that g.id = g.val[0], then f.val[g.id] = f.val[g.val[0]]
      vec <- e + R[g.id] // + equal xor in this context
      R[fogid] <- vec
      ψ[fogid] <- vec
  e <- e << 1
```

As we can see, `R` and `ψ` play the same role when it is necessary to
store the information about the value of a function in a basis. In order
to remove `R`, we need to store the information about which function is
in `R` or not and the number of elements in `R`. We replace `R` by
respectively an array `inR` which stores at index `i` if the function of
identifier `i` would be in `R` or not, and by an integer `nbR` the
number of elements in `R`. The array `s` is also not necessary, we can
take the function in the order of how they are stored in `S` (in order
to mimic [[Tol](#Tol)], we store at index `i` of `S` the function of
identifier `i`). To cover all the elements of `S`, it suffices to loop on
all the possible integers in [0, 2^8[, then `nbR` is not useful.  Still
to mimic [[Tol](#Tol)], we change a bit the way to construct `e`.

We then have

```txt
ψ[0] <- 0
inR[0] <- true

for i in [1, 2^8[:
  ψ[i] <- 0
  inR[i] <- false

k <- 0
for i in [0, 2^8[:
  if not inR[i]:
    k <- k + 1
    ψ[i] <- 1 << (k - 1)
    inR[i] <- true
    for j in [0, 2^8[:
      if inR[j]:
        fogid <- S[i].val[j] // remember that g.id = g.val[0], then f.val[g.id] = f.val[g.val[0]]
        ψ[fogid] <- ψ[i] + ψ[j] // + equal xor in this context
        inR[fogid] <- true
```

Since we know that we will run through all the function, we can change
the "`for i in [0, 2^8[`" into a while loop on `k` and with a few minor
changes, we get the [[Tol](#Tol)] algorithm.

#### Algorithm in [[Tol](#Tol)]

We just rewrite the algorithm

```
ψ[0] <- 0
inR[0] <- true

for i in [1, 2^8[:
  ψ[i] <- 0
  inR[i] <- false

k <- 0 // in the article, k = j
while k != 8:
  i <- 1
  while inR[i]:
    i <- i + 1

  k <- k + 1
  ψ[i] = 1 << (k - 1)
  inR[i] <- true
  for j in [1, 2^8[:
    if inR[j]:
      fogid <- S[i].val[j] // in the article, fogid = m, but we keep the previous notation
      ψ[fogid] <- ψ[i] + ψ[j]
      inR[fogid] <- true
```

Note that the index `i` in the while loop seems to be initialized only
one time to `1` outside of the loop. This gives

```
ψ[0] <- 0
inR[0] <- true

for i in [1, 2^8[:
  ψ[i] <- 0
  inR[i] <- false

k <- 0 // in the article, k = j
i <- 1
while k != 8:
  while inR[i]:
    i <- i + 1

  k <- k + 1
  ψ[i] = 1 << (k - 1)
  inR[i] <- true
  for j in [1, 2^8[:
    if inR[j]:
      fogid <- S[i].val[j] // in the article, fogid = m, but we keep the previous notation
      ψ[fogid] <- ψ[i] + ψ[j]
      inR[fogid] <- true
```

This is what is implemented in
[`approximateencoding.hpp`](src/approximateencoding.hpp), computing the
inverse of `ψ` at the same time.

## Relation between affine parasites, Section 3.2

We choose to not use the linear dependencies that may be found between the
affine encodings and directly use the result of Section 3.3 on each independent
encodings. We, however, need to implement part of Section 3.2 for a result in
Section 3.3.

#### Find `L[i,j,x]` and `c[i,j,x]`

Remember that `y[i)` and `y[j]` are 8-bit length. Then, the matrix
`M(L[i,j,0])` which represents `L[i,j,0]` has 64 entries in `GF(2)` and
`c[i,j,x]` 8 entries. We know `y[i]` and `y[j]` by their values.

Implementing an algorithm to solve a tiny linear system on GF(2) is not
the only way to get `M(L[i,j,0])` and `c[i,j,0]`. This requires an
inversion of an 8×8 matrix over GF(2) (which needs around `O(n^3)` with
classical algorithm) constructed by requiring at least 8 times `y[j]` to
get an invertible matrix (see
[wiki](https://en.wikipedia.org/wiki/General_linear_group#Over_finite_fields)
to estimate probability), and a bit more to get the non-linear part. We
propose an algorithm that works in almost `O(2^n)`, but when `n=8`, we
are just below the crossing point.

We will use the inverse of the function in order to recover the unknown
parts.  Let us begin with `c[i,j,0]`: we have `y[i](y[j]^(-1)(0), 0, 0,
0) = c[i,j,0]`. Let us consider in the following that the functions `yi`
and `yj` output a bit-vector.  Then, we can write `y[i](x[0], 0, 0, 0) =
y[j](x[0], 0, 0, 0) * M(L[i,j,0]) + c[i,j,0]`. Like the computation of
`c[i,j,0]`, we can choose specific values of `y[j](x[0], 0, 0, 0)`.
Choosing values such that `y[j](x, 0, 0, 0) = 2^k`, where `k` is a power
of 2 less than `2^n`, allows to find the row coefficients of a line. For
example, to find the last line of `M(L[i,j,0])`, it suffices to find `x`
such that `y[j](x, 0, 0, 0) = 0b00…01`, that is `x = y[j]^(-1)(1)`. The
previous line is found by setting `x = y[j]^(-1)(2)`. We then recover
`M(L[i,j,0])` and `c[i,j,0]` only with evaluation of functions, instead
of solving linear system. This is what is implemented in
[`relationencoding.hpp`](src/relationencoding.hpp).

## Recovering the affine parasites, Section 3.3

### About `Q[0]`'s linear part

#### Computation of the set B

We adapt a bit the list B of the possible β, as computed in the `find_B()`
function in the following [Sage](https://www.sagemath.org/) file. Indeed, the
`α[i,j]` come from the `MixColumn` coefficients, `β` takes its values in `B
= {'7b', 'a4', '8d', '8c', '46', 'f5', '06', '02', 'f6', '8f', 'f7', '03'}`. In
the article, `B` contains more values since the authors consider the set of
`(α[i,x] * α[j,y]) / (α[i,y] * α[j,x])`, where `x ≠ y` and `(i, j, x, y)` in
`[0, 4)^4`. In the set we describe, we only consider `(α[i,x] * α[j,y]) /
(α[i,y] * α[j,x])`, where `x ≠ y` and `(i, j, x, y)` in `[0, 4)^2 ×
[0, 2)^2`. Note that in this set, there is no square, on the contrary to the set
described in the original article.

```py
Q.<y> = GF(2)[]
P.<x> = GF(2^8, name='x', modulus=y^8+y^4+y^3+y+1)
assert((x^6+x^4+x^2+x+1) * (x^7+x+1) == x^7+x^6+1)

# p is a polynomial of P
def poly_to_bin(p):
    L = p.polynomial().coefficients(sparse=False)
    L.reverse()
    L = [0 for i in range(8 - len(L))] + L
    return "".join([str(i) for i in L])

# p is a polynomial of P
def poly_to_hex(p):
    L = hex(Integer(poly_to_bin(p), 2)).split('x')[1]
    return "".join(["0" for i in range(2- len(L))]) + L

# s is a string with just the binary representation
def bin_to_poly(s):
    k = 7
    p = P(0)
    for i in s:
        p = p + GF(2)(i) * x^k
        k = k - 1
    return p

# s is a string with just the hexadecimal representation
def hex_to_poly(s):
    L = bin(Integer(s, 16)).split("b")[1]
    L = "".join(["0" for i in range(8 - len(L))]) + L
    return bin_to_poly(L)

P2  = [P(0)] + [P(1)]
P4  = [P(0)] + [hex_to_poly("bc")^i for i in range(1,  4)]
P16 = [P(0)] + [hex_to_poly("0d")^i for i in range(1, 16)]

# Assert that all the elements of B are not in a subfield of P
def assert_not_in_any_subfield(B):
    for b in B:
        assert((b not in [poly_to_hex(k) for k in P2]) and (b not in [poly_to_hex(k) for k in P4]) and (b not in [poly_to_hex(k) for k in P16]))

# B as it is written in the article (do not respect the FIPS notation)
B_original = ['02', 'd8', '03', '6f', '04', 'bc', '06', 'b7', '05', '25', '4a', 'f8', '7f', 'c8', '64', '5f']
B = []
for i in B_original:
    if i[0] == '0':
        B.append(i)
    else:
        B.append(i[1] + i[0])

assert_not_in_any_subfield(B)

# Recompute B by our own way
def find_B():
    B = []
    L = [[('02', '03'), ('01', '02'), ('01', '01'), ('03', '01')], # should be sufficient to compute L on 3.3
         # [('02', '01'), ('01', '03'), ('01', '02'), ('03', '01')],
         # [('02', '01'), ('01', '01'), ('01', '03'), ('03', '02')],
         # [('03', '01'), ('02', '03'), ('01', '02'), ('01', '01')],
         # [('03', '01'), ('02', '01'), ('01', '03'), ('01', '02')],
         # [('01', '01'), ('03', '01'), ('02', '03'), ('01', '02')]
         ]

    for l in L:
        for (a00, a01) in l:
            for (a10, a11) in l:
                if (a00 == a10 and a01 == a11):
                    continue
                B.append(poly_to_hex(hex_to_poly(a00) * hex_to_poly(a11) / (hex_to_poly(a01) * hex_to_poly(a10))))
    return list(set(B))


B_sufficient = find_B()
for b in B_sufficient:
    assert(b in B)
```

From the first phase, we know that we can assume that all encodings, the `Q[i]`
and the `P[i]` are affine, especially, we will denote `Q[i](x) = A[i](x) + q[i]`. From
the second phase, we have the ability to get `M(Lijx)` and `cijx`. Let consider
that we want all the linear relation involving `i = 0`, we are especially able
to compute `M(L[0,j,x])`. The goal here is to compute `Atilde[0](x) = A[0](Λ[𝛾](x))`.

Let now consider `L[i,j,0,1] = L[i,j,0](L[i,j,1]^(-1)(x))`, then `M(L[i,j,0,1])
= M(A[0])^(-1) * M(Λ[β]) * M(A[0])`. Since `M(L[i,j,0,1])` and `M(Λ[β])` share
the same characteristic polynomial: a way to identify `β` is to compute the
characteristic polynomial of `M(L[i,j,0,1])`.

```py
# h is a string with just the hexadecimal representation
def get_matrice_multiplication(h):
    M = []
    ph = hex_to_poly(h)

    s = "10000000"
    for i in range(0, 8):
        p = bin_to_poly(s)
        r = p * ph
        M.append([GF(2)(k) for k in poly_to_bin(r)])
        s = "0" + s[0:7]
    return matrix(M)

for b in find_B():
    "".join([str(i) for i in get_matrice_multiplication(b).characteristic_polynomial().coefficients(sparse=False)[::-1]])
```

Given `M(L[i,j,0,1])`, we look for `Atilde[0](x) = A[0](Λ[𝛾](x))`, where
`𝛾 ≠ 0`. A way to find this function is to look for
`L[i,j,x,y](Atilde[0](x)) = Atilde[0](Λ[β](x))`. Since we now the
matrices of `L[i,j,x,y]`, and then `β` and therefore `Λ[β]`, we need to
solve `M(Atilde[0]) * M(L[i,j,x,y]) = M(Λ[β]) * M(Atilde[0])`. To solve that, a way is to develop each product keeping
the elements of `M(Atilde[0])` unknown, giving a system of 64 equations.
On an example of 3×3 matrices instead of 8×8, it gives

```txt
(s[0,0], s[0,1], s[0,2], s[1,0], s[1,1], s[1,2], s[2,0], s[2,1], s[2,2]) * matrix(
  [[L[0,0] + Λ[0,0], L[0,1]      , L[0,2]      , Λ[1,0]       , 0       , 0        , Λ[2,0]      , 0        , 0],
   [L[1,0]      , L[1,1] + Λ[0,0], L[1,2]      , 0         , Λ[1,0]     , 0        , 0        , Λ[2,0]      , 0],
   [L[2,0]      , L[2,1]      , L[2,2] + Λ[0,0], 0         , 0       , Λ[1,0]      , 0        , 0        , Λ[2,0]],
   [Λ[0,1]      , 0        , 0        , L[0,0] + Λ[1,1], L[0,1]      , L[0,2]      , Λ[2,1]      , 0        , 0],
   [0        , Λ[0,1]      , 0        , L[1,0]      , L[1,1] + Λ[1,1], L[1,2]      , 0        , Λ[2,1]      , 0],
   [0        , 0        , Λ[0,1]      , L[2,0]      , L[2,1]      , L[2,2] + Λ[1,1], 0        , 0        , Λ[2,1]],
   [Λ[0,2]      , 0        , 0        , Λ[1,2]      , 0        , 0        , L[0,0] + Λ[2,2], L[0,1]      , L[0,2]],
   [0        , Λ[0,2]      , 0        , 0        , Λ[1,2]      , 0        , L[1,0]      , L[1,1] + Λ[2,2], L[1,2]],
   [0        , 0        , Λ[0,2]      , 0        , 0        , Λ[1,2]      , L[2,0]      , L[2,1]      , L[2,2] + Λ[2,2]]])
 = (s[0,0], s[0,1], s[0,2], s[1,0], s[1,1], s[1,2], s[2,0], s[2,1], s[2,2]) * S = 0.
```

We then look for the kernel of `S` and recover the matrix `Atilde[0] =
matrix([[s[0,0], s[0,1], s[0,2]], [s[1,0], s[1,1], s[1,2]], [s[2,0], s[2,1],
s[2,2]]])`.

## Support of Shuffled States and of Decryption

The [first version](https://github.com/SideChannelMarvels/BlueGalaxyEnergy/releases/tag/v1.0.1)
of BlueGalaxyEnergy tools implements BGE for encryption whiteboxes with unshuffled
internal states. For the next version, we improved the tools by handling
shuffled states and decryption whiteboxes, with minimal changes to the
previous implementation.

### Shuffled States

In order to support shuffled states, we needed to add three features :

- Before the inputs generation, finding a permutation that imitates the
  byte propagation of an AES round.
- During recovery of affine parasites (Section 3.3), manipulating the
  coefficient of the MixColumn. We need to find the impact of each input byte of a
  column on each output byte. We also need to support `β` for every characteristic
  polynomial of this phase, as we skip the support of some `β` that will never
  appear in an unshuffled state.
- After the key bytes recovery, using the impact coefficient found
  during the affine parasites recovery to finish unshuffling the
  state. Even with these coefficients, we still have 16 possibilities: we
  use the key schedule to select the correct key.

The two last features followed the Phase 4 of [[DMRP](#DMRP)] with local
improvement.

#### First permutation

When generating the inputs for the BGE attack, we use optimized inputs (like
`(x,0,0,0,x,0,0,0,x,0,0,0,x,0,0,0)`) to lower the number of inputs needed.
These optimized inputs allow us to reduce by 4 the number of calls to the whitebox.
Indeed, for a given round, changing one input byte will only impact the 4
output bytes of the same column. By using inputs that target the four columns at
once, we reduce the load by 4. However, these optimized inputs can only be used
if each input byte has the expected impact on the output byte. Otherwise, an
optimized input may fault the same columns twice, while another column will be
unchanged.

We need to find a permutation that emulates the impact of each input byte
on the expected output bytes.

| Input bytes (encrypt) | Input bytes (decrypt)  | Impacted Output Bytes |
|-----------------------|------------------------|-----------------------|
| 0, 5, 10, 15          | 0, 7, 10, 13           | 0, 1, 2, 3            |
| 3, 4, 9, 14           | 1, 4, 11, 14           | 4, 5, 6, 7            |
| 2, 7, 8, 13           | 2, 5, 8, 15            | 8, 9, 10, 11          |
| 1, 6, 11, 12          | 3, 6, 9, 12            | 12, 13, 14, 15        |

For each round, we need a permutation on both the inputs and the outputs.
However, for consecutive rounds, the permutation between them must be the
same (both as the output of the previous round and as the input of the
next round). Once we found a permutation for each state, BGE input generation
is handled on a modified whitebox that applies the permutation in each round.

The permutations chosen at this step are kept and will be updated when the
MixColumn coefficients are known for each round.

#### Support of all β and MixColumn coefficient extraction during affine parasites recovery

Although columns remain mixed after the previous step,
implementing all possible `β` values associated with their
corresponding characteristic polynomials of `L` is necessary.
However, during this listing process, we discovered that
two distinct `β` values share same characteristic polynomial.

However, most collisions can be avoided by grouping multiple characteristic
polynomials with shared characteristics. Consider two input
bytes `(i, j)` in `[0, 4)^2` with `i ≠ j` and four output bytes
`(u, v, x, y)` in `[0, 4)^4` that are all distinct.
With this setup, we can compute the affine parasite relationships for
four different `β` values:

- `(α[i,u] * α[j,v]) / (α[i,v] * α[j,u])`,
- `(α[i,v] * α[j,u]) / (α[i,u] * α[j,v])`,
- `(α[i,x] * α[j,y]) / (α[i,y] * α[j,x])`,
- `(α[i,y] * α[j,x]) / (α[i,x] * α[j,y])`.

Grouping `β` in sets of four resolved most collisions on characteristic polynomials.
Only one collision remained with polynomial 471, encountered
in shuffled decryption whiteboxes (where `β` can be 54 or 102). Here,
we iterate on the MixColumn coefficient recovery algorithm (explained later),
which succeeds only with correct coefficients.

In addition to recovering the correct `β` for each characteristic polynomial of `L`,
we also need to recover every `α[i,x]` for every `(i, x)` in `[0, 4)^2`.
We known that `α[i,x]` takes the value of the MixColumns and, for encryption:

- `{α[i,0], α[i,1], α[i,2], α[i,3]} == {1, 1, 2, 3}`
- `{α[0,x], α[1,x], α[2,x], α[3,x]} == {1, 1, 2, 3}`

(for decryption, use the decryption MixColumns coefficients, which are `{9, 11, 13, 14}`).
When we recover a `β` value, we can limit the possible coefficients associated with each
`α[i,x]`. If we find the value
`β0 = (α0 * α1) / (α2 * α3) = (α[i,x] * α[j,y]) / (α[i,y] * α[j,x])`, then we know that
`α[i,x]` is either `α0` or `α1`. to determine the correct value, we resolve
the same system a second time with the following relationship of affine parasites:

- `(α[i,u] * α[j,y]) / (α[i,y] * α[j,u])`,
- `(α[i,y] * α[j,u]) / (α[i,u] * α[j,y])`,
- `(α[i,x] * α[j,v]) / (α[i,v] * α[j,x])`,
- `(α[i,v] * α[j,x]) / (α[i,x] * α[j,v])`.

With this system, we can find
`β1 = (α0 * α4) / (α5 * α3) = (α[i,x] * α[j,v]) / (α[i,v] * α[j,x])`.
This allows us to select `α[i,x]` that exists in both
the set `{α0, α1}` and the set `{α0, α4}`. However, the
coefficients of the encryption MixColumn have twice the coefficient 1.
This means it's possible to have `α1` and `a4` both equal to 1,
while the expected `α[i,x] == α0` has
another value. As this is the only corner case, if the intersection of
`{α0, α1}` and `{α0, α4}` returns two values, and one of them is 1, we select
the other value.

By solving this equation, we can retrieve all the MixColumn
coefficients for the input bytes `(i, j)`. We then repeat the same algorithm
to find the coefficients associated with the two other input bytes.

The script [precomputeBetaTable.sage](src/precomputeBetaTable.sage) precomputes
all possible characteristic polynomials, along with their associated `β`
values and coefficients.

#### Final unshuffling of the internal state

For each round and its columns, the previous steps allowed us to retrieve the MixColumn coefficients.
We also know which bytes interact with each column.
By associating these coefficients with the permutation chosen in the first step,
we can determine the MixColumn coefficient for each input byte impacting each output byte.

Out of the 576 possible permutations (`4! * 4!`) for a round's column,
only 4 will have coefficients positioned correctly to resemble a MixColumn.
For each round's column, if we pick a random byte and place it first, only one permutation remains.

After eliminating all impossible permutations for each round's column, we still
have 6144 possible permutations for a given round (`4^4 * 4!`).  This is because
we still don't know the relative order of bytes within each column, and the
relative order of columns within the round.

At this step, [[DMRP](#DMRP)] proposes iterating through the 6144 possibilities.
However, we can eliminate more incompatible permutations by using the next round, reducing the number
of possibilities to 16 without yet using the keyschedule.

To achieve this, we choose one input byte from round N and place it first.
Its column is considered the first of the inputs, and the byte itself is the first within the column.
This constraint fixes four bytes in the round's outputs.

With these four output bytes fixed, we move to round N+1.
Due to the ShiftRow operation, each fixed byte will occupy a different column in round N+1,
further fixing the order of both columns and bytes within each column.
This completely defines the permutation between N and N+1,
as well as the output permutation for N+1. We can repeat this process for subsequent rounds.

For the input permutation of round N, we can perform the reverse process with all output bytes now fixed.
By the end of this algorithm, only 16 permutations remain, each with a distinct byte at the first position of the first round input.

To find the correct permutation, we iterate through each one, extract the associated roundkeys,
and verify if they correspond to the actual AES key. For AES-128 and AES-256,
this increases the minimum required number of rounds by 1.

### Support of decryption

To support decryption whiteboxes, we need to rearrange the AES steps.

During encrypt, the round operations order is as follows:

- Xor the round key
- Perform the SBox
- Perform the ShiftRow
- Perform the MixColumn

The last round differs with:

- Xor the round key
- Perform the SBox
- Perform the ShiftRow
- Xor the last round key

To achieve a similar structure for decryption, we consider the following AES order.

First round:

- Xor the last round key
- Perform the InvShiftRow
- Perform the InvSBox
- Xor the round key
- Perform the InvMixColumn

Intermediary rounds:

- Perform the InvShiftRow
- Perform the InvSBox
- Xor the round key
- Perform the InvMixColumn

Last round:

- Perform the InvShiftRow
- Perform the InvSBox
- Xor the round key

Due to internal properties of AES steps, we can rearrange the rounds as :

First round:
- Xor the last round key
- Perform the InvSBox
- Perform the InvShiftRow
- Perform the InvMixColumn

Intermediary rounds:

- Xor the InvMixColumn of the previous round key
- Perform the InvSBox
- Perform the InvShiftRow
- Perform the InvMixColumn

Last round:

- Xor the InvMixColumn of the previous round key
- Perform the InvSBox
- Perform the InvShiftRow
- Xor the round key

With this final representation, we get a similar structure for encryption and
decryption, with each operation replaced by its inverse. The main difference lies
in the requirement to perform an InvMixColumn on most round keys during decryption.
Consequently, we need to perform a MixColumn on all keys recovered through the BGE
attack before feeding them into the key schedule to find the decryption key.

While most steps of the BGE attack can be readily adapted for decryption, we
discovered that Proposition 3 of [[BGE](#BGE)] no longer holds true when the
SBox is replaced with it inverse.

Proposition 3 states that:

```txt
There exist unique pairs (δi, ci) i=0,...,3 of elements in GF(2^8), δi being non-zero, such that

Ptilde0 : x -> (InvSbox ◦ Λδ0 ◦ Atilde0^{-1} ) ( y0(x, ‘00’, ‘00’, ‘00’) ⊕ c0)
Ptilde1 : x -> (InvSbox ◦ Λδ1 ◦ Atilde0^{-1} ) ( y0(‘00’, x, ‘00’, ‘00’) ⊕ c1)
Ptilde2 : x -> (InvSbox ◦ Λδ2 ◦ Atilde0^{-1} ) ( y0(‘00’, ‘00’, x, ‘00’) ⊕ c2)
Ptilde3 : x -> (InvSbox ◦ Λδ3 ◦ Atilde0^{-1} ) ( y0(‘00’, ‘00’, ‘00’, x) ⊕ c3)

are affine mappings. Any pair (δi, ci) can be computed with time complexity 2^{24}.
Moreover, those mappings are exactly Ptildei = Pi(x) ⊕ ki.
```

Replacing InvSbox with SBox during the decryption whitebox still yields unique
`ci i=0,...,3` values that create affine mappings for each Ptildei.
However, with these `ci` values, any non-zero `δi` value also produces an affine mapping.

We hypothesize that this discrepancy arises from the SBox's internal structure.
Notably, the SBox can be represented as a composition of two functions:
an affine mapping (`aff_1f_63`)
and the multiplicative inverse operation within the finite field GF(2^8) (MultInv).
Symbolically, this is expressed as `SBox = aff_1f_63 ◦ MultInv` and
its inverse as `InvSbox = MultInv ◦ aff_1f_63^{-1}`.

Proposition 3 is valid for encryption whiteboxes because the affine mapping `aff_1f_63`
is positioned between `MultInv` and `Λδi`. However, in the decryption setting,
`aff_1f_63` no longer occupies the same relative position.

When the correct `ci` values are used, all operations between the `SBox`
and its inverse within the `y0` function become multiplications within GF(2^8).
Since multiplication in this field is commutative, the two `MultInv` operations
(one inside `y0` and the other within `Sbox`) effectively cancel each other out.

Consequently, the entire composition lacks any multiplicative inverse operations,
rendering all resulting compositions affine mappings.
Unfortunately, this absence of MultInv eliminates the capability to uniquely identify the correct `δi` values.

While we cannot determine the exact `δi` values due to the limitations mentioned earlier,
we still have 255 possible candidates.

The MixColumn coefficients provide us with the following relationship for a given output byte `x`:
`δ0 * α[0,x] = δ1 * α[1,x] = δ2 * α[2,x] = δ3 * α[3,x]`
Furthermore, fixing `δ0` to an arbitrary value for a specific output byte `x`
uniquely defines the corresponding `Ptildei`, which are also shared with the other output bytes in the same column.
Therefore, for a given column, we only have 255 possibilities,
which can be enumerated by iterating through all possible values of `δ0` for the first output byte.

To identify the correct `δ0` value, we leverage the key extraction equation.
In encryption whiteboxes, the key can be extracted by composing a `Ptilde` mapping
from round N+1 with the corresponding `Q` mapping from round N.
This composition should result in an affine function with a linear coefficient of 1.

The composition depends on both the `δ0` value of the `Ptilde` mapping for a specific column in round N+1
and the `δ0` value of the `Q` mapping for a column in round N.
Since we lack the exact values, we need to analyze all possible combinations.
With 255 possibilities for each `δ0`, this amounts to `255 * 255 ≈ 2^{16}` possible combinations.

Through this analysis, we identified only one pair of `δ0` values that yields
an affine mapping with a linear coefficient of 1.
This unique pair represents the correct `δ0` values for both the `Ptilde` and `Q` mappings.

The brute-force approach requiring `2^{16}` iterations only needs to be performed once for a given decryption whitebox.
Once completed, we know the `δ0` values for a specific column in round N and a column in round N+1.
Due to the ShiftRow operation, each `Q` mapping from round N's column will be paired with a `Ptilde` from a different column in round N+1.

Therefore, we only need to repeat the 255-iteration brute force for each remaining pair of columns between rounds N and N+1.
Notably, by iteratively processing subsequent rounds (N+2, N+3, ...) in order, we'll always have a fixed `δ0` value for the `Q` mapping.

Based on this observation, we note that the overall complexity of the BGE attack on decryption whiteboxes is lower than that of encryption whiteboxes.

## Bibliography

<a name="BGE">[BGE]</a> Olivier Billet, Henri Gilbert and Charaf
Ech-Chatbi, Cryptanalysis of a White Box AES Implementation, SAC 2004.
[pdf](https://link.springer.com/content/pdf/10.1007%2F978-3-540-30564-4_16.pdf)

<a name="LRDMRP">[LRDMRP]</a> Tancrède Lepoint, Matthieu Rivain, Yoni De Mulder,
Peter Roelse and Bart Preneel, Two Attacks on a White-Box AES
Implementation, SAC 2013. [pdf](https://link.springer.com/content/pdf/10.1007%2F978-3-662-43414-7_14.pdf)

<a name="Tol">[Tol]</a> Ludo Tolhuizen, Improved cryptanalysis of an AES
implementation. 33rd WIC Symposium on Information Theory in the Benelux,
2012. [pdf](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.913.5807&rep=rep1&type=pdf)

<a name="">[DMRP]</a> Yoni De Mulder, Peter Roelse and Bart Preneel,
Revisiting the BGE Attack on a White-Box AES Implementation,
Cryptology ePrint Archive, 2013.
[pdf](http://eprint.iacr.org/2013/450.pdf)
