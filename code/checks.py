"""Checks for the cancelation norm on a free group and the drift of the simple random walk.

Letters of the free group F_r are encoded as 0, ..., 2r-1, the inverse of a being a ^ 1.

1. The dynamic program for L agrees with the definition (minimum over all noncrossing
   inverse matchings, enumerated directly) for every word of length at most 7 over F_2.
2. For every word of length at most 8 over F_2: L(w) = L(red(w)), L(w) = 0 iff w reduces
   to the empty word, L(w^{-1}) = L(w). For all pairs u, v of length at most 3:
   L(uv) <= L(u) + L(v) and L(u v u^{-1}) = L(v).
3. The tail bound of Lemma 3.1 for r = 32 and c = 1/4, summed exactly, against 2^{-n/2},
   and the exponent of Lemma 3.2 for r = 3, 4, 8, 32.
4. Monte Carlo estimates of E L(X_n)/n for r = 2, 3, 32 (not used in the proofs).
Output: checks_out.txt (stdout) and checks.json.
"""
import itertools, json, math, pathlib, random, time
from fractions import Fraction
from functools import lru_cache

out = {}


@lru_cache(maxsize=200000)
def norm(w):
    """L(w) by the first-letter recursion: leave w[0] unpaired, or pair it with an inverse letter."""
    if not w:
        return 0
    v = 1 + norm(w[1:])
    for j in range(1, len(w)):
        if w[j] == w[0] ^ 1:
            v = min(v, norm(w[1:j]) + norm(w[j + 1:]))
    return v


def matchings(positions):
    """All noncrossing partial matchings of a sorted tuple of positions."""
    if not positions:
        yield ()
        return
    first, rest = positions[0], positions[1:]
    for m in matchings(rest):
        yield m
    for idx in range(len(rest)):
        inside, outside = rest[:idx], rest[idx + 1:]
        for a in matchings(inside):
            for b in matchings(outside):
                yield ((first, rest[idx]),) + a + b


def norm_by_definition(w):
    best = len(w)
    for m in matchings(tuple(range(len(w)))):
        if all(w[i] == w[j] ^ 1 for i, j in m):
            best = min(best, len(w) - 2 * len(m))
    return best


def red(w):
    s = []
    for a in w:
        if s and s[-1] == a ^ 1:
            s.pop()
        else:
            s.append(a)
    return tuple(s)


def inv(w):
    return tuple(a ^ 1 for a in w[::-1])


# 1. dynamic program against the definition
t0 = time.time()
count = 0
for n in range(8):
    for w in itertools.product(range(4), repeat=n):
        assert norm(w) == norm_by_definition(w), w
        count += 1
out["dp_vs_definition"] = dict(words=count, max_length=7, rank=2)
print(f"1. dynamic program equals the matching definition on all {count} words of length <= 7 in F_2")

# 2. metric properties on small words
count = 0
for n in range(9):
    for w in itertools.product(range(4), repeat=n):
        v = norm(w)
        assert v == norm(red(w)), w
        assert (v == 0) == (red(w) == ()), w
        assert v == norm(inv(w)), w
        count += 1
small = [w for n in range(4) for w in itertools.product(range(4), repeat=n)]
pairs = 0
for u in small:
    for v in small:
        assert norm(u + v) <= norm(u) + norm(v)
        assert norm(u + v + inv(u)) == norm(v)
        pairs += 1
out["metric_properties"] = dict(words=count, pairs=pairs)
print(f"2. reduction invariance, nondegeneracy, inverse symmetry on {count} words; "
      f"triangle inequality and conjugation invariance on {pairs} pairs")

# 3. the tail bounds
rows = []
for n in (8, 16, 32, 64, 128, 256):
    total = Fraction(0)
    for k in range(0, n // 4 + 1):
        if (n - k) % 2:
            continue
        m = (n - k) // 2
        catalan = math.comb(2 * m, m) // (m + 1)
        total += math.comb(n, k) * catalan * Fraction(1, 64 ** m)
    ratio = total / Fraction(1, 2 ** (n // 2)) if n % 2 == 0 else None
    rows.append(dict(n=n, bound=float(total), bound_over_2_pow_minus_half_n=float(ratio)))
    print(f"3a. r=32 n={n}: sum_k binom(n,k) C_m 64^-m over k <= n/4 = {float(total):.3e}, "
          f"ratio to 2^(-n/2) = {float(ratio):.3e}")
out["tail_r32"] = rows


def best_exponent(r, c):
    """min over s in (0,1) of -c log s + log(y + s), y = sqrt(2/r); negative means decay."""
    y = math.sqrt(2 / r)
    return min(-c * math.log(s) + math.log(y + s) for s in (i / 10000 for i in range(1, 10000)))


rows = []
for r in (3, 4, 8, 32):
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if best_exponent(r, mid) < 0:
            lo = mid
        else:
            hi = mid
    rows.append(dict(r=r, largest_c=lo))
    print(f"3b. r={r}: Lemma 3.2 gives exponential decay of P(L <= cn) for every c < {lo:.4f}")
out["chernoff"] = rows


# 4. Monte Carlo
def norm_interval(w):
    n = len(w)
    D = [[0] * (n + 1) for _ in range(n + 2)]
    pos = {}
    for i, a in enumerate(w):
        pos.setdefault(a, []).append(i)
    for i in range(n - 1, -1, -1):
        partners = [k for k in pos.get(w[i] ^ 1, []) if k > i]
        Di1 = D[i + 1]
        row = D[i]
        for j in range(i + 1, n + 1):
            v = 1 + Di1[j]
            for k in partners:
                if k >= j:
                    break
                x = Di1[k] + D[k + 1][j]
                if x < v:
                    v = x
            row[j] = v
    return D[0][n]


for w in itertools.product(range(4), repeat=6):
    assert norm_interval(list(w)) == norm(w)
rng = random.Random(20261105)
rows = []
for r in (2, 3, 32):
    for n, samples in ((40, 200), (80, 100), (160, 40)):
        vals = [norm_interval([rng.randrange(2 * r) for _ in range(n)]) / n for _ in range(samples)]
        mean = sum(vals) / samples
        rows.append(dict(r=r, n=n, samples=samples, mean_L_over_n=mean, min_L_over_n=min(vals)))
        print(f"4. r={r} n={n}: mean L(X_n)/n = {mean:.4f} over {samples} walks (smallest {min(vals):.4f})")
out["monte_carlo"] = rows
out["elapsed_seconds"] = round(time.time() - t0, 1)
pathlib.Path(__file__).with_name("checks.json").write_text(json.dumps(out, indent=1))
