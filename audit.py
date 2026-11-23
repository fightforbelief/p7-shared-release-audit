"""Sharp changes in a two-group rate gap under anchored site replacement."""
from dataclasses import dataclass
from itertools import combinations
from math import isclose, isfinite, sqrt


@dataclass(frozen=True)
class Group:
    share: tuple       # P(group | site)
    low: tuple         # feasible rate at each site
    high: tuple


@dataclass(frozen=True)
class Witness:
    value: float
    first: int
    second: int
    fraction_first: float


@dataclass(frozen=True)
class Envelope:
    lower: Witness
    upper: Witness
    single_lower: float
    single_upper: float


def dot(x, y):
    return sum(a * b for a, b in zip(x, y, strict=True))


def check(anchor, groups, rho):
    if not anchor or not all(isfinite(b) and b >= 0 for b in anchor):
        raise ValueError("anchor must contain finite nonnegative weights")
    if not isclose(sum(anchor), 1, abs_tol=1e-12, rel_tol=0) or not 0 <= rho <= 1:
        raise ValueError("invalid anchor sum or replacement fraction")
    if len(groups) != 2:
        raise ValueError("two groups required")
    for g in groups:
        if any(len(a) != len(anchor) for a in (g.share, g.low, g.high)):
            raise ValueError("inconsistent site dimensions")
        if not all(0 < q < 1 and 0 <= l <= h <= 1
                   for q, l, h in zip(g.share, g.low, g.high)):
            raise ValueError("invalid group proportions or rate bounds")
    if any(not isclose(a+b, 1, rel_tol=0, abs_tol=1e-12)
           for a, b in zip(groups[0].share, groups[1].share)):
        raise ValueError("the two group proportions must sum to one")


class OrientedGroup:
    """Maximize a signed contrast, with one default corner for absent sites."""
    def __init__(self, anchor, group, sign, rho):
        self.b = anchor
        self.q = group.share
        self.base_mass = dot(anchor, self.q)
        # Work with signed outcomes, so every component is maximized.
        if sign == 1:
            self.low, self.high = group.low, group.high
        else:
            self.low = tuple(-v for v in group.high)
            self.high = tuple(-v for v in group.low)
        self.base_num = sum(b*q*l for b, q, l in zip(anchor, self.q, self.low))
        self.rho = rho
        self.fixed = self.low == self.high

    def signs(self, i, j, t):
        mix_mass = t*self.q[i] + (1-t)*self.q[j]
        return (t*self.base_mass-self.b[i]*mix_mass,
                (1-t)*self.base_mass-self.b[j]*mix_mass)

    def changes(self, i, j):
        if self.fixed:
            return []
        roots = []
        for left, right in zip(self.signs(i, j, 0), self.signs(i, j, 1)):
            if left*right < 0:
                roots.append(left/(left-right))
        return roots

    def curve(self, i, j, midpoint):
        """Numerator and denominator at t=0 and t=1 of one fractional term."""
        si, sj = self.signs(i, j, midpoint)
        ri = self.high[i] if si > 0 else self.low[i]
        rj = self.high[j] if sj > 0 else self.low[j]
        ni = self.b[i]*self.q[i]*(ri-self.low[i])
        nj = self.b[j]*self.q[j]*(rj-self.low[j])
        baseline = (self.base_num+ni+nj)/self.base_mass
        anchor_mass = (1-self.rho)*self.base_mass
        mi, mj = self.rho*self.q[i], self.rho*self.q[j]
        return (mj*(rj-baseline), mi*(ri-baseline), anchor_mass+mj, anchor_mass+mi)


def lerp(a, b, t):
    return (1-t)*a+t*b


def stationary(curves):
    # Derivative = c0/D0^2 + c1/D1^2. Positive denominators
    # permit taking a square root without introducing a spurious root.
    c, d = [n1*d0-n0*d1 for n0, n1, d0, d1 in curves]
    if c*d >= 0:
        return None
    x, y = curves
    left = sqrt(abs(c))*y[2]-sqrt(abs(d))*x[2]
    right = sqrt(abs(c))*y[3]-sqrt(abs(d))*x[3]
    return left/(left-right) if left*right < 0 else None


def maximize(anchor, groups, rho, direction):
    if rho == 0 or len(anchor) == 1:
        return Witness(0, 0, 0, 1), 0
    terms = [OrientedGroup(anchor, g, sign, rho)
             for g, sign in zip(groups, (direction, -direction))]
    best = Witness(float("-inf"), 0, 0, 1)
    single = float("-inf")
    for i, j in combinations(range(len(anchor)), 2):
        cuts = sorted({0, 1, *(v for term in terms for v in term.changes(i, j))})
        for a, b in zip(cuts, cuts[1:]):
            curves = [term.curve(i, j, (a+b)/2) for term in terms]
            root = stationary(curves)
            candidates = [a, b]
            if root is not None and a < root < b:
                candidates.append(root)
            for t in candidates:
                val = sum(lerp(n0,n1,t)/lerp(d0,d1,t) for n0,n1,d0,d1 in curves)
                if t == 0 or t == 1:
                    single = max(single, val)
                if val > best.value:
                    best = Witness(val, i, j, t)
    return Witness(direction*best.value, best.first, best.second, best.fraction_first), direction*single


def envelope(anchor, groups, rho):
    check(anchor, groups, rho)
    lower, single_lower = maximize(anchor, groups, rho, -1)
    upper, single_upper = maximize(anchor, groups, rho, 1)
    return Envelope(lower, upper, single_lower, single_upper)


def independent(anchor, groups, rho):
    """Range the target and anchor gaps over separate latent outcome tables."""
    check(anchor, groups, rho)
    width = sum(sum(b*q*(h-l) for b,q,l,h in zip(anchor,g.share,g.low,g.high))
                / dot(anchor,g.share) for g in groups)
    bounds = []
    for direction in (-1, 1):
        fixed = []
        for g, sign in zip(groups, (direction, -direction)):
            r = g.high if sign == 1 else g.low
            fixed.append(Group(g.share, r, r))
        change, _ = maximize(anchor, fixed, rho, direction)
        bounds.append(change.value + direction*width)
    return tuple(bounds)


@dataclass(frozen=True)
class Bracket:
    inner_lower: float
    inner_upper: float
    allowance: float

    @property
    def lower(self):
        return self.inner_lower-self.allowance

    @property
    def upper(self):
        return self.inner_upper+self.allowance


def linear_certificate(anchor, groups, rho):
    """Linear-time outer endpoints from single-recipient corners plus curvature."""
    check(anchor, groups, rho)
    correction = 0
    for g in groups:
        base = dot(anchor,g.share)
        d_lo = (1-rho)*base + rho*min(g.share)
        d_hi = (1-rho)*base + rho*max(g.share)
        # Stable form of (sqrt(d_hi)-sqrt(d_lo))/(sqrt(d_hi)+sqrt(d_lo)).
        curvature = rho*(max(g.share)-min(g.share))/(sqrt(d_hi)+sqrt(d_lo))**2
        correction += curvature*(rho*max(g.share)/d_hi)*(max(g.high)-min(g.low))
    endpoints = []
    for direction in (-1,1):
        terms = [OrientedGroup(anchor,g,s,rho) for g,s in zip(groups,(direction,-direction))]
        best = 0.0  # the anchor is attainable at every radius
        for i in range(len(anchor)):
            value = 0
            for term in terms:
                # Under replacement by site i its own conditional weight increases.
                upper = term.high[i]
                base = (term.base_num+anchor[i]*term.q[i]*(upper-term.low[i]))/term.base_mass
                mass = rho*term.q[i]
                value += mass*(upper-base)/((1-rho)*term.base_mass+mass)
            best = max(best,value)
        endpoints.append(direction*best)
    return Bracket(*endpoints, correction)
