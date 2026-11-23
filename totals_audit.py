"""Sharp group-gap changes when each site's two admitted counts share a sum band."""
from dataclasses import dataclass
from itertools import combinations, product
from audit import Group, Witness, check, lerp, linear_certificate, stationary


@dataclass(frozen=True)
class Site:
    women: int
    men: int
    women_bounds: tuple
    men_bounds: tuple
    sum_bounds: tuple

    def __post_init__(self):
        if min(self.women, self.men) <= 0:
            raise ValueError("positive group denominators required")
        for (low, high), cap in zip((self.women_bounds, self.men_bounds, self.sum_bounds),
                                    (self.women, self.men, self.women+self.men)):
            if not (isinstance(low,int) and isinstance(high,int) and 0 <= low <= high <= cap):
                raise ValueError("invalid integer count interval")
        if not self.vertices():
            raise ValueError("empty site polygon")

    @property
    def total(self):
        return self.women+self.men

    def contains(self, x, y):
        a,b = self.women_bounds
        c,d = self.men_bounds
        l,h = self.sum_bounds
        return a <= x <= b and c <= y <= d and l <= x+y <= h

    def vertices(self):
        candidates = set(product(self.women_bounds, self.men_bounds))
        for total in self.sum_bounds:
            candidates.update((x,total-x) for x in self.women_bounds)
            candidates.update((total-y,y) for y in self.men_bounds)
        return sorted((x,y) for x,y in candidates if self.contains(x,y))


def projections(sites):
    groups = []
    for g in (0,1):
        q, low, high = [], [], []
        for s in sites:
            counts = [v[g] for v in s.vertices()]
            n = (s.women,s.men)[g]
            q.append(n/s.total)
            low.append(min(counts)/n)
            high.append(max(counts)/n)
        groups.append(Group(tuple(q),tuple(low),tuple(high)))
    return groups


def linear(anchor, sites, rho):
    # The two opposite-coordinate corners are jointly feasible at every site.
    return linear_certificate(anchor, projections(sites), rho)


def endpoints(anchor, sites, rho):
    groups = projections(sites)
    check(anchor,groups,rho)
    if rho == 0 or len(sites)==1:
        zero = Witness(0,0,0,1)
        return zero,zero
    vertices = [s.vertices() for s in sites]
    masses = [sum(b*q for b,q in zip(anchor,g.share)) for g in groups]
    anchor_scale = [b/s.total for b,s in zip(anchor,sites)]
    bounds = []
    for sign in (-1,1):
        # At non-recipients the two coefficients have opposite signs.
        default = [min(v,key=lambda point:(sign*point[0],-sign*point[1])) for v in vertices]
        baseline = [sum(f*v[g] for f,v in zip(anchor_scale,default)) for g in (0,1)]
        best = Witness(float("-inf"),0,0,1)
        for i,j in combinations(range(len(sites)),2):
            for vi,vj in product(vertices[i],vertices[j]):
                curves = []
                for g,orientation in ((0,sign),(1,-sign)):
                    updated = (baseline[g] + anchor_scale[i]*(vi[g]-default[i][g])
                               + anchor_scale[j]*(vj[g]-default[j][g]))
                    mean = updated/masses[g]
                    qi,qj = groups[g].share[i],groups[g].share[j]
                    denominator = (1-rho)*masses[g]
                    curves.append((orientation*rho*(vj[g]/sites[j].total-qj*mean),
                                   orientation*rho*(vi[g]/sites[i].total-qi*mean),
                                   denominator+rho*qj,denominator+rho*qi))
                root = stationary(curves)
                for t in ((0,1) if root is None else (0,root,1)):
                    val = sum(lerp(a,b,t)/lerp(c,d,t) for a,b,c,d in curves)
                    if val>best.value:
                        best = Witness(val,i,j,t)
        bounds.append(Witness(sign*best.value,best.first,best.second,best.fraction_first))
    return tuple(bounds)
