"""Shared-table outer bounds from single-ratio, TV and Bernstein inequalities."""
from math import sqrt
from audit import OrientedGroup, check

LINEAR = ("groupwise", "tv_transfer", "reciprocal")
METHODS = (*LINEAR, "linear_intersection", "bernstein")


def tv_factor(low, high):
    return (high-low)/(sqrt(high)+sqrt(low))**2


class VertexGroup:
    def __init__(self, anchor, group, sign, rho):
        term = OrientedGroup(anchor, group, sign, rho)
        self.low = term.low
        self.width = [h-l for l,h in zip(term.low, term.high)]
        self.base = term.base_num/term.base_mass
        self.share = [b*q/term.base_mass for b,q in zip(anchor,group.share)]
        self.mass = [rho*q for q in group.share]
        self.den = [(1-rho)*term.base_mass+m for m in self.mass]
        self.upper = [m*(l-self.base+(1-s)*width)/d for m,l,s,width,d in
                      zip(self.mass,self.low,self.share,self.width,self.den)]

    def cross_numerator(self, i, j, other_den):
        # Signed rate coefficients at i,j can switch sign. All others are
        # nonpositive after the common baseline mean has been subtracted.
        p = self.mass[i]*other_den[j]
        q = self.mass[j]*other_den[i]
        return (p*(self.low[i]-self.base) + q*(self.low[j]-self.base)
                + max(0.,p-(p+q)*self.share[i])*self.width[i]
                + max(0.,q-(p+q)*self.share[j])*self.width[j])


def reciprocal_bound(ref, other, center):
    ratios = [a/b for a,b in zip(ref.den,other.den)]
    low,high = min(ratios),max(ratios)
    # D_ref / D_other is a weighted average of the vertex ratios.
    return max(a+center+max(low*(b-center), high*(b-center))/r
               for a,b,r in zip(ref.upper,other.upper,ratios))


def oriented_bounds(terms, reverse, methods):
    a,b = terms
    inner = max(0.,max(x+y for x,y in zip(a.upper,b.upper)))
    bounds = {}
    if "groupwise" in methods:
        bounds["groupwise"] = max(a.upper)+max(b.upper)
    if "tv_transfer" in methods or "reciprocal" in methods:
        ranges = [(-max(neg.upper),max(pos.upper)) for pos,neg in zip(terms,reverse)]
        spans = [h-l for l,h in ranges]
        if "tv_transfer" in methods:
            ratios = [x/y for x,y in zip(a.den,b.den)]
            to_group = tv_factor(min(ratios),max(ratios))*min(spans)
            to_mix = sum(tv_factor(min(t.den),max(t.den))*span for t,span in zip(terms,spans))
            bounds["tv_transfer"] = inner+min(to_group,to_mix)
        if "reciprocal" in methods:
            bounds["reciprocal"] = min(
                reciprocal_bound(terms[1-j],terms[j],center)
                for j,(l,h) in enumerate(ranges) for center in (0.,(l+h)/2))
    if "bernstein" in methods:
        outer = inner  # diagonal coefficients are exactly the vertex changes
        for i in range(len(a.den)):
            for j in range(i):
                denom = a.den[i]*b.den[j]+a.den[j]*b.den[i]
                num = a.cross_numerator(i,j,b.den)+b.cross_numerator(i,j,a.den)
                outer = max(outer,num/denom)
        bounds["bernstein"] = outer
    return inner,bounds


def bounds(anchor, groups, rho, methods=METHODS):
    check(anchor,groups,rho)
    wanted = set(methods)
    if "linear_intersection" in wanted:
        wanted.update(LINEAR)
    both = [tuple(VertexGroup(anchor,g,sign,rho) for g,sign in zip(groups,(d,-d)))
            for d in (-1,1)]
    results = {}
    for d,terms,reverse in ((-1,both[0],both[1]),(1,both[1],both[0])):
        inner,outer = oriented_bounds(terms,reverse,wanted)
        if "linear_intersection" in wanted:
            outer["linear_intersection"] = min(outer[k] for k in LINEAR)
        for name,v in {"inner":inner,**outer}.items():
            results.setdefault(name,[]).append(d*v)
    return results
