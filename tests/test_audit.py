import itertools
import math
import random
import unittest

from audit import Group, envelope, independent
from study import panels, rounded


def mean(weights, group, rates):
    q = group.share
    return sum(b*p*r for b,p,r in zip(weights,q,rates))/sum(b*p for b,p in zip(weights,q))


def direct(anchor, u, groups, rho, direction):
    weights = [(1-rho)*b+rho*v for b,v in zip(anchor,u)]
    total = 0
    for g, sign in zip(groups,(direction,-direction)):
        f0 = sum(b*p for b,p in zip(anchor,g.share))
        f1 = sum(b*p for b,p in zip(weights,g.share))
        rates = [h if sign*(w/f1-b/f0)>0 else l
                 for b,w,l,h in zip(anchor,weights,g.low,g.high)]
        total += sign*(mean(weights,g,rates)-mean(anchor,g,rates))
    return direction*total


def max_quadratic_triangle(f):
    """Fit a quadratic in barycentric coordinates then inspect all KKT faces."""
    c = f(0,0)
    xx = 2*(f(1,0)+c-2*f(.5,0))
    yy = 2*(f(0,1)+c-2*f(0,.5))
    x, y = f(1,0)-xx-c, f(0,1)-yy-c
    xy = 4*(f(.5,.5)-c-.25*(xx+yy)-.5*(x+y))
    values = [c,f(1,0),f(0,1)]
    for at in (lambda t:f(t,0),lambda t:f(0,t),lambda t:f(t,1-t)):
        a = 2*(at(0)+at(1)-2*at(.5))
        b = at(1)-at(0)-a
        if a<0 and 0 < -b/(2*a) < 1:
            values.append(at(-b/(2*a)))
    det = 4*xx*yy-xy*xy
    if det:
        s, t = (xy*y-2*yy*x)/det, (xy*x-2*xx*y)/det
        if s>0 and t>0 and s+t<1:
            values.append(f(s,t))
    return max(values)


class TestAudit(unittest.TestCase):
    def test_global_oracle_and_witness(self):
        rng = random.Random(20261123)
        for rep in range(20):
            q = tuple(rng.uniform(.04,.96) for _ in range(3))
            groups = []
            for share in (q,tuple(1-p for p in q)):
                pairs = [sorted((rng.random(),rng.random())) for _ in q]
                groups.append(Group(share,*(tuple(v) for v in zip(*pairs))))
            anchor = (.25,.55,.2) if rep%2 else (0,.55,.45)
            rho = (1,.01,.3,.8)[rep%4]
            res = envelope(anchor,groups,rho)
            for orient,witness in ((1,res.upper),(-1,res.lower)):
                weights = [0]*3
                weights[witness.first] += witness.fraction_first
                weights[witness.second] += 1-witness.fraction_first
                self.assertAlmostEqual(direct(anchor,weights,groups,rho,orient),witness.value,places=12)
                # Fix every corner separately. After multiplication by both
                # positive group masses the excess over the candidate is quadratic.
                # This oracle includes interiors of the full three-site simplex.
                max_excess = -math.inf
                for bits in itertools.product((0,1),repeat=6):
                    a,b = groups
                    r,s = [tuple((l,h)[bit] for l,h,bit in zip(g.low,g.high,choice))
                           for g,choice in zip(groups,(bits[:3],bits[3:]))]
                    baseline = mean(anchor,a,r)-mean(anchor,b,s)
                    def excess(x,y):
                        w = [(1-rho)*bi+rho*ui for bi,ui in zip(anchor,(x,y,1-x-y))]
                        d0 = sum(wi*qi for wi,qi in zip(w,a.share))
                        d1 = sum(wi*qi for wi,qi in zip(w,b.share))
                        return orient*d0*d1*(mean(w,a,r)-mean(w,b,s)-baseline-witness.value)
                    max_excess = max(max_excess,max_quadratic_triangle(excess))
                self.assertAlmostEqual(max_excess,0,places=11)

    def test_interior_required(self):
        q = (.9,.1)
        g = [Group(p,(1,0),(1,0)) for p in (q,q[::-1])]
        a = (1,0)
        res = envelope(a,g,1)
        self.assertAlmostEqual(res.upper.value,.8)
        self.assertAlmostEqual(res.single_upper,0)

    def test_zero_and_one_site(self):
        g = [Group((q,),(.1,),(.9,)) for q in (.3,.7)]
        for rho in (0,.1,1):
            e = envelope((1,),g,rho)
            self.assertEqual((e.lower.value,e.upper.value),(0,0))
            self.assertAlmostEqual(independent((1,),g,rho)[1],1.6)

    def test_known_equal_composition_and_invalid(self):
        g = [Group((.5,.5),(0,.1),(.8,1))]*2
        e = envelope((.5,.5),g,.25)
        self.assertAlmostEqual(e.upper.value,e.single_upper)
        with self.assertRaises(ValueError): envelope((1,1),g,.1)
        with self.assertRaises(ValueError): envelope((.5,.5),g,1.5)

    def test_rounding_preimage(self):
        for step,n in itertools.product((1,5,10,100),(7,19,209)):
            buckets = {}
            for y in range(n+1):
                buckets.setdefault((2*y+step)//(2*step),[]).append(y)
            for values in buckets.values():
                for y in values:
                    self.assertEqual(rounded(y,n,step),(min(values),max(values)))

    def test_source_selection(self):
        for year in (2022,2023):
            selected,flow = panels(year)
            self.assertEqual(len(selected),30)
            ids = [r[0] for _,_,rows in selected for r in rows]
            self.assertEqual(len(set(ids)),240)
            self.assertEqual(sum(v for k,v in flow.items() if k!="raw"),flow["raw"])

class TestCurvature(unittest.TestCase):
    def test_second_order_is_necessary(self):
        from audit import linear_certificate
        a = (0,0,1)
        rates = (1,0,.5)
        g = [Group(q,rates,rates) for q in ((.9,.1,.5),(.1,.9,.5))]
        for rho in (1,.25,.1,.01,.001):
            res = envelope(a,g,rho)
            self.assertAlmostEqual(res.upper.value,.8*rho)
            self.assertAlmostEqual(res.single_upper,.8*rho*(1-rho)/(1-.64*rho*rho))
            bracket = linear_certificate(a,g,rho)
            lo,hi = bracket.lower,bracket.upper
            self.assertAlmostEqual(bracket.inner_upper,res.single_upper)
            self.assertGreaterEqual(hi+1e-12,res.upper.value)
            self.assertLessEqual(lo-1e-12,res.lower.value)
