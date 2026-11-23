import itertools
import random
import unittest
from public_totals import direct, site_from_release
from totals_audit import Site, endpoints, linear
from test_audit import max_quadratic_triangle


class TotalChecks(unittest.TestCase):
    def test_polygon_linear_support(self):
        rng = random.Random(20261123)
        for _ in range(60):
            a,b = sorted(rng.sample(range(9),2))
            c,d = sorted(rng.sample(range(8),2))
            l,h = sorted(rng.choices(range(a+c,b+d+1),k=2))
            s = Site(8,7,(a,b),(c,d),(l,h))
            feasible = [(x,y) for x in range(9) for y in range(8) if s.contains(x,y)]
            for px,py in ((1,1),(1,-1),(-1,1),(2,3),(-3,-1),(0,1)):
                self.assertEqual(max(px*x+py*y for x,y in feasible),
                                 max(px*x+py*y for x,y in s.vertices()))

    def test_total_and_missing_groups(self):
        row=dict(APPLCN="42",APPLCNW="20",APPLCNM="15",APPLCNAN="",APPLCNUN="",
                 ADMSSN="27",ADMSSNW="14",ADMSSNM="10",ADMSSNAN="",ADMSSNUN="")
        self.assertEqual(site_from_release(row,True).sum_bounds,(20,27))
        row.update(APPLCNAN="3",APPLCNUN="4",ADMSSNAN="2",ADMSSNUN="1")
        self.assertEqual(site_from_release(row,True).sum_bounds,(20,27))
        row.update(APPLCN="60",APPLCNUN="22",ADMSSNUN="6")
        self.assertEqual(site_from_release(row,True).sum_bounds,(10,22))
        self.assertEqual(site_from_release(row,False).sum_bounds,(2,27))
        with self.assertRaises(ValueError):
            Site(2,2,(0,1),(0,1),(3,4))

    def test_full_simplex_integer_oracle(self):
        # Independent oracle checks the whole simplex and every feasible
        # integer table, including nonvertices of degenerate polygons.
        sites = [Site(4,2,(1,3),(0,2),(2,3)),Site(2,5,(0,2),(2,4),(4,4)),
                 Site(3,3,(1,3),(0,1),(2,3))]
        options = [[(x,y) for x in range(s.women+1) for y in range(s.men+1)
                    if s.contains(x,y)] for s in sites]
        for rho,anchor in ((0,(.2,.3,.5)),(.1,(0,.5,.5)),(.7,(.2,.3,.5)),(1,(.2,.3,.5))):
            low,high = endpoints(anchor,sites,rho)
            scan = linear(anchor,sites,rho)
            self.assertLessEqual(scan.lower,low.value+1e-12)
            self.assertGreaterEqual(scan.upper,high.value-1e-12)
            for sign,w in ((-1,low),(1,high)):
                self.assertAlmostEqual(direct(anchor,sites,rho,w,sign),w.value,places=12)
                best_excess = -float("inf")
                for table in itertools.product(*options):
                    def components(weights):
                        g = sum(v*s.women/s.total for v,s in zip(weights,sites))
                        x = sum(v*p[0]/s.total for v,p,s in zip(weights,table,sites))
                        y = sum(v*p[1]/s.total for v,p,s in zip(weights,table,sites))
                        return g,x,y
                    g,x,y = components(anchor)
                    base = x/g-y/(1-g)
                    def excess(a,b):
                        weights = [(1-rho)*v+rho*u for v,u in zip(anchor,(a,b,1-a-b))]
                        m,x,y = components(weights)
                        return sign*(x*(1-m)-y*m-(base+w.value)*m*(1-m))
                    best_excess=max(best_excess,max_quadratic_triangle(excess))
                self.assertAlmostEqual(best_excess,0,places=12)
