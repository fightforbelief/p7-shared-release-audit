import unittest
from random import Random
from audit import Group, envelope
from shared_baselines import METHODS, bounds


class SharedBaselines(unittest.TestCase):
    def test_signed_enclosures(self):
        rng = Random(1203)
        for k in (1,2,4,7):
            for same_share in (True,False):
                q = tuple(.22 if same_share else rng.uniform(.01,.99) for _ in range(k))
                groups = []
                for shares in (q,tuple(1-x for x in q)):
                    rates = [sorted((rng.random(),rng.random())) for _ in q]
                    groups.append(Group(shares,*map(tuple,zip(*rates))))
                for anchor in ([1/k]*k,[1]+[0]*(k-1)):
                    for rho in (0,.01,.5,1):
                        exact = envelope(anchor,groups,rho)
                        result = bounds(anchor,groups,rho)
                        for v,w in zip(result["inner"],(min(0,exact.single_lower),max(0,exact.single_upper))):
                            self.assertAlmostEqual(v,w)
                        for name in METHODS:
                            lo,hi = result[name]
                            self.assertLessEqual(lo,exact.lower.value+1e-11)
                            self.assertGreaterEqual(hi,exact.upper.value-1e-11)
                            if same_share and name in ("tv_transfer","reciprocal","bernstein"):
                                self.assertAlmostEqual(lo,exact.lower.value)
                                self.assertAlmostEqual(hi,exact.upper.value)
