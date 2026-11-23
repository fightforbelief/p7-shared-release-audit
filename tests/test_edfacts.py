import unittest
from edfacts_audit import count_band, fixed_bounds
from totals_audit import Site


class EDFactsChecks(unittest.TestCase):
    def test_bands_and_ties(self):
        self.assertEqual(count_band(200,"42"),(83,85))
        self.assertEqual(count_band(200,"42",True),(83,85))
        self.assertEqual(count_band(1000,"42",True),(415,429))
        self.assertEqual(count_band(20,"LT50"),(0,9))
        self.assertEqual(count_band(20,"GE50"),(10,20))
        self.assertEqual(count_band(380,"PS"),(0,380))
        self.assertEqual(count_band(100,"25-29"),(25,29))
        with self.assertRaises(ValueError):
            count_band(100,".")

    def test_fixed_policy_at_anchor(self):
        sites=[Site(10,40,(0,8),(0,40),(10,12)),Site(25,5,(0,25),(0,5),(2,9))]
        self.assertEqual(fixed_bounds(sites,[.3,.7],[.3,.7]),(0,0,0,0))
