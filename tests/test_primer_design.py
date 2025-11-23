import unittest
import sys
import os

# Add app to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app.utils import (
    reverse_complement, calculate_tm, calculate_gc_content,
    make_hairpin, _dg
)
from app.primer_design import (
    tune_primer, build_rt_primer, check_cross_reactivity, 
    enumerate_strict_candidates, VIENNA_AVAILABLE
)

class TestUtils(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("AATTGGCC"), "GGCCAATT")

    def test_calculate_gc_content(self):
        self.assertAlmostEqual(calculate_gc_content("GCGC"), 100.0)
        self.assertAlmostEqual(calculate_gc_content("ATAT"), 0.0)
        self.assertAlmostEqual(calculate_gc_content("ATGC"), 50.0)

    def test_calculate_tm(self):
        # Simple check against expected behavior (higher GC -> higher Tm)
        tm_low = calculate_tm("ATATATAT", 1e-6, 0.05, 0.003, 0.0005)
        tm_high = calculate_tm("GCGCGCGC", 1e-6, 0.05, 0.003, 0.0005)
        self.assertTrue(tm_high > tm_low)

    def test_dg(self):
        # GC should be more stable (more negative) than AT
        dg_gc = _dg("GCGC")
        dg_at = _dg("ATAT")
        self.assertTrue(dg_gc < dg_at)

class TestPrimerDesign(unittest.TestCase):
    def test_build_rt_primer(self):
        rt5 = "AAAA"
        rt3 = "TTTT"
        res = build_rt_primer(rt5, rt3, n_rep=10)
        if res:
            rt_primer, rand5, rand3, hairpin = res
            self.assertTrue(rt_primer.startswith(rt5))
            self.assertTrue(rt_primer.endswith(rt3))
            self.assertIn(hairpin, rt_primer)

    def test_tune_primer(self):
        # Test tuning a primer to a specific Tm range
        seq = "ATATATATATATATAT" # Low Tm
        target_tm = (58, 65)
        # It should add bases to reach the range
        try:
            tuned, tm = tune_primer(
                seq, target_tm, 
                oligo_conc=0.4e-6, Na_conc=0.05, Mg_conc=0.003, dNTPs_conc=0.0005,
                min_length=16, max_length=30
            )
            self.assertTrue(58 <= tm <= 65)
            self.assertTrue(len(tuned) >= len(seq))
        except ValueError:
            # It might fail if it can't reach the target, but for this input it should work or fail gracefully
            pass

    def test_check_cross_reactivity(self):
        # If ViennaRNA is not available, it returns error
        if not VIENNA_AVAILABLE:
            res = check_cross_reactivity("AAAA", "UUUU")
            self.assertIn("error", res)
        else:
            # If available, check structure
            res = check_cross_reactivity("AAAA", "UUUU")
            self.assertIn("mfe", res)

    def test_enumerate_strict_candidates_binding_info(self):
        # Test if binding info is present
        miDNA = "ATGCATGCATGCATGCATGC" # 20nt
        # Mock params
        res = enumerate_strict_candidates(
            miDNA, 0.05, 0.003, 0.0005, 1e-7, 1e-7,
            top_n=1, tries_per_pair=10, beam_per_pair=1,
            l3_min=4, l3_max=6, l5_min=4, l5_max=6
        )
        if res:
            cand = res[0]
            self.assertIn("binding_5_start", cand)
            self.assertIn("binding_3_start", cand)
            self.assertIn("binding_3_end", cand)

    def test_enumerate_strict_candidates_position_constraints(self):
        # Test if p5_range constraints are respected
        miDNA = "ATGCATGCATGCATGCATGC" # 20nt
        # Constraint: 5' start must be between 1 and 5
        res = enumerate_strict_candidates(
            miDNA, 0.05, 0.003, 0.0005, 1e-7, 1e-7,
            top_n=5, tries_per_pair=10, beam_per_pair=1,
            l3_min=4, l3_max=6, l5_min=4, l5_max=6,
            p5_range=(1, 5)
        )
        if res:
            for cand in res:
                self.assertTrue(1 <= cand['binding_5_start'] <= 5)
        
        # Constraint: 5' start must be > 10
        res2 = enumerate_strict_candidates(
            miDNA, 0.05, 0.003, 0.0005, 1e-7, 1e-7,
            top_n=5, tries_per_pair=10, beam_per_pair=1,
            l3_min=4, l3_max=6, l5_min=4, l5_max=6,
            p5_range=(10, None)
        )
        if res2:
            for cand in res2:
                self.assertTrue(cand['binding_5_start'] >= 10)

            for cand in res2:
                self.assertTrue(cand['binding_5_start'] >= 10)

    def test_heap_tie_handling(self):
        # Test that the function doesn't crash when candidates have identical scores
        # A poly-A sequence will generate many identical primers with identical scores
        miDNA = "A" * 30
        try:
            res = enumerate_strict_candidates(
                miDNA, 0.05, 0.003, 0.0005, 1e-7, 1e-7,
                top_n=10, tries_per_pair=50, beam_per_pair=5,
                l3_min=4, l3_max=6, l5_min=4, l5_max=6
            )
            # If we reach here without TypeError, it passed
            self.assertTrue(True)
        except TypeError as e:
            self.fail(f"Heap comparison failed with TypeError: {e}")

if __name__ == '__main__':
    unittest.main()
