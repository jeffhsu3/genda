
"""Testing for overlap intervals
"""
import unittest
from genda.transcripts.exon_utils import calcOverlap, collideIntervals, \
        collapseIntervals

class TestOverlapFunctions(unittest.TestCase):

    def setUp(self):
        # Simple Overlap
        self.simple = [(1,10), (6,15)]
        # One interval enclosed in another
        self.enclosed = [(100,200), (110,150)]
        # Partial overlap
        self.partial = [(150,300), (160,300), (170,330)]
        # No overlap
        self.no = [(150,300), (10,30)]
        # Equal
        self.equal = [(1,15), (1,5)]
        #Complex interval list 
        self.full = [(7,20), (1,5), (8,11), (18,50), (100,150)]

    def test_bpOverlap(self):
        # Make sure overlaps are calculated correctly
        self.assertEqual(calcOverlap(self.simple), 4)
        self.assertEqual(calcOverlap(self.enclosed), 40)
        self.assertEqual(calcOverlap(self.partial),400)

    def test_collideIntervals(self):
        self.assertEqual(collideIntervals(self.simple[0], self.simple[1]),
                [(1,15)])
        self.assertEqual(collideIntervals(self.enclosed[0], self.enclosed[1]),
                [(100,200)])
        self.assertEqual(collideIntervals(self.no[0], self.no[1]),self.no)

    def test_collapseIntervals(self):
        self.assertEqual(collapseIntervals(self.simple), [(1,15)])
        print(self.partial)
        self.assertEqual(collapseIntervals(self.partial), [(150,330)])
        print(self.full)
        self.assertEqual(collapseIntervals(self.full), [(1,5),(7,50),(100,150)])

    def test_unique_bp(self):
        self.assertEqual(sum(map(lambda x \
                :x[1]-x[0],collapseIntervals(self.partial))) -
                         calcOverlap(self.partial),330-150)


if __name__ == '__main__':
    unittest.main()
