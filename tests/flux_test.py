import unittest

import sys
import os

#Add src directory to the path so python can find the files
sys.path.append(os.path.abspath("../src"))
import flux

class FluxTest(unittest.TestCase):
    '''Unit test for the Flux object'''    
    def test(self):

        self.flux = flux.Flux("test_flux.txt")

if __name__ == '__main__':

    unittest.main()
