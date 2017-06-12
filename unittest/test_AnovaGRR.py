""" UNIT TEST ON ANOVA GAGE R&R PYTHON MODULE
# Description:
    This is the unit test for ANOVA Gage R&R module.
# Dependencies: Numpy
# Author: Shin-Fu (Kelvin) Wu
# Date: 2017/06/08
# Reference:
    * https://www.rdocumentation.org/packages/qualityTools/versions/1.31.1/topics/gageRRDesign
"""

import os
import sys
import unittest
import numpy as np

root = os.path.join(os.path.dirname(__file__), '..')
sys.path.append(root)
from algorithms.GageRR.AnovaGRR import AnovaGRR

class Test(unittest.TestCase):
    
    def __init__(self, methodName='runTest'):
        super().__init__(methodName)
        M = [23,22,22,22,22,25,23,22,23,22,20,22,22,22,24,25,27,28,23,24,23,24,24,\
             22,22,22,24,23,22,24,20,20,25,24,22,24,21,20,21,22,21,22,21,21,24,27,\
             25,27,23,22,25,23,23,22,22,23,25,21,24,23]
        measurement = {}
        for app in ['A', 'B', 'C']:
            measurement[app] = {}
        
        for app in ['A', 'B', 'C']:
            for i in range(1,11):
                measurement[app][i]=[]
        
        n=1
        for i in range(0, len(M), 3):
            measurement['A'][n].append(M[i])
            measurement['B'][n].append(M[i+1])
            measurement['C'][n].append(M[i+2])
            n+=1
            if n > 10:
                n -= 10
        self.grr = AnovaGRR(measurement, sig=6, tolerance=10)
    def testANOVA(self):
        np.testing.assert_equal(self.grr.ANOVAtable.loc['Operator', 'DF'], 2)
        np.testing.assert_almost_equal(self.grr.ANOVAtable.loc['Operator', 'SS'], 20.63, 2)
        np.testing.assert_almost_equal(self.grr.ANOVAtable.loc['Operator', 'MS'], 10.317, 3)        
        np.testing.assert_equal(self.grr.ANOVAtable.loc['Part', 'DF'], 9)
        np.testing.assert_almost_equal(self.grr.ANOVAtable.loc['Part', 'SS'], 107.07, 2)
        np.testing.assert_almost_equal(self.grr.ANOVAtable.loc['Part', 'MS'], 11.896, 3)
        np.testing.assert_equal(self.grr.ANOVAtable.loc['Operator by Part', 'DF'], 18)
        np.testing.assert_almost_equal(self.grr.ANOVAtable.loc['Operator by Part', 'SS'], 22.03, 2)
        np.testing.assert_almost_equal(self.grr.ANOVAtable.loc['Operator by Part', 'MS'], 1.224, 3)
        np.testing.assert_equal(self.grr.ANOVAtable.loc['Equipment', 'DF'], 30)
        np.testing.assert_almost_equal(self.grr.ANOVAtable.loc['Equipment', 'SS'], 36, 2)
        np.testing.assert_almost_equal(self.grr.ANOVAtable.loc['Equipment', 'MS'], 1.2, 3)
    def testGRR(self):
        desired = np.array([[1.664, 0.483, 1.290, 7.74, 0.695, 0.774],
                            [1.209, 0.351, 1.100, 6.60, 0.592, 0.660],
                            [0.455, 0.132, 0.675, 4.05, 0.364, 0.405],
                            [0.455, 0.132, 0.675, 4.05, 0.364, 0.405],
                            [0.000, 0.000, 0.000, 0.00, 0.000, 0.000],
                            [1.781, 0.517, 1.335, 8.01, 0.719, 0.801],
                            [3.446, 1.000, 1.856, 11.14, 1.000, 1.114]])
        np.testing.assert_array_almost_equal( np.round(self.grr.GRRtable.values[:, -1], 2), desired[:, -1], 1)
    
if __name__ == '__main__':
    unittest.main(verbosity=1)    