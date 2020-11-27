# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 09:00:53 2020

@author: jihon
"""

import numpy as np
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import numpy2ri

numpy2ri.activate()
robjects.r('''source('R/calculateApexMaxBoundary.R')''')
robjects.r('''source('R/calculateFWHM.R')''')
robjects.r('''source('R/calculateGaussianSimilarity.R')''')
robjects.r('''source('R/calculateJaggedness.R')''')
robjects.r('''source('R/calculateModality.R')''')
robjects.r('''source('R/calculateSharpness.R')''')
robjects.r('''source('R/calculateSymmetry.R')''')
robjects.r('''source('R/calculateTPASR.R')''')
robjects.r('''source('R/calculateZigZagIndex.R')''')

calculateApexMaxBoundaryRatio = robjects.globalenv['calculateApexMaxBoundaryRatio']
calculateFWHM = robjects.globalenv['calculateFWHM']
calculateGaussianSimilarity = robjects.globalenv['calculateGaussianSimilarity']
calculateJaggedness = robjects.globalenv['calculateJaggedness']
calculateModality = robjects.globalenv['calculateModality']
calculateSharpness = robjects.globalenv['calculateSharpness']
calculateSymmetry = robjects.globalenv['calculateSymmetry']
calculateTPASR = robjects.globalenv['calculateTPASR']
calculateZigZagIndex = robjects.globalenv['calculateZigZagIndex']

class PeakEval:
    
    def __init__(self, x, rtmin, rtmax):
        self.x = x
        self.rtmin = rtmin
        self.rtmax = rtmax
        
    def ApexMaxBoundaryRatio(self):
        return calculateApexMaxBoundaryRatio(self.x, self.rtmin, self.rtmax)[0]
    
    def FWHM(self):
        return calculateFWHM(self.x, self.rtmin, self.rtmax)[0]
    
    def GaussianSimilarity(self):
        return calculateGaussianSimilarity(self.x, self.rtmin, self.rtmax)[0]
    
    def Jaggedness(self):
        return calculateJaggedness(self.x, self.rtmin, self.rtmax)[0]
    
    def Modality(self):
        return calculateModality(self.x, self.rtmin, self.rtmax)[0]
    
    def Sharpness(self):
        return calculateSharpness(self.x, self.rtmin, self.rtmax)[0]
    
    def Symmetry(self):
        return calculateSymmetry(self.x, self.rtmin, self.rtmax)[0]
    
    def TPASR(self):
        return calculateTPASR(self.x, self.rtmin, self.rtmax)[0]    

    def ZigZagIndex(self):
        return calculateZigZagIndex(self.x, self.rtmin, self.rtmax)[0]   
    
    
if __name__ == '__main__':
    
    x = np.array(pd.read_csv('Example/XIC.csv'))
    # x[:,1] = x[:,1] / max(x[:,1])
    rtmin, rtmax = 22.41, 36.42
    
    peak_eval = PeakEval(x, rtmin, rtmax)
    print(peak_eval.ApexMaxBoundaryRatio())
    print(peak_eval.FWHM())
    print(peak_eval.GaussianSimilarity())
    print(peak_eval.Jaggedness())
    print(peak_eval.Modality())
    print(peak_eval.Sharpness())
    print(peak_eval.Symmetry())
    print(peak_eval.TPASR())
    print(peak_eval.ZigZagIndex())
    