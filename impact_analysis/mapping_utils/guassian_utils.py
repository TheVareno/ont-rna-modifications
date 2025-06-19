
"""
author: Hadi Vareno
e-mail: mohammad.noori.vareno@uni-jena.de
github: https://github.com/TheVareno
"""

import numpy as np
        
def gaussian(x, mu, sigma, amplitude):
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def cal_kurtosis(data, mean, std):
    return (np.sum(((data - mean) / std) ** 4) / len(data)) - 3
      
def cal_skewness(data, mean, std): 
    return np.sum(((data - mean) / std) ** 3) / len(data)   
