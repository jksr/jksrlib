import numpy as np
#import pymannkendall as mk

def is_binarizable(array, diff_thresh=None, bin_thresh=None):
    _array = np.sort(array)
    between = np.diff(_array)
    within = np.empty(between.shape)
    
    for i in range(len(_array)-1):
        #within[i] = max( np.ptp(_array[:i+1]), np.ptp(_array[i+1:]) )
        within[i] = np.ptp(_array[:i+1]) + np.ptp(_array[i+1:])
    
    _binarizability = between / within
    binary_diff = between[ np.argmax(_binarizability) ]
    binarizability = max(_binarizability)
    
    if diff_thresh is not None and bin_thresh is not None:
        return binary_diff > diff_thresh and binarizability > bin_thresh
    else:
        return binary_diff, binarizability
    
# def has_mk_trend(array, p=None):
#     trend, *_ = mk.original_test(array, p)
#     return trend

def has_pearson_trend(array, rho=None):
    _rho = np.corrcoef( np.arange(len(array)), array )[0,1]
    if rho is None:
        return _rho
    
    rho = abs(rho)
    if _rho > rho:
        return 'increasing'
    elif _rho < -rho:
        return 'decreasing'
    else:
        return 'no trend'
