import numpy as np

def mols():
    return {"CuC": {"atoms": ["Cu", "C"], "M": [2, 4]}}

def rvals():
    return ["{:.2f}".format(r) for r in np.arange(1.8, 5.1, 0.2)]
