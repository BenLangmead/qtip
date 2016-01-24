import math
import numpy as np


def pcor_to_mapq_np(pcor):
    old = np.seterr(divide='ignore')
    ret = np.abs(-10.0 * np.log10(1.0 - pcor))
    np.seterr(**old)
    return ret


def mapq_to_pcor_np(mapq):
    return 1.0 - 10.0 ** (-0.1 * mapq)


def round_pcor_np(pcor):
    return mapq_to_pcor_np(np.round(pcor_to_mapq_np(pcor)))


def pcor_to_mapq(p):
    """ Convert probability correct (pcor) to mapping quality (MAPQ) """
    return abs(-10.0 * math.log10(1.0 - p)) if p < 1.0 else float('inf')


def mapq_to_pcor(p):
    """ Convert mapping quality (MAPQ) to probability correct (pcor) """
    return (1.0 - (10.0 ** (-0.1 * p))) if p < float('inf') else 1.0


def round_pcor(p):
    return mapq_to_pcor(round(pcor_to_mapq(p)))
