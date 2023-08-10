# General
import os
import sys
import argparse
import itertools
import math
from array import *
import numpy


################################################################################
def eachJetPtBinDef(valKind, centBin):
    # mcGen: determined to become lager range than detector level range
    # mcDet: Low edge is determined from 5 sigma of delta pT distribution
    #        Upper edge is determined cause the number of entry for unfold bins
    # reported: determined by kinematic efficiency (over 80% efficiency)
    # 0: ptMin, 1:ptMax 
    ptRangeDict = None
    binArrayGen = None
    binArrayDet = None
    binArrayDetL = None
    binArrayDetU = None
    binArrayFullDet = None
    
    if valKind == 0:
        if centBin== 0:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,120]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 100, 120])
            # binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120])
            # binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([60, 70, 80, 100, 120])
            binArrayGen = ([20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200])
            binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250])
            binArrayReport = ([30, 40, 50, 60, 70, 80, 100, 120])
        elif centBin == 1:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,120]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200])
            # binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 140])
            # binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140])
            # binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120, 140])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([50, 60, 70, 80, 100, 120])
            binArrayGen = ([20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200])
            binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250])
            binArrayReport = ([30, 40, 50, 60, 70, 80, 100, 120])
        elif centBin == 2:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,120]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            # binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250])
            # binArrayReport = ([35, 40, 50, 60, 70, 80, 100, 120, 140])
            binArrayGen = ([20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200])
            binArrayDet = ([30, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetL = ([25, 30, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250])
            binArrayReport = ([30, 40, 50, 60, 70, 80, 100, 120])
        elif centBin == 3:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,120]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 140])

            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])

            # binArrayGen = ([15, 20, 30, 40, 60, 80, 100, 140, 200])
            # binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            # binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            # binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            # binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([30, 35, 40, 50, 60, 70, 80, 100, 140])

            binArrayGen = ([15, 20, 30, 40, 50, 60, 70, 80, 95, 120, 140])
            binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 95, 120])
            binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120])
            binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 95, 120, 140, 200, 250])
            # binArrayReport = ([40, 50, 60,  80, 100, 120])
            binArrayReport = ([30, 40, 50, 60, 70, 80, 95, 120])
        elif centBin == 4:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,140]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 140])
            binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            binArrayReport = ([30, 35, 40, 50, 60, 70, 80, 100, 120])
    elif valKind == 1:
        if centBin== 0:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,120]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 100, 120])
            # binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120])
            # binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([60, 70, 80, 100, 120])
            binArrayGen = ([20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200])
            binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250])
            binArrayReport = ([30, 40, 50, 60, 70, 80, 100, 120])
        elif centBin == 1:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,120]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200])
            # binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 140])
            # binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140])
            # binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120, 140])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([50, 60, 70, 80, 100, 120])
            binArrayGen = ([20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200])
            binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250])
            binArrayReport = ([30, 40, 50, 60, 70, 80, 100, 120])
        elif centBin == 2:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,120]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            # binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250])
            # binArrayReport = ([35, 40, 50, 60, 70, 80, 100, 120, 140])
            binArrayGen = ([20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200])
            binArrayDet = ([30, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetL = ([25, 30, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120, 150])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250])
            binArrayReport = ([30, 40, 50, 60, 70, 80, 100, 120])
        elif centBin == 3:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,120]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 140])

            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])

            # binArrayGen = ([15, 20, 30, 40, 60, 80, 100, 140, 200])
            # binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            # binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            # binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            # binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([30, 35, 40, 50, 60, 70, 80, 100, 140])

            binArrayGen = ([15, 20, 30, 40, 50, 60, 70, 80, 95, 120, 140])
            binArrayDet = ([30, 35, 40, 50, 60, 70, 80, 95, 120])
            binArrayDetL = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120])
            binArrayDetU = ([35, 40, 50, 60, 70, 80, 100, 120])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 95, 120, 140, 200, 250])
            # binArrayReport = ([40, 50, 60,  80, 100, 120])
            binArrayReport = ([30, 40, 50, 60, 70, 80, 95, 120])
        elif centBin == 4:
            ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[30,140]}
            # binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200])
            # binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            # binArrayReport = ([30, 35, 40, 50, 60, 70, 80, 100, 120, 140])
            binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 140])
            binArrayFullDet = ([0, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 200, 250])
            binArrayReport = ([30, 35, 40, 50, 60, 70, 80, 100, 120])
    else:
        ptRangeDict = {'mcGen':[15,200], 'mcDet':[20,150], 'meas':[40,120], 'reported':[40,120]}
        binArrayGen = ([10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 190, 250])
        binArrayDet = ([20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
        binArrayDetL = ([15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
        binArrayDetU = ([25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150])
        binArrayFullDet = ([0, 20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 250])
        binArrayReport = ([40, 50, 60, 70, 80, 100, 120])

    binArrayGen = array('d',binArrayGen)
    binArrayDet = array('d',binArrayDet)
    binArrayDetU = array('d',binArrayDetU)
    binArrayDetL = array('d',binArrayDetL)
    binArrayFullDet = array('d',binArrayFullDet)
    binArrayReport = array('d',binArrayReport)

    ptBinArrayDict = {'mcGen':binArrayGen, 'mcDet':binArrayDet, \
        'mcDetL':binArrayDetL, 'mcDetU':binArrayDetU, \
            'mcFullDet':binArrayFullDet, 'reported':binArrayReport} 

    return ptRangeDict, ptBinArrayDict
################################################################################