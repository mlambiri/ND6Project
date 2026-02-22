#!/usr/bin/env python3
# written 20210527 by Jack
# edited 20211211

import argparse
import numpy as np
# from math import science as chem

gas_constant = 8.3144626 / 4184  # kcal/K/mol


def bar_analysis(forwardfilepath, reversefilepath, maxiter, convergencecriteria, temp):
    denergylists = []
    for filepath in forwardfilepath, reversefilepath:
        with open(filepath) as file:
            startdata = False
            denergylist = []
            for line in file:
                if line.startswith('#STARTING COLLECTION OF ENSEMBLE AVERAGE'):
                    startdata = True
                elif line.startswith('#Free energy change for lambda window'):
                    break
                elif startdata:
                    denergylist.append(float(line.split()[6]))
        denergylists.append(denergylist)
    cold = 0
    for i in range(maxiter):
        ave = [0, 0]
        for j in range(2):
            for denergy in denergylists[j]:
                ave[j] += 1 / (1 + np.exp((denergy - (-1) ** j * cold) / gas_constant / temp))  # maybe
            ave[j] /= len(denergylists[j])
        cnew = cold - gas_constant * temp * np.log(ave[0]/ave[1])
        if abs(cnew - cold) < convergencecriteria:
            break
        else:
            cold = cnew
    else:
        raise Exception('oof')

    return cnew


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='do a Bennett Acceptance Ratio analysis to get delta G')

    parser.add_argument('forwardfilepath', type=str, help='path to forward fep file')
    parser.add_argument('reversefilepath', type=str, help='path to reverse fep file')
    parser.add_argument('-m', '--maxiter', type=int, default=100, help='maximum iterations in BAR analysis')
    parser.add_argument('-c', '--criteria', type=float, default=.01, help='convergence criteria in kcal/mol')
    parser.add_argument('-t', '--temperature', type=int, default=298, help='temperature at which simulations ran')

    args = parser.parse_args()

    dg = bar_analysis(args.forwardfilepath, args.reversefilepath, args.maxiter, args.criteria, args.temperature)

    print(dg)
