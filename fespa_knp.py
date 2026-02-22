#!/usr/bin/env python3
# use kth nearest point to calculate kT*ln(rho)
# copied from fespameta.py 20211209
# edited 20211215

import argparse
from math import lgamma
import numpy as np
import os
import os.path

gas_constant = 8.3144626 / 4184


def digamma(x):
    # digamma(x) = L_x-1 - gamma
    # L_x-1 = sum(i = 1 to x - 1)[1 / i] if x > 1, or 0 if x = 1
    result = -np.euler_gamma
    for i in range(1, x):
        result += 1/i
    return result


def main(jobname, temp, kmax, summary, numreps, convergence):

    jobname = jobname.split('.')[0]
    center = True
    rotate = True
    with open(f'{jobname}.colvars.conf') as conffile:
        for line in conffile:
            line = line.lstrip()
            if line.startswith('colvarstrajfrequency'):
                snapshotfreq = int(line.split()[1])
            elif line.startswith('centerreference'):
                if 'off' in line:
                    center = False
            elif line.startswith('rotatereference'):
                if 'off' in line:
                    rotate = False

    rlist = []
    trajlength = 0
    reps = 0
    for repnum in range(1, numreps + 1):
        if repnum == 1:
            repind = ''
        else:
            repind = f'-{repnum}'
        trajfilename = f'{jobname}{repind}.colvars.traj'
        if os.path.isfile(trajfilename):
            reps = repnum
            if not os.path.getsize(trajfilename):
                raise Exception(f'Empty trajectory file: {trajfilename}')
            with open(trajfilename) as trajfile:
                for line in trajfile:
                    if not line.startswith('#') and not line.split()[0] == '0':  # skip initial structure
                        trajlength += 1
                        line = line.split()
                        rlist.append((int(int(line[0]) / snapshotfreq), float(line[1])))
        else:
            if not reps:
                raise Exception('no completed runs')
            else:
                break
    rlist.sort(key=lambda x: x[1])

    ndim = 0
    with open(f'{jobname}.pdb') as pdbfile:
        for line in pdbfile:
            if line.startswith('ATOM'):
                if float(line[60:66]):
                    ndim += 3
    if center:
        ndim -= 3
        if rotate:
            ndim -= 3

    rtlnps = []
    for k in range(1, kmax + 1):
        # ln(pdf(r)) ~ ln(pdf_kNP(r)) - ln(k) + digamma(k + 1) if r is arbitrary point
        # pdf_kNP = (k / npoints) / volume
        # volume = (sqrt(pi) * radius) ^ n_dimensions / gamma(n_dimensions / 2 + 1)
        rtlnp = gas_constant * temp * (
                digamma(k + 1) - np.log(len(rlist))
                - ndim * (np.log(rlist[k - 1][1]) + np.log(np.pi) / 2)
                + lgamma(ndim / 2 + 1)
        )
        rtlnps.append(rtlnp)

    if convergence:
        resultline = ''
        if reps == numreps:
            for rtlnp in rtlnps:
                resultline += f'{rtlnp:<12.2f} '
        print(resultline)
    else:
        resultline = ''
        for rtlnp in rtlnps:
            resultline += f'{rtlnp:<12.2f} '
        if summary:
            print(resultline)
        else:
            print(f'{jobname}: {len(rlist)}/{trajlength} total')
            print('{:>24} {:>6} {:>24} {:>18}'.format(
                    'RT*ln(pdf) (kcal/mol)', 'k', 'kth NP distance (Å)', 'kth NP index'))
            for k in range(1, kmax + 1):
                print('{:>24.2f} {:>6} {:>24.2f} {:>18}'.format(rtlnps[k - 1], k, rlist[k - 1][1], rlist[k - 1][0]))
            print('Summary: ' + resultline)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='use kth nearest point approximation to calculate kT*ln(rho)')

    parser.add_argument('jobname', type=str, help='name of job (with or without extension)')
    parser.add_argument('-k', '--kmax', type=int, default=8, help='max value of k for which to compute RT*ln(p)')
    parser.add_argument('-s', '--summary', action='store_true', help='whether to only print final values')
    parser.add_argument('-n', '--numreps', type=int, default=1, help='max number of submissions to analyze')
    parser.add_argument('-c', '--convergence', action='store_true', help='test for convergence')
    parser.add_argument('-t', '--temperature', type=float, help='temperature at which simulations ran')

    args = parser.parse_args()

    main(args.jobname, args.temperature, args.kmax, args.summary, args.numreps, args.convergence)
