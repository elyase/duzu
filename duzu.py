"""
Python port of the Duflo and Zucker 10 parameter formula:

     http://amdc.in2p3.fr/web/dz.html


The MIT License

Copyright (c) 2013 Yaser Martinez Palenzuela

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from __future__ import division
import numpy as np
from math import sqrt

def duzu10(Z=82, N=126, coeffs=np.array([0.7043, 17.7418, 16.2562, 37.5562, 53.9017,
                                    0.4711, 2.1307, 0.0210, 40.5356, 6.0632], order='F')):
    """Calculates """
    term = np.zeros(10, order='F')
    onp = np.zeros((9, 2, 2), order='F')
    noc = np.zeros((18, 2), order='F')
    op, n2, dx, qx, os, oei, dei, pp, y = [np.zeros(2) for _ in range(9)]

    nuclei = [N, Z]
    a = sum(nuclei)
    t = abs(N - Z)
    r = a ** (1 / 3)
    rc = r * (1 - 0.25 * (t / a)**2)                            # Charge radius
    ra = rc**2 / r
    z2 = Z * (Z - 1)                                                              
    term[0] = (-z2 + 0.76 * z2**(2 / 3)) / rc                   # Coulomb energy

    for deformed in [0, 1]:                                     # deformed=0  spherical, deformed=1  deformed
        ju = 4 if deformed else 0                               # nucleons associated to deform.
        map(lambda x: x.fill(0), [term[1:], noc, onp, os, op])  # init with 0
        for I3 in [0, 1]:                                       # I3->isospin(0->N or 1->Z)
            n2[I3] = 2 * (nuclei[I3] // 2)                      # (for pairing calculation)
            ncum = i = 0
            while True:
                i += 1                                          # sub-shells (ssh) j and r filling
                idd = (i + 1) if i % 2 else (i * (i - 2) // 4)
                ncum += idd
                if ncum < nuclei[I3]:
                    noc[i-1, I3] = idd                           # nb of nucleons in each ssh
                else:
                    break
            imax = i + 1                                        # imax = last subshell nb
            ip = (i - 1) // 2                                   # HO number (p)
            ipm = i // 2
            pp[I3] = ip
            moc = nuclei[I3] - ncum + idd
            noc[i-1, I3] = moc - ju                             # nb of nucleons in last ssh
            noc[i, I3] = ju
            if i % 2:                                           # 'i' is odd, ssh j
                oei[I3] = moc + ip * (ip - 1)                   # nb of nucleons in last EI shell
                dei[I3] = ip * (ip + 1) + 2                     # size of the EI shell
            else:                                               # ssh r
                oei[I3] = moc - ju                              # nb of nucleons in last EI shell
                dei[I3] = (ip + 1) * (ip + 2) + 2               # size of the EI shell
            qx[I3] = oei[I3] * (dei[I3] - oei[I3] - ju) / dei[I3]    # n*(D-n)/D        S3(j)
            dx[I3] = qx[I3] * (2 * oei[I3] - dei[I3])           # n*(D-n)*(2n-D)/D  Q
            if deformed:
                qx[I3] /= sqrt(dei[I3])                         # scaling for deformed
            for i in range(1, imax+1):                          # Amplitudes
                ip = (i - 1) // 2
                fact = sqrt((ip + 1) * (ip + 2))
                onp[ip, 0, I3] += noc[i - 1, I3] / fact         # for FM term
                vm = -1.0
                if i % 2:
                    vm = 0.5 * ip                               # for spin-orbit term
                onp[ip, 1, I3] += noc[i - 1, I3] * vm
            for ip in range(0, ipm + 1):                        # FM and SO terms
                den = ((ip + 1) * (ip + 2))**(3. / 2)
                op[I3] = op[I3] + onp[ip, 0, I3]                # FM
                os[I3] = (os[I3] + onp[ip, 1, I3] * (1 + onp[ip, 0, I3]) * (ip * ip / den)
                               + onp[ip, 1, I3] * (1 - onp[ip, 0, I3]) * ((4 * ip - 5) / den))       # SO
            op[I3] = op[I3] * op[I3]

        term[1] = sum(op)                                       # Master term (FM): volume
        term[2] = -term[1] / ra                                 # surface
        term[1] += sum(os)                                      # FM + SO
        term[3] = -t * (t + 2) / r**2                           # isospin term : volume
        term[4] = -term[3] / ra                                 # surface
        if not deformed:                                        # spherical
            term[5] = sum(dx)                                   # S3  volume
            term[6] = -term[5] / ra                             # surface
            px = sqrt(pp[0]) + sqrt(pp[1])
            term[7] = qx[0] * qx[1] * (2 ** px)                 # QQ sph.
        else:                                                   # deformed
            term[8] = qx[0] * qx[1]                             # QQ deform.
        term[4] = t * (1 - t) / (a * ra**3) + term[4]           # Wigner term
        condition = (N > Z, n2[0] == nuclei[0], n2[1] == nuclei[1])     # PAIRING
        term[9] = {
            (True , True , False): 1 - t / a,
            (True , False, True ): 1,
            (False, False, False): 1,
            (True , False, False): t / a,
            (False, False, True ): 1 - t / a,
            (False, True , True ): 2 - t / a,
            (True , True , True ): 2 - t / a,
        }.get(condition, 0)
        term[1:] /= ra
        y[deformed] = np.dot(term, coeffs)
    binding_energy = y[0] if Z < 50 else max(y)                 # y[0]->spherical, y[1]->deformed
    mass_excess = Z * 7.28903 + N * 8.07138 - binding_energy
    return mass_excess

if __name__ == '__main__':
    print duzu10(100, 150)
