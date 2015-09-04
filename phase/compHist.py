#
# Filename: compHist.py
#
# Purpose: 
#   Compare two histograms to determine if they are significantly
#   different 



# Function: chi2GOF
#  Compares two histograms using the chi-squared goodness-of-fit 
#  test, assuming the first histogram passed in is the 'expected' 
#  and the second is the 'observed'

import scipy.special as ss
from math import sqrt
import scipy.stats as st


def chi2GOF (eHist, oHist):

    chi2 = 0.

    for i in range(0, len(eHist[:,0])):

        for j in range(0, len(eHist[0,:])):

            e = eHist[i,j]
            o = oHist[i,j]

#            print e, o

            if e!=0 and o!=0:
                chi2 += (o-e)*(o-e) / e

    return chi2


def chi2Matrix (pHist, qHist):
    
    chi2 = 0.

    for i in range(0, len(pHist[:,0])):

        for j in range(0,len(pHist[0,:])):

            p = pHist[i,j]
            q = qHist[i,j]

#            print p, q
            if q!=0 or p!=0:
                chi2 += (p-q)*(p-q) / (p+q)

    return chi2 / 2.0



# From Numerical Recipies
# Section 14.3
# Using Eqn 14.3.3
def chstwo (pHist, qHist):

    nbins = len(pHist[:,0])         # Number of bins
    df = nbins*nbins               # Degrees of freedom

    chi2 = 0.
    
    # Get the total of each hist
    pTot = 0.
    qTot = 0.
    for i in range(0,nbins):
        for j in range(0,nbins):
            pTot += pHist[i,j]
            qTot += qHist[i,j]

#    print pTot, qTot

    for i in range(0,nbins):
        for j in range(0,nbins):

            p = pHist[i,j]
            q = qHist[i,j]

            if p==0 and q==0:
                df -= 1
            else:
                top = sqrt(pTot/qTot)*q - sqrt(qTot/pTot)*p
                chi2 += top*top / (p+q)

    prob = ss.gammainc(0.5*df, 0.5*chi2)

    prob = st.chi2.cdf( chi2, df)
    return chi2, prob, df
