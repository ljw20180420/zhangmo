# -*- coding: utf-8 -*-
"""
Created on Tue May 31 15:47:34 2016

@author: wxt
"""

from __future__ import division
import numpy as np

np.seterr(divide = "ignore")

def extract_matrix(clr, chrom, correct):

    if correct == False:
        M = clr.matrix(balance=correct, sparse=True).fetch(chrom).tocsr()
    else:
        M = clr.matrix(balance=correct, sparse=True).fetch(chrom).tocsr()
        with clr.open('r') as grp: 
            if 'scale' in grp['bins']['weight'].attrs: 
                scale = grp['bins']['weight'].attrs.get('scale')
            else:
                raw = clr.matrix(balance=False, sparse=True).fetch(chrom).tocsr()
                marg = np.array(raw.sum(0)).ravel()
                scale = marg[marg > 0].mean()
        M = M * scale
    
    return M


class Chrom(object):

    def __init__(self, chrom, res, hicdata):

        self.chrom = chrom
        self.res = res
        self.chromLen = hicdata.shape[0]
        self.hmm = None

        x, y = hicdata.nonzero()
        mask = x < y
        x, y = x[mask], y[mask]
        mat_ = hicdata[x, y]
        if isinstance(mat_, np.matrix):
            IF = np.array(mat_).ravel()
        else:
            # mat_ is a sparse matrix
            IF = np.array(mat_.todense()).ravel()

        IF[np.isnan(IF)] = 0
        self.IF =IF
        self.x, self.y = x, y

        del hicdata

    def calDI(self, window=2000000): # this is rewritten to behave the same as matlab version
        """
        Calculate DI for each bin.
        """
        ws = window // self.res
        # Perform filtering according to window size
        mask = self.y - self.x <= ws
        x, y = self.x[mask], self.y[mask]
        idata = self.IF[mask]
        
        Len = y.max() + 1
        # Downstream
        downs = np.bincount(x, weights = idata)
        # Upstream
        ups = np.bincount(y, weights = idata)
        ## Correct for length
        cdowns = np.zeros(Len)
        cdowns[:downs.size] = downs
        cups = np.zeros(Len)
        cups[(Len-ups.size):] = ups
        ## Formula
        numerators = cdowns - cups
        denominators = cdowns + cups
        self.DIs = numerators**3 / np.abs(numerators) / denominators
        self.DIs = np.nan_to_num(self.DIs, posinf=0, neginf=0)

    def pipe(self, seq, probs, minsize): # this is modified because the Map seems to be wrong and the regionStart is not necessary now
        """
        Estimate the median posterior probability of a region(a stretch of same
        state). We believe in a region only if it has a median posterior
        probability >= 0.99, or its size surpass 2 bins.
        
        TADs always begin with a single downstream biased state, and end with
        a last HMM upstream biased state.
        """
        path = [int(s.name) for i, s in self.hmm.viterbi(seq)[1][1:-1]]
        state_probs = self.hmm.predict_proba(seq)

        # determine 0 1 2 which state is left boundary, right boundary, no bias
        idx = np.argsort(np.array([self.hmm.states[i].parameters[0] for i in range(3)]))
        state_names = [int(self.hmm.states[id].name) for id in idx]

        # hmm.end.

        # Stretch consecutive same state  -->  Region
        mediate = []
        start = 0
        end = 1
        cs = path[0] # Current State
        prob_pool = [state_probs[0][cs]]
        for i in range(1, len(path)):
            state = path[i]
            if state != cs:
                mediate.append([start, end, cs, np.median(prob_pool)])
                start = i
                end = i + 1
                cs = state
                prob_pool = [state_probs[i][cs]]
            else:
                end = i + 1
                prob_pool.append(state_probs[i][cs])
        mediate.append([start, end, cs, np.median(prob_pool)])

        dawn = []
        # Calibrate the first and the last line
        if (mediate[0][1] - mediate[0][0]) <= 3:
            mediate[0][2] = mediate[1][2]
        if (mediate[-1][1] - mediate[-1][0]) <= 3:
            mediate[-1][2] = mediate[-2][2]
        
        dawn.append([mediate[0][0], mediate[0][1], mediate[0][2]])
        # Two criteria
        for i in range(1, len(mediate)-1):
            temp = mediate[i]
            if ((temp[1] - temp[0]) >= minsize) or (temp[-1] >= probs):
                dawn.append([temp[0], temp[1], temp[2]])
            else:
                Previous = mediate[i-1]
                Next = mediate[i+1]
                if Previous[2] == Next[2]:
                    dawn.append([temp[0], temp[1], Previous[2]])
                else:
                    dawn.append([temp[0], temp[1], 1])
        
        dawn.append([mediate[-1][0], mediate[-1][1], mediate[-1][2]])

        ## Infer TADs
        preTADs = []
        # Artificial Chromosome Size
        genome_size = dawn[-1][1]
        temp = []
        for i in range(len(dawn)):
            start = dawn[i][0]
            end = dawn[i][1]
            state = dawn[i][2]
            if i == 0:
                pre_state = state
                pre_end = end
                continue
            if state != pre_state:
                if pre_state == state_names[1]:
                    temp.append(start)
                if state == state_names[1]:
                    temp.extend([pre_end, pre_state])
                    preTADs.append(temp)
                    temp = []
                if (pre_state != state_names[1]) and (state != state_names[1]):
                    temp.extend([pre_end, pre_state])
                    preTADs.append(temp)
                    temp = [start]
            
            pre_end = end
            pre_state = state
        
        if pre_state != state_names[1]:
            temp.extend([genome_size, pre_state])
            preTADs.append(temp)
        
        TADs = []
        pre_state = -1
        temp = []
        for i in range(len(preTADs)):
            if pre_state == -1:
                if (preTADs[i][-1] != state_names[2]) or (len(preTADs[i]) < 3):
                    continue
                
            start = preTADs[i][0]
            end = preTADs[i][1]
            state = preTADs[i][2]
            
            if state != pre_state:
                if (state == state_names[2]) and (pre_state == -1):
                    temp.append(start)
                if (state == state_names[2]) and (pre_state == state_names[0]):
                    temp.append(pre_end)
                    TADs.append(temp)
                    temp = [start]
            
            pre_state = state
            pre_end = end
            
        if (pre_state == state_names[0]) and (len(temp) == 1):
            temp.append(pre_end)
            TADs.append(temp)

        return TADs

    
    def minCore(self, probs, minsize): # this is rewritten to abandon regionDIs, the whole DIs is used instead
        
        domains = self.pipe(self.DIs, probs, minsize)
        tmpDomains = []
        for domain in sorted(domains):
            domain[0] = domain[0] * self.res
            domain[1] = domain[1] * self.res
            tmpDomains.append(domain)

        return tmpDomains
    

    def callDomains(self, model, window=2000000, probs=0.99, minsize=3):
        
        self.hmm = model
        self.calDI(window=window)
        self.domains = self.minCore(minsize, probs, minsize)
    
    

