import sys
import threading, multiprocessing
from utils.helperfunctions import predictcrosspoints, find_runs, openfile
import numpy as np
from typing import Dict, Tuple, List, Set


def predictCrossoverRegions(kid, sireobj, damobj,maternal_strand=0, paternal_strand=1, 
                            min_flank_support=0, min_flank_support_fraction=0.01) -> Dict[int,Dict[int,Tuple]]:
    detectedCrossoverRegions: Dict[int,Dict[int,Tuple]] = {}
    for chrom in kid.data.keys():
        patcr, patcalls  = predictcrosspoints(
            kid.data[chrom][paternal_strand], #paternal strand
            sireobj.data[chrom][paternal_strand], #paternal paternal
            sireobj.data[chrom][maternal_strand], ignorevalue=9, 
            paternal_strand=paternal_strand,maternal_strand=maternal_strand) #paternal maternal
        runvalue, runstart, runlength = find_runs(patcr)
        supportForBlock = np.array([np.sum(patcalls[s:s+e] == v) for s,e,v in zip(runstart, runlength, runvalue)])
        supportForBlock_fraction = np.array([np.sum(patcalls[s:s+e] == v)/e for s,e,v in zip(runstart, runlength, runvalue)])
        
        exactTransition = [True if i > 0 and runvalue[i-1] < 3 and x < 3 else False for i,x in enumerate(runvalue)]
        
        xovers_i = np.where(np.logical_or(runvalue == 3,exactTransition))[0]
        xovers_i = xovers_i[np.logical_and(xovers_i > 0, xovers_i < len(runvalue)-1)]
        starts = runstart[xovers_i]
        xoverlength = runlength
        xoverlength[exactTransition] = 0
        lengths = xoverlength[xovers_i]
        before_length = supportForBlock[xovers_i-1]
        after_length = supportForBlock[xovers_i+1]
        before_fraction = supportForBlock_fraction[xovers_i-1]
        after_fraction = supportForBlock_fraction[xovers_i+1]
        if min_flank_support > 0:
            qualify = np.bitwise_and(np.bitwise_and(after_length >= min_flank_support, before_length >= min_flank_support),
                                     np.bitwise_and(before_fraction >= min_flank_support_fraction, after_fraction >= min_flank_support_fraction))
            patcregion = (starts[qualify], lengths[qualify])
        else:
            patcregion = (starts, lengths)
        
        matcr, patcalls = predictcrosspoints(
            kid.data[chrom][maternal_strand],#maternal strand        
            damobj.data[chrom][paternal_strand], #maternal paternal
            damobj.data[chrom][maternal_strand], ignorevalue=9,
            paternal_strand=paternal_strand,maternal_strand=maternal_strand) #maternal maternal
        runvalue, runstart, runlength = find_runs(matcr)
        supportForBlock = np.array([np.sum(patcalls[s:s+e] == v) for s,e,v in zip(runstart, runlength, runvalue)])
        supportForBlock_fraction = np.array([np.sum(patcalls[s:s+e] == v)/e for s,e,v in zip(runstart, runlength, runvalue)])
        
        exactTransition = [True if i > 0 and runvalue[i-1] < 3 and x < 3 else False for i,x in enumerate(runvalue)]
        xovers_i = np.where(np.logical_or(runvalue == 3,exactTransition))[0]
        xovers_i = xovers_i[np.logical_and(xovers_i > 0, xovers_i < len(runvalue)-1)]
        
        starts = runstart[xovers_i]
        xoverlength = runlength
        xoverlength[exactTransition] = 0
        lengths = xoverlength[xovers_i]
        before_length = supportForBlock[xovers_i-1]
        after_length = supportForBlock[xovers_i+1]
        before_fraction = supportForBlock_fraction[xovers_i-1]
        after_fraction = supportForBlock_fraction[xovers_i+1]
        
        if min_flank_support > 0:
            qualify = np.bitwise_and(np.bitwise_and(after_length >= min_flank_support, before_length >= min_flank_support),
                                     np.bitwise_and(before_fraction >= min_flank_support_fraction, after_fraction >= min_flank_support_fraction))
            matcregion = (starts[qualify], lengths[qualify])
        else:
            matcregion = (starts, lengths)
        detectedCrossoverRegions[chrom] = {maternal_strand:matcregion, paternal_strand:patcregion}
    return(detectedCrossoverRegions)

def writeCrossoverRegionsShortFormat(detectedregions:Dict[int,Dict[int,Tuple]], idkid, file, writemode="wt"):
    with openfile(file, writemode) as fout:
        for chrom, xovers in detectedregions.items():
            #print("ShortFormat %s %s" % (chrom, xovers))
            xoversstr = []
            if len(xovers[0][0]) > 0: #you can't do len on a zip so this is easier
                patcregion = list(zip(xovers[0][0], xovers[0][1])) #structure is ([starts], [lengths]) so we zip to start,length tuples
                xoversstr.extend(["p\t%s\t%s" % (r_starts,r_length) for r_starts,r_length in patcregion])
            if len(xovers[1][0]) > 0:
                matcregion = list(zip(xovers[1][0], xovers[1][1]))
                xoversstr.extend(["p\t%s\t%s" % (r_starts,r_length) for r_starts,r_length in matcregion])
            line = "%s\t%s\t%s\n" % (idkid, chrom, "\t".join(xoversstr))
            fout.write(line)

def writeCrossoverRegionsBedFormat(detectedregions:Dict[int,Dict[int,Tuple]], idkid, file, writemode="wt"):
    with openfile(file, writemode) as fout:
        for chrom, xovers in detectedregions.items():
            #print("BedFormat %s %s" % (chrom, xovers))
            if len(xovers[0][0]) > 0: #you can't do len on a zip so this is easier
                patcregion = list(zip(xovers[0][0], xovers[0][1])) #structure is ([starts], [lengths]) so we zip to start,length tuples
                for r_starts,r_length in patcregion:
                        fout.write("%s\t%s\t%s\tpat_%s\t%.3f\n" % (chrom, r_starts, r_starts+r_length, idkid, 1/r_length))
            if len(xovers[1][0]) > 0:
                matcregion = list(zip(xovers[1][0], xovers[1][1]))
                for r_starts,r_length in matcregion:    
                        fout.write("%s\t%s\t%s\tmat_%s\t%.3f\n" % (chrom, r_starts, r_starts+r_length, idkid, 1/r_length))


