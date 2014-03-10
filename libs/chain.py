#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013, 2014 Junko Tsuji

# Chain detected hits with indel, chaining, and overlap threshoulds.


# merge organelle coordinates
def mergeOrg(d, hitDict):
    unchange = 0
    for chrom in hitDict:
        intact = hitDict[chrom][:]
        for hit in hitDict[chrom]:
            hitBeg = map(int, hit[3].split(","))
            hitEnd = map(int, hit[4].split(","))
            L, i = len(hitBeg), 1
            if L == 1: continue
            while i < L:
                span = hitBeg[i] - hitEnd[i-1]
                if abs(span) <= d:
                    hitBeg.pop(i), hitEnd.pop(i-1)
                    L = len(hitBeg)
                    i -= 1
                if L == 1: break
                i += 1
            hit[3] = ",".join(map(str, hitBeg))
            hit[4] = ",".join(map(str, hitEnd))
        if intact == hitDict[chrom]:
            unchange += 1
    return unchange


# chain organelle coordinates
def chainOrg(oldBeg, oldEnd, newBeg, newEnd, strand):
    if strand == "-":
        beg, end = newBeg, oldEnd
        a = oldBeg.split(",")[1:]
        b = newEnd.split(",")[:-1]
        if a: beg += "," + ",".join(oldBeg.split(",")[1:])
        if b: end  = ",".join(newEnd.split(",")[:-1]) + "," + oldEnd
    else:
        beg, end = oldBeg, newEnd
        a = newBeg.split(",")[1:]
        b = oldEnd.split(",")[:-1]
        if a: beg += "," + ",".join(newBeg.split(",")[1:])
        if b: end  = ",".join(oldEnd.split(",")[:-1]) + "," + newEnd
    return beg, end


# calculate span between front and rear organelle hits
def calcOrgSpan(oldBeg, oldEnd, newBeg, newEnd, strand):
    if strand == "-":
        orgSpan = int(oldBeg.split(",")[0]) - int(newEnd.split(",")[-1])
    else:
        orgSpan = int(newBeg.split(",")[0]) - int(oldEnd.split(",")[-1])
    return abs(orgSpan)


# chain and merge nuclear coordinates
def mergeNuc(d, i, c, o, hitDict):
    unchange = 0
    for chrom in hitDict:
        hitlen = len(hitDict[chrom])

        index = 1
        chained = []
        old = hitDict[chrom][0]
        while index < hitlen:
            new = hitDict[chrom][index]
            if old[5] != new[5]:
                chained.append(old)
                old = new
                index += 1
                continue
            nucSpan = int(new[0].split(",")[0]) - int(old[1].split(",")[-1])
            if nucSpan > i:
                chained.append(old)
                old = new
                index += 1
                continue

            orgSpan = calcOrgSpan(old[3], old[4], new[3], new[4], old[5])

            # deletion and chaining
            if nucSpan <= (c + d):
                old[1] = ",".join(old[1].split(",")[:-1] + [new[1]])
                if orgSpan <= c:
                    old[3], old[4] = chainOrg(old[3], old[4], new[3], new[4], old[5])
                else:
                    if old[5] == "-":
                        old[3] = new[3] + "," + old[3] 
                        old[4] = new[4] + "," + old[4]
                    else:
                        old[3] = old[3] + "," + new[3]
                        old[4] = old[4] + "," + new[4]
            # insertion and overlap
            elif nucSpan <= i and orgSpan <= o:
                old[0] += "," + new[0]
                old[1] += "," + new[1]
                old[3], old[4] = chainOrg(old[3], old[4], new[3], new[4], old[5])
            # other
            else:
                chained.append(old)
                old = new
            index += 1

        chained.append(old)
        if hitDict[chrom] == chained: unchange += 1
        else: hitDict[chrom] = chained
    return unchange


# core function
def chain(opts, hitDict):
    entries = len(hitDict.keys())

    # merge organelle and nuclear coordinates
    sameOrg = 0
    sameNuc = 0
    while sameOrg + sameNuc != entries * 2:
        if sameOrg != entries:
            sameOrg = mergeOrg(opts.deletion, hitDict)
        if sameNuc != entries:
            sameOrg = 0
            sameNuc = mergeNuc(opts.deletion, opts.insertion,
                               opts.concat, opts.overlap, hitDict)
    hitCounts = sum(map(len, hitDict.values()))
    return hitCounts
