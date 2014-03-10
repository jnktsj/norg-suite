#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013, 2014 Junko Tsuji

# Filter hits derived from nuclear segmental duplications.

# built-in modules
import commands, string

# modlues in 'libs'
from align import evalueSimulation
from io import readFasta, writeFasta, extractSeq


# generate flanking sequence coordinates
def makeFlankCoord(pos, flanklen):
    flankpos = []
    for p in pos:
        beg = str(max(p[0] - flanklen, 0)) + "," + str(p[1])
        end = str(p[0]) + "," + str(p[1] + flanklen)
        flankpos.append([beg, end])
    return flankpos


# extract flanking sequences
def extractFlankSeq(flankSeq, flankRev, flanklen, hitDict, nucleus):
    fs = open(flankSeq, "w")
    fr = open(flankRev, "w")
    wrote = 0
    mark = string.maketrans(":", "^")
    for name, seq in readFasta(nucleus):
        coord = makeFlankCoord(hitDict.get(name, []), flanklen)
        if not coord: continue
        for n, s in extractSeq(name, coord, seq, True):
            beg, end = n.split(":")[1:]
            n = ">" + name.translate(mark) + ":" + \
                beg.split(",")[-1] + ":" + end.split(",")[0]
            if len(s):
                writeFasta(n, s, fs)
                writeFasta(n, s[::-1], fr)
                wrote += 1
    fs.close()
    fr.close()
    return wrote


# format results and find nuclear duplicates
def formatResult(flankResult, coverage):
    nuclearDups = set()
    for line in flankResult.split("\n"):
        line = line.split("\t")
        if line[1] != line[6]:
            a = float(line[3])/float(line[5])
            b = float(line[8])/float(line[10])
            if max(a, b) >= coverage:
                nuclearDups.add(line[1])
                nuclearDups.add(line[6])
    return nuclearDups


# calculate overlapped fraction
def overlap(beg, end, hitBeg, hitEnd):
    hitlen = hitEnd - hitBeg
    if beg <= hitBeg and hitEnd <= end:
        return 1
    elif hitBeg < beg and end < hitEnd:
        return float(end - beg)/hitlen
    elif (hitBeg < beg and beg < hitEnd) and hitEnd <= end:
        return float(hitEnd - beg)/hitlen
    elif beg <= hitBeg and (hitBeg < end and end < hitEnd):
        return float(end - hitBeg)/hitlen
    return 0


# check nuclear duplicates with segmental duplication database
def checkSegDup(database, hitDict, coverage):
    segmentalDups = set()
    for line in open(database, "r"):
        line = line.rstrip("\n").split("\t")
        chrom = line[0]
        beg, end = map(int, line[1:3])
        for hit in hitDict.get(chrom, []):
            if overlap(beg, end, hit[0], hit[1]) >= coverage:
                title = chrom + ":" + ":".join(map(str, hit[:2]))
                segmentalDups.add(title)
    return segmentalDups


# filter duplicate hits
def filterDups(nuclearDups, hitDict):
    dupDict = {}
    mark = string.maketrans("^", ":")
    for line in nuclearDups:
        line = line.split(":")
        chrom = line[0].translate(mark)
        beg, end = map(int, line[1:3])

        i = 0
        hit = hitDict[chrom]
        hitlen = len(hit)
        while i < hitlen:
            if beg == hit[i][0] and end == hit[i][1]:
                dupDict.setdefault(chrom, []).append(hit[i])
                hitDict[chrom].pop(i)
                break
            i += 1
    # clean up empty arrays and sort existing entries
    delCard = []
    for chrom in hitDict:
        hitDict[chrom].sort()
        dupDict.get(chrom, []).sort()
        if not hitDict[chrom]:
            delCard.append(chrom)
    for chrom in delCard:
        del hitDict[chrom]
    return dupDict


# core funciton
def rmdup(prog, opts, wd, hitDict, organelle, nucleus):

    dupDict = {}
    dupCounts = 0
    nuclearDups = set()

    # command line format strings
    lastdb = opts.last + "/lastdb %s %s %s"
    lastal = opts.last + "/lastal -e%s -j4 -f0 %s %s | grep -v '#'"
    lastex = opts.last + "/lastex -E%s %s.prj %s.prj | sed -n 4p - | cut -f1"

    # flanking sequence file names
    flankSeq = wd + "/flankSeq.fa"
    flankRev = wd + "/flankRev.fa"

    # index file name
    flankIndex = wd + "/flankIndex"

    # compare flanking sequence similarity
    if opts.verbose:
        print prog + ": filter nuclear duplications"
        print prog + ": compare flanking sequence similarities"

    wrote = extractFlankSeq(flankSeq, flankRev,
                            opts.flank_len, hitDict, nucleus)
    if not wrote:
        hitCounts = sum(map(len, hitDict.values()))
        return hitCounts, dupDict, dupCounts
    
    commands.getoutput(lastdb % ("-c", flankIndex, flankSeq))
    score = commands.getoutput(lastex % ('1e-25', flankIndex, flankIndex))
    score, evalue = evalueSimulation(lastex, lastal, flankIndex,
                                     flankIndex, flankRev, score, '1e-25')
    flankResult = commands.getoutput(lastal % (score, flankIndex, flankSeq))
    nuclearDups = formatResult(flankResult, opts.dup_coverage)
    # crosscheck with segmental duplication database
    if opts.segdup_db:
        if opts.verbose:
            print prog + ": crosscheck with segmental duplication database"
        segmentalDups = checkSegDup(opts.segdup_db, hitDict, opts.dup_coverage)
        nuclearDups = nuclearDups.union(segmentalDups)

    # remove duplicate hits
    dupDict = filterDups(nuclearDups, hitDict)
    hitCounts = sum(map(len, hitDict.values()))
    dupCounts = sum(map(len, dupDict.values()))

    return hitCounts, dupDict, dupCounts
