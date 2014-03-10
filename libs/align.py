#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013, 2014 Junko Tsuji

# Perform local alignments between organelle and nuclear genomes.


# built-in modules
import commands, os.path

# modules in 'libs'
from io import readFasta, writeFasta


# simulate suitable e-value
def evalueSimulation(lastex, lastal, d, f, s, score, evalue):
    falsePositives = 1
    while True:
        alignCommand = lastal % (score, d, s) + "| wc -l"
        falsePositives = int(commands.getoutput(alignCommand))
        if falsePositives == 0: break
        evalue = "%.2e" % (float(evalue) * 0.1)
        score = commands.getoutput(lastex % (evalue, d, f))
    return score, evalue


# correct doubled sequence positions
def doubleToSingleCoord(beg, end, length):
    netlen = int(length/2)
    if beg >= netlen: beg -= netlen
    if end >  netlen: end -= netlen
    if beg > end:
        return str(beg)+",0", str(netlen)+","+str(end)
    else:
        return str(beg), str(end)


# merge overlapped organelle coordinates
def correctOrgCoord(old, new, strand):
    pos = ",".join(map(str, old))
    if new:
        if strand == "-":
            pos = ",".join(map(str, new[::-1])) + "," + pos
        else:
            pos = pos + "," + ",".join(map(str, new))
    return pos


# merge overlapped nuclear coordinates
def mergeOverlap(lines):
    merged = []
    index = 1
    entries = len(lines)
    old = lines[0]
    while index < entries:
        new = lines[index]
        if old[1] < new[0] or old[5] != new[5]:
            merged.append(old)
            old = new
        else:
            if old[1] < new[1]: old[1] = new[1]
            oldBegs = map(int, old[3].split(","))
            oldEnds = map(int, old[4].split(","))
            newBegs = map(int, new[3].split(","))
            newEnds = map(int, new[4].split(","))
            for i in range(len(oldBegs)):
                j = 0
                while j < len(newBegs):
                    if not newBegs or not newEnds: break
                    if oldBegs[i] <= newBegs[j] and newEnds[j] <= oldEnds[i]:
                        newBegs.pop(j), newEnds.pop(j)
                        j -= 1
                        continue
                    if newBegs[j] <= oldBegs[i] and oldBegs[i] < newEnds[j]:
                        oldBegs[i] = newBegs[j]
                        if newEnds[j] > oldEnds[i]: oldEnds[i] = newEnds[j]
                        newBegs.pop(j), newEnds.pop(j)
                        j -= 1
                        continue
                    if newBegs[j] < oldEnds[i] and oldEnds[i] <= newEnds[j]:
                        oldEnds[i] = newEnds[j]
                        if newBegs[j] < oldBegs[i]: oldBegs[i] = newBegs[j]
                        newBegs.pop(j), newEnds.pop(j)
                        j -= 1
                        continue
                    j += 1
            old[3] = correctOrgCoord(oldBegs, newBegs, old[5])
            old[4] = correctOrgCoord(oldEnds, newEnds, old[5])
        index += 1
    merged.append(old)
    return merged


# format alignment result and create dictionary
def formatResult(alignResult):
    hitDict = {}
    for line in alignResult.split("\n"):
        tab = [''] * 7
        line = line.split("\t")
        orglen = int(line[5])
        tab[0], tab[3], tab[6] = line[6], line[1], line[9]
        if tab[6] == "-":
            tab[1] = int(line[10]) - int(line[7]) - int(line[8])
            tab[2] = int(line[10]) - int(line[7])
        else:
            tab[1] = int(line[7])
            tab[2] = int(line[7]) + int(line[8])
        tab[4] = int(line[2])
        tab[5] = int(line[2]) + int(line[3])
        tab[4], tab[5] = doubleToSingleCoord(tab[4], tab[5], orglen)
        hitDict.setdefault(tab[0], []).append(tab[1:])
    return hitDict


# core function
def align(prog, opts, args, wd):

    hitDict = {}
    hitCounts = 0

    # command line format strings
    tantan = opts.tantan + "tantan %s > %s"
    lastdb = opts.last   + "lastdb %s %s %s"
    lastal = opts.last   + "lastal -e%s -j4 -f0 %s %s | grep -v '#'"
    lastex = opts.last   + "lastex -E%s %s.prj %s.prj | sed -n 4p - | cut -f1"

    # (soft-masked) sequence file names
    organelle = wd + "/orgSeq"
    nucleus = wd + "/nucSeq"

    # index file names
    nucBaseFreq = wd + "/nucBaseFreq"
    orgIndex = wd + "/orgIndex"

    # simulation file names
    orgRev = wd + "/orgRev"
    orgRevIndex = wd + "/orgRevIndex"

    # double circular organelle genome
    if int(commands.getoutput("grep -c '>' " + args[0])) != 1:
        raise Exception("there must be exactly 1 sequence: " + args[0])
    commands.getoutput("cp " + args[0] + " " + organelle)
    if opts.dnaform == 'circular':
        if opts.verbose:
            print prog + ": double organelle sequence"
        commands.getoutput("grep -v '>' " + args[0] + " >> " + organelle)

    # soft-mask genomes with 'tantan'
    if opts.rmsk_organelle:
        if opts.verbose:
            print prog + ": soft-mask organelle genome"
        tmpfile = wd + "/tmp"
        commands.getoutput(tantan % (organelle, tmpfile))
        commands.getoutput("mv " + tmpfile + " " + organelle)
    if opts.rmsk_nucleus:
        if opts.verbose:
            print prog + ": soft-mask nuclear genome"
        commands.getoutput(tantan % (args[1], nucleus))
    else: nucleus = args[1]

    # index genomes with 'lastdb'
    if opts.rmsk_organelle or opts.rmsk_nucleus:
        lastdbOption = " "
    else:
        lastdbOption = "-c"
    commands.getoutput(lastdb % ("-x", nucBaseFreq, nucleus))
    commands.getoutput(lastdb % (lastdbOption, orgIndex, organelle))

    # calculate (and simulate) E-value
    if opts.evalue > 0: opts.greedy = False
    else              : opts.evalue = 0.01
    evalue = "%.2e" % opts.evalue
    score = commands.getoutput(lastex % (evalue, orgIndex, nucBaseFreq))
    if opts.greedy:
        if opts.verbose:
            print prog + ": simulate E-value threshold"
        f = open(orgRev, "w")
        for name, seq in readFasta(organelle):
            writeFasta(">" + name, seq[::-1], f)
        f.close()
        commands.getoutput(lastdb % (lastdbOption, orgRevIndex, orgRev))
        score, evalue = evalueSimulation(lastex, lastal, orgRevIndex,
                                         nucBaseFreq, nucleus, score, evalue)
        commands.getoutput("rm " + orgRev + " " + orgRevIndex + "*")

    # align organelle and nuclear genomes
    if opts.verbose:
        print prog + ": set: e-value=" + evalue + " / score=" + score
        print prog + ": start alignment"
    alignResult = commands.getoutput(lastal % (score, orgIndex, nucleus))
    if not alignResult:
        raise Exception("no hits are found")
    
    # format result
    hitDict = formatResult(alignResult)
    
    # merge overlapped hits
    if opts.verbose:
        print prog + ": cull overlapped alignments"
    for chrom in hitDict.keys():
        hitDict[chrom].sort()
        hitDict[chrom] = mergeOverlap(hitDict[chrom])
    hitCounts = sum(map(len, hitDict.values()))
    
    return hitDict, hitCounts, organelle, nucleus

