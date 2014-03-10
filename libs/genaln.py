#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013, 2014 Junko Tsuji

# Verify ages with genome alignments

# built-in modules
import commands, string

# modlues in 'libs'
from align import evalueSimulation
from io import readFasta, writeFasta, extractSeq


# generate coordinates of hit and its peripheral regions
def editBaitCoord(pos, flanklen):
    baitPos = []
    for p in pos:
        b = int(p[0].split(",")[0])
        e = int(p[1].split(",")[-1])
        beg = ",".join([str(max(b - flanklen, 0)), p[0], str(e)])
        end = ",".join([str(b), p[1], str(e + flanklen)])
        baitPos.append([beg, end])
    return baitPos


# extract hit and its peripheral regions
def extractBaitSeq(baitSeq, baitRev, genome, coord):
    fs = open(baitSeq, "w")
    fr = open(baitRev, "w")
    for name, seq in readFasta(genome):
        c = editBaitCoord(coord.get(name, []), 200)
        if not c: continue
        for n, s in extractSeq(name, c, seq, True):
            beg, end = n.split(":")[1:]
            n = ">" + name + ":" + \
                ",".join(beg.split(",")[1:-1]) + ":" + \
                ",".join(end.split(",")[1:-1])
            writeFasta(n, s, fs)
            writeFasta(n, s[::-1], fr)
    fs.close()
    fr.close()


# format alignment results and verify ages
def formatResult(i, alignResult, homolog):
    for line in alignResult.split("\n"):
        line = line.split("\t")
        name, c = line[6], line[1]
        if line[9] == "-":
            beg = int(line[10]) - int(line[7]) - int(line[8])
            end = int(line[10]) - int(line[7])
        else:
            beg = int(line[7])
            end = int(line[7]) + int(line[8])
        unit = [0] * homolog[name][-1]
        homolog[name][i].setdefault(c, unit)
        for j in range(beg, end):
            homolog[name][i][c][j] = 1

    # reduce alignment information
    for name in homolog.keys():
        hitlist = []
        if not homolog[name][i].keys():
            homolog[name][i] = 0
            continue
        for c in homolog[name][i].keys():
            upper = sum(homolog[name][i][c][100:200])
            lower = sum(homolog[name][i][c][-200:-100])
            # flank sequences should aligned 40%
            if upper + lower > 180:
                hit = sum(homolog[name][i][c][200:-199])
                hitlist.append(hit)
        hitlen = float(max(hitlist))
        if hitlen/homolog[name][-1] > 0.1:
            homolog[name][i] = 1
        else:
            homolog[name][i] = 0


# align hit and the peripheral sequences
def checkGenomes(prog, coord, ref, labels, gens, spes, spesV, opts, wd):

    # command line format strings
    lastdb = opts.last + "/lastdb %s %s %s"
    lastal = opts.last + "/lastal -e%s -j4 -f0 %s %s | grep -v '#'"
    lastex = opts.last + "/lastex -E%s %s.prj %s.prj | sed -n 4p - | cut -f1"

    # index file name
    bBaseFreq = wd + "/bBaseFreq"

    # hit and its flanking sequences
    baitSeq = wd + "/baitSeq.fa"
    baitRev = wd + "/baitRev.fa"

    i = [k for k in spes.keys() if spes[k] == ref][0]
    extractBaitSeq(baitSeq, baitRev, gens[i], coord)
    commands.getoutput(lastdb % ("-x", bBaseFreq, baitSeq))
    spes.pop(i)

    L = spes.keys()
    L.sort()

    homolog = {}
    for chrom in coord.keys():
        for obj in coord[chrom]:
            key = ":".join([chrom, obj[0], obj[1]])
            beg = map(int, obj[0].split(","))
            end = map(int, obj[1].split(","))
            l = sum([end[i] - beg[i] for i in range(len(beg))])
            homolog.setdefault(key, [{} for j in L])
            homolog[key].append(l + 400)

    for j in range(len(L)):
        if opts.verbose:
            print prog + ": ...... surveying: " + spesV[L[j]]
        gIndex = wd + "/" + spes[L[j]] + "Index"
        commands.getoutput(lastdb % ("-c", gIndex, gens[L[j]]))
        score = commands.getoutput(lastex % ('1e-10', gIndex, bBaseFreq))
        score, e = evalueSimulation(lastdb, lastal, gIndex,
                                    bBaseFreq, baitRev, score, '1e-10')
        alignResult = commands.getoutput(lastal % (score, gIndex, baitSeq))
        formatResult(j, alignResult, homolog)

    return homolog


# revise ages if there are contradictions
def reviseAges(c, lab, hom, spes):
    L = spes.keys()
    L.sort()
    for chrom in c:
        N = len(c[chrom])
        for i in range(N):
            n = ":".join([chrom, c[chrom][i][0], c[chrom][i][1]])
            ages = [lab[spes[L[j]]] for j in range(len(L)) if hom[n][j]]
            if not ages: ages.append("a")
            else:        ages.append(c[chrom][i][-1])
            c[chrom][i][-1] = max(ages)


# verify ages with genome alignments
def genomeAlignments(prog, coord, ref, labels, opts, wd):
    # load data from file list
    index = 0
    genomes = []
    species, spVerbose = {}, {}
    marks = string.maketrans(" \t()[]:;,", "________|")
    for line in open(opts.genome_seqs, "r"):
        l = line.rstrip("\n").split()
        genomes.append(l[1])
        species.setdefault(index, l[0].translate(marks))
        spVerbose.setdefault(index, l[0])
        index += 1
    homolog = checkGenomes(prog, coord, ref, labels, genomes,
                           species, spVerbose, opts, wd)
    reviseAges(coord, labels, homolog, species)
