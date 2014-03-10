#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013, 2014 Junko Tsuji

# Verify ages with UCSC alignment nets.

# built-in module
import string, re


# store homologous regions in net alignments
def storeNets(seq, pos, homolog, coord):
    L = len(coord)
    for l in range(L):
        length = 0
        beg = map(int, coord[l][0].split(","))
        end = map(int, coord[l][1].split(","))
        for i in range(len(beg)):
            length += float(end[i]-beg[i])
            homolog[l][pos] += sum(seq[beg[i]:end[i]])
        homolog[l][pos] /= length


# check net alignemnts whether there are homologous chunks
def checkNets(prog, coord, labels, nets, spes, ignore, opts):
    homolog = {}
    N = len(nets)
    pat = re.compile("^\s+", re.I)
    ignoreCard = ignore.split(",")

    for key in coord.keys():
        L = len(coord[key])
        homolog.setdefault(key, [[0] * N for i in range(L)])

    for i in range(N):
        f = False
        name, seq = "", []
        if opts.verbose:
            print prog + ": ...... surveying " + spes[i]
        for line in open(nets[i], "r"):
            if line.startswith("net"):
                l = line.rstrip("\n").split()
                if name and seq:
                    storeNets(seq, i, homolog[name], coord[name])
                if l[1] in coord:
                    f = True
                    name, seq = l[1], [0] * int(l[2])
                else:
                    f = False
                    name, seq = "", []
                continue
            if not f: continue
            headSpaces = pat.findall(line)
            if headSpaces:
                # level filter: hits should be above level 2
                level = len(headSpaces[0])
                if level > 2: continue
                l = line.lstrip()
                if l.split()[3] in ignoreCard: continue
                chunk = int(l.split()[2])
                beg = int(l.split()[1])
                end = chunk + beg
                if l.startswith("fill"):
                    seq[beg:end] = [1] * chunk
                elif l.startswith("gap"):
                    seq[beg:end] = [0] * chunk
        if name and seq:
            storeNets(seq, i, homolog[name], coord[name])
    return homolog


# revise ages if there are contradictions
def reviseAges(coord, labels, homolog, species):
    N = len(species)
    for chrom in coord:
        h = homolog[chrom]
        L = len(coord[chrom])
        for i in range(L):
            ages = [labels[species[j]] for j in range(N) if h[i][j] >= 0.1]
            if not ages: ages.append("a")
            else:        ages.append(coord[chrom][i][-1])
            coord[chrom][i][-1] = max(ages)


# verify ages with UCSC alignment nets
def alignmentNets(prog, coord, labels, opts):
    index = 0
    nets = []
    species = {}
    marks = string.maketrans(" \t()[]:;,", "________|")
    for line in open(opts.alignment_net, "r"):
        l = line.rstrip("\n").split()
        nets.append(l[1])
        species.setdefault(index, l[0].translate(marks))
        index += 1
    homolog = checkNets(prog, coord, labels, nets,
                        species, opts.ignore_card, opts)
    reviseAges(coord, labels, homolog, species)

