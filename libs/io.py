#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013, 2014 Junko Tsuji

# Input and Output functions.

import sys, os.path, string, re
from colorsys import hsv_to_rgb

# check whether path exists or not
def pathExists(path):
    if not os.path.exists(path):
        raise Exception("can't find "+path)


# read fasta file(s)
def readFasta(file):
    name = ""
    for line in open(file, "r"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    if name: yield name, "".join(seq)


# write fasta file
def writeFasta(name, seq, f=sys.stdout):
    linesize = 50
    beg = 0
    seqlen = len(seq)
    f.write(name + "\n")
    while beg < seqlen:
        end = beg + linesize
        f.write(seq[beg:end] + "\n")
        beg = end


# function to generate RGB color
def rgbColors(num):
    colors = []
    unit = 1/float(num)
    for i in range(num):
        t = tuple(int(j * 255) for j in hsv_to_rgb(unit * i, 1, 1))
        colors.append('%s,%s,%s' % t)
    return colors


# write coordinate file
def writeCoord(hitDict, file, s=None, ages=None):
    f = open(file, "w")
    chroms = hitDict.keys()
    chroms.sort()

    # write definition and bed strings
    palette = {}
    prefix = os.path.basename(file).split(".")[0]
    bedString = 'track name="' + prefix + '" visibility=2 itemRgb="On"\n'
    if s:
        f.write(s + bedString)
        flag = '1'
        ages.sort()
        palette = dict(zip(ages, rgbColors(len(ages))))
    else:
        f.write(bedString)
        flag = '0'

    # output coordinate in bed format
    for chrom in chroms:
        for hit in hitDict[chrom]:
            b = map(int, hit[0].split(","))
            e = map(int, hit[1].split(","))
            if s: name = ":".join(hit[2:-2]) + ":" + hit[-1]
            else: name = ":".join(hit[2:-1])
            color = palette.get(hit[-1], "0,0,0")
            count = str(len(b))
            sizes = ",".join([str(e[i] - b[i]) for i in range(len(b))])
            starts = ",".join([str(b[i] -b[0]) for i in range(len(b))])
            line = [chrom, str(b[0]), str(e[-1]), name, flag, hit[5],
                    hit[0], hit[1], color, count, sizes, starts]
            f.write("\t".join(line) + "\n")
    f.close()


# extract sequences with coordinates
def extractSeq(chrom, pos, seq, useMark=False):
    mark = string.maketrans(":", "_")
    for p in pos:
        if useMark: c = chrom.translate(mark)
        else:       c = chrom
        name = c + ":" + ":".join(p)
        subseq = ""
        beg = map(int, p[0].split(","))
        end = map(int, p[1].split(","))
        chunk = len(beg)
        for i in range(chunk):
            subseq += seq[beg[i]:end[i]]
        yield name, subseq


# read newick files
def readTree(file):
    blank = string.maketrans("", "")
    upper = re.compile("\(")
    lower = re.compile("\)")

    # load newick file as string object
    treeStr = ""
    for line in open(file, "r"):
        line = line.translate(blank, "\t\n ")
        treeStr += line
    if treeStr[-1] != ";":
        raise Exception(file + " should be end with ';'")
    if len(upper.findall(treeStr)) != len(lower.findall(treeStr)):
        raise Exception("tree format error: " + file)

    tree = {}
    treeElem = treeStr.split(",")
    stack, bifur, child, parent = [], [], [], []
    
    for e in treeElem:
        U = upper.findall(e)
        if U: stack.extend(U)
        stack.append(e.split(")")[0].replace("(",""))
        L = lower.findall(e)
        if not L: continue

        for l in L:
            while True:
                tmp = stack.pop()
                if tmp == '(': break
                child.append(tmp)
            stacklen = len(stack) - 1
            for i in range(stacklen, 0, -1):
                if stack[i] == '(': break
                parent.append(stack[i])
            if child and not parent:
                bifur.append(tuple(child))
            elif not child and parent:
                bifur.append(tuple(parent))
                bifurlen = len(bifur)
                for i in range(bifurlen):
                    for j in range(bifurlen):
                        if i != j:
                            tree.setdefault(bifur[i], []).append(bifur[j])
            else:
                tree.setdefault(tuple(child), []).append(tuple(parent))
                tree.setdefault(tuple(parent), []).append(tuple(child))
            child, parent = [], []
    return tree
