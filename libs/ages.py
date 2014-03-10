#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013, 2014 Junko Tsuji

# Estimate insertion ages with phylogenetic tree and
# genome alignments.

# built-in module
import string, os.path
from StringIO import StringIO

# modules in 'libs'
from io import readTree
from netaln import alignmentNets
from genaln import genomeAlignments


# search path to a specific node
def searchPath(goal, path, graph):
    if goal == path[-1]:
        return
    for key in graph[path[-1]]:
        if key not in path:
            path.append(key)
            searchPath(goal, path, graph)


# sort reference organelle names by ages
def sortKeys(keys, path, refTree):
    sorted = []
    keylen = len(keys)
    index = [i for i in range(keylen)]
    for i in range(keylen):
        if keys[i] in path:
            sorted.append(keys[i])
            index.pop(i)
            break
    if keylen > 2:
        for i in index:
            if keys[i] in refTree[sorted[0]]:
                sorted.append(keys[i])
                index.pop(i)
                break
    sorted.append(keys[index[0]])
    return sorted[::-1]


# define ages of reference organelle phylogenetic tree
def ageLabels(tree, outgroup, organelle):
    labels = {}
    cmplab = {}

    targetKeys, storedKeys = [], []
    species = 0
    for names in tree:
        species += len(names)
        if outgroup  in names: storedKeys.append(names)
        if organelle in names: targetKeys.append(names)
    alphabets = [chr(i) for i in range(97, 123) if chr(i) != 'o'][:species-1]
    for ch in alphabets: cmplab.setdefault(ch, [])
    viewed = set()
    cmplab.setdefault('o', []).append(outgroup)
    cmplab[alphabets[0]].append(organelle)

    labels.update({outgroup:'o'})
    labels.update({organelle:alphabets.pop(0)})
    
    path = [targetKeys[0]]
    searchPath(storedKeys[0], path, tree)
    if len(targetKeys[0]) > 1:
        for name in targetKeys[0]:
            if name not in labels:
                labels.update({name:alphabets[0]})
                cmplab[alphabets[0]].append("(" + name + ")")
    viewed.add(targetKeys[0])
    
    targetKeys = tree[targetKeys[0]]
    while True:
        storedKeys = []
        if len(targetKeys) > 1:
            targetKeys = sortKeys(targetKeys, path, tree)
        for key in targetKeys:
            for name in key:
                if name != outgroup:
                    cmplab[alphabets[0]].append(name)
                    labels.update({name:alphabets.pop(0)})
            viewed.add(key)
            for newKey in tree[key]:
                storedKeys.append(newKey)
        if len(labels.keys()) == species: break
        targetKeys = [key for key in storedKeys if key not in viewed]
    return labels, cmplab


# generate age definition string
def ageString(labels, treefile):
    blank = string.maketrans("", "")
    treeinfo = [l.translate(blank, "\n\t ") for l in open(treefile,"r")]
    index = list(labels.keys())
    index.sort()
    ageStr = StringIO()
    ageStr.write("#\n")
    ageStr.write("# Reference tree:\n")
    ageStr.write("# " + "".join(treeinfo) + "\n")
    ageStr.write("#\n")
    ageStr.write("# Lables:\n")
    for key in index:
        if labels[key]:
            ageStr.write("# " + key + ": ")
            ageStr.write(",".join(labels[key][::-1]) + "\n")
    ageStr.write("#\n")
    return ageStr.getvalue()


# define insertion age with phylogenetic tree
def defineAge(f, labels):
    age = "a"
    tree = readTree(f)
    targetKey = ()
    for key in tree:
        for name in key:
            if not labels.get(name):
                targetKey = key
                break
        if targetKey: break
    if len(targetKey) == 1:
        newTargetKey = []
        for key in tree[targetKey]:
            for name in key:
                newTargetKey.append(name)
        targetKey = newTargetKey
    refs = [labels[n] for n in targetKey if labels.get(n)]
    refs.append(age)
    age = max(refs)
    return age


# label insertion ages with phylogenetic trees
def nucAges(prog, reftree, treefile, files, org, out, opts, wd):
    if opts.verbose:
        print prog + ": label ages with phylgenetic tree"
    labels, cmplab = ageLabels(reftree, out, org)
    ageStr = ageString(cmplab, treefile)
    coordDict = {}
    for f in files:
        c = os.path.basename(f).split(".")[0].split(":")
        c.append(defineAge(f, labels))
        coordDict.setdefault(c[0], []).append(c[1:])

    if opts.alignment_net:
        if opts.verbose:
            print prog + ": verify ages with UCSC alignment nets"
        alignmentNets(prog, coordDict, labels, opts)

    if opts.genome_seqs:
        if opts.verbose:
            print prog + ": verify ages with genome alignments"
        ref = cmplab["a"][0]
        genomeAlignments(prog, coordDict, ref, labels, opts, wd)

    return coordDict, ageStr, cmplab.keys()
