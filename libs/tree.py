#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2012, 2013, 2014 Junko Tsuji

# Build phylogenetic trees of organelles and nuclear insertions.

# built-in module
import commands, string, os.path
from time import gmtime, strftime

# modules in 'libs'
from io import readFasta, writeFasta


# execute mafft
def exeMafft(path, a, b, c, l):
    if l < 5000000: p = "--maxiterate 1000"
    else:           p = "--maxiterate 2"
    mafft = path + "mafft --quiet --preservecase %s %s %s > %s"
    commands.getoutput(mafft % (p, a, b, c))
    if not os.path.exists(c):
        raise Exception("mafft calculation failed")


# execute raxml
def exeRaxml(path, a, o, n, wd, guide=None, target=None):
    raxml1 = path + \
        "raxmlHPC -f a -m GTRGAMMA -x 12345 -# 500 -n T1" + n + " -s "
    raxml2 = path + \
        "raxmlHPC -m GTRGAMMA -J MRE -z RAxML_bootstrap.T1" + n + " -n T2" + n
    commands.getoutput(raxml1 + a)
    commands.getoutput(raxml2)
    commands.getoutput("mv RAxML_MajorityRuleExtendedConsensusTree.T2" + n
                       + " " + o)
    if not os.path.exists(o):
        if guide:
            ob = os.path.basename(o).split(".")[0].replace(":","_")
            pat = string.maketrans("", "")
            tree = "".join([l.translate(pat, "\t\n ") for l in open(guide)])
            tree = tree.replace(target, target+","+ob)
            f = open(o, "w")
            f.write(tree+"\n")
            f.close()
        else:
            raise Exception("raxml computation failed")
    commands.getoutput("mv RAxML_*" + n + "* " + wd)


# write phylip file with fasta file
def fastaToPhylip(infile, outfile):
    f = open(outfile, "w")
    marks = string.maketrans(" \t()[]:;,", "________|")
    orgs = [[name, seq] for name, seq in readFasta(infile)]
    maxlen = max(map(len, zip(*orgs)[0]))
    seqlen = len(orgs[0][1])
    f.write(str(len(orgs)) + " " + str(seqlen) + "\n")
    for o in orgs:
        space = " " + " " * (maxlen - len(o[0]) + 3)
        f.write(o[0].translate(marks) + space + o[1] + "\n")
    f.close()


# find index of aligned reference organelle genome
def findRefIndex(beg, end, seq):
    refBeg, refEnd = [], []
    chunk = len(beg)
    for i in range(chunk):
        b = len(seq[:beg[i]].replace("-", ""))
        e = len(seq[:end[i]].replace("-", ""))
        begOffset, endOffset = beg[i], end[i]
        while b < beg[i]:
            if seq[begOffset] != "-": b += 1
            begOffset += 1
        while e < end[i]:
            if seq[endOffset] != "-": e += 1
            endOffset += 1
        refBeg.append(begOffset)
        refEnd.append(endOffset)
    return refBeg, refEnd


# compute length of first sequence in multi-fasta file
def seqLen(file):
    end = False
    seq = ""
    for line in open(file, "r"):
        if line.startswith(">"):
            if end: return len(seq)
            end = True
        else:
            seq += line.rstrip("\n")


# core function for reference organelle genomes
def orgTree(prog, opts, seq, wd):
    refFasta = opts.ref_msa
    refPhylip = wd + "/refPhylip"
    refTree = opts.ref_tree

    # calculate multiple alignment file
    if not refFasta:
        refFasta = wd + "/refFasta"
        if opts.verbose:
            print prog + ": calculate organelle multiple alignment"
        exeMafft(opts.mafft, "--globalpair", seq, refFasta, seqLen(seq))
    fastaToPhylip(refFasta, refPhylip)

    # compute phylogenetic tree
    if not refTree:
        refTree = wd + "/refTree"
        if opts.verbose:
            print prog + ": compute organelle phylogenetic tree"
        id = strftime("%H%M%S_org", gmtime())
        exeRaxml(opts.raxml, refPhylip, refTree, id, wd)
    return refFasta, refTree


# core function for organelle insertions
def nucTree(prog, opts, f, refSeq, wd, reftree, organelle):

    orgGenomes = len(refSeq)

    files = []
    id = strftime("%H%M%S_nuc", gmtime())
    revcmp = string.maketrans("ACGTNSWRYKMBDHVacgtnswrykmbdhv",
                              "TGCANSWYRMKVHDBtgcanswyrmkvhdb")
    # generate multi-fasta sequence files
    if opts.verbose:
        print prog + ": create fasta file(s): insertions and corresponding organelle subsequences"

    for name, seq in readFasta(f):
        coord = name.split(":")
        beg = map(int, coord[4].split(","))
        end = map(int, coord[5].split(","))
        refBeg, refEnd = findRefIndex(beg, end, refSeq[0][1])
        if coord[-1] == "-":
            seq = seq[::-1].translate(revcmp)
        file = wd + "/" + name + ".mfa"
        fp = open(file, "w")
        chunk = len(refBeg)
        writeFasta(">" + name, seq, fp)
        for i in range(orgGenomes):
            subseq = ""
            for j in range(chunk):
                subseq += refSeq[i][1][refBeg[j]:refEnd[j]]
            if subseq:
                writeFasta(">" + refSeq[i][0], subseq.replace("-", ""), fp)
        fp.close()
        files.append(file)

    # compute multiple sequence alignments
    flen = len(files)
    tmp = wd + "/" + id
    if opts.verbose:
        print prog + ": compute multiple sequence alignments of insertion sequences"
    for i in range(flen):
        msa = files[i].split(".")[0] + ".msa"
        exeMafft(opts.mafft, "--localpair", files[i], tmp, seqLen(files[i]))
        fastaToPhylip(tmp, msa)
        files[i] = msa
    commands.getoutput("rm " + tmp)

    # calculate phylogenetic trees
    if opts.verbose:
        print prog + ": calculate phylogenetic trees of insertion sequences"
    for i in range(flen):
        ptr = files[i].split(".")[0] + ".ptr"
        exeRaxml(opts.raxml, files[i], ptr, id, wd, reftree, organelle)
        files[i] = ptr
    return files

