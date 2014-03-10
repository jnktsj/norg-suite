norg-suite
==========

This package contains two programs:

* norg-seq: finds organelle genome insertions in the nuclear
            genomes with local alignments between the organelle
            and the nuclear genomes.

* norg-age: estimates evolutionary ages of organelle insertions
            that have been defined by norg-seq, using phylogenetic
            trees and genome alignments.


Requirement
-----------

norg-suite requires the following environment to run:

1) Python 2.5, or more latest

2) Third party tools

   * norg-seq needs:
      - tantan (repeat masking): http://www.cbrc.jp/tantan/
      - LAST (local alignments): http://last.cbrc.jp/

   * norg-age needs:
      - MAFFT (multiple alignments):
               http://mafft.cbrc.jp/alignment/software/
      - RAXML (phylogenetic trees) :
               http://sco.h-its.org/exelixis/software.html

To prepare the above tools easily, there is a shell script
which downloads the newest versions of the tools automatically.
Under "tools" directory, you can use the script like this:

   chmod +x fetch.sh && ./fetch.sh -a

After fetching the tools, please compile them by following
each software instruction.


Usage
-----

Prease see "doc" directory for the details. There is two types
of manual formats, html and txt. We recommend to view html
version for getting more detailed descriptions of norg-suite.

Manuals enclosed "html" or "text" named:

   * norg-suite general description: norg-suite.html   (html only)
   * norg-seq user's manual        : norg-seq.html|txt (html, txt)
   * norg-age user's manual        : norg-age.html|txt (html, txt)

