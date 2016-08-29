**DIAMOND_XXL v0.8.1 by Benjamin Albrecht** - https://github.com/BenjaminAlbrecht84/Diamond_XXL

DIAMOND_XXL is an extension of the BLAST-compatible local aligner DIAMOND for mapping **long** DNA query sequences (up to 10kb) against a protein database. In contrast to DIAMOND, DIAMOND_XLL only reports high scoring alignments of query sequences covering **whole** protein sequences. 

Usually, due to some typical properties of the query sequences, such as high numbers of insertions and deletions leading to *Frameshifts* within the corresponding protein sequence, the underlying algorithm of DIAMOND is not suitable for this purpose. However, in the beginning DIAMOND_XXL still makes use of DIAMOND for calculating an initial set of protein sequences, which is then further examined. Note that since DIAMOND is a pretty fast mapping tool (approx. 20,000 times faster than BLAST), this initial step, which is of high computational complexity, is performed in a very efficient way. 

Download & Installation
=======================

For running the program DIAMOND_XLL you have to perform the following two main steps:

1. Download and install DIAMOND (v0.8.17 or higher); see http://github.com/bbuchfink/diamond.
2. Get the runnable jar file ``diamond_xxl-0.8.2.jar`` from https://github.com/BenjaminAlbrecht84/Diamond_XXL/releases/download/v0.8.2/diamond_xxl-0.8.2.jar

Basic command line use
======================
We assume to have a protein database file in FASTA format named ``nr.faa`` and a file of DNA reads that we want to align named ``reads.fna``. In order to map the reads in ``reads.fna`` against the proteins in ``nr.faa`` you have to perform the following two steps:

1. Set up a reference database for DIAMOND_XXL
----------------------------------------------

In order to set up a reference database for DIAMOND_XXL, first DIAMOND has to be executed by running the ``makedb`` command as follows::

    $ diamond makedb --in nr.faa -d nr

This will create a DIAMOND database file with name ``nr.dmnd``. 

2. Aligning reads against the reference database
------------------------------------------------

The alignment task can then be initiated by running DIMAOND_XXL like this::

    $ java -Xmx 10g -jar diamond_xxl.jar -exe <diamond_binary> -d nr.dmnd -q reads.fna -a <output_folder>

where the ``-exe`` option specifies the path to the DIAMOND binary file and the ``-a`` option specifies the path to the output folder. Please get sure that you have installed DIAMOND v0.8.17 or higher, otherwise the program will not be able to read the ``nr.dmnd`` file. Moreover, please do not forget to ensure that the JVM gets enough memory. 

The output is automatically written into the output folder specified by the ``-a`` option. It consists of the two files ``reads.sam`` and ``reads.runs``.

- The file ``reads.sam`` listing all matches in typical *SAM* format.
- The file ``reads.runs`` listing all runs of matches covering whole protein sequences. Note that, since now sequences of *HSPs* are considered, *BitScores* are replaced by *SumScores* and *eValues* are replaced by *SumProbabilities*.

Options
=======

========== ======= ===========
Option     Default Description
========== ======= ===========
-exe               Path to DIAMOND binary file.
-d                 Path to DIAMOND database file.
-q                 Path to query input file in FASTA or FASTQ format.
-a                 Path to output folder.
-p         max     Number of CPU threads.
-l         1000    Length of the shredded reads. 
-s         30      Minimum sumScore to keep a run of matches.
-e         0.001   Maximum sumProbability to keep a run of matches. 
========== ======= ===========

DIAMOND_XXL initially calls DIAMOND with its default parameters scoring an alignment with the *BLOSUM62* Matrix, a *gap open penalty* of 11 and a *gap extension penalty* of 1. 
