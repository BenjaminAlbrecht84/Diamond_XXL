**DIAMOND_XXL v0.8.3 by Benjamin Albrecht** - https://github.com/BenjaminAlbrecht84/Diamond_XXL




For running the program DIAMOND_XXL you have to perform the following two main steps:

1. Download and install DIAMOND (v0.8.17 or higher); see http://github.com/bbuchfink/diamond.
2. Get the runnable jar file ``diamond_xxl-0.8.5.jar`` from https://github.com/BenjaminAlbrecht84/Diamond_XXL/releases/download/v0.8.5/diamond_xxl.-0.8.5.jar

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

    $ java -Xmx 10g -jar diamond_xxl.jar -exe <diamond_binary> -d nr.dmnd -q reads.fna -o <output_folder>

where the ``-exe`` option specifies the path to the DIAMOND binary file and the ``-o`` option specifies the path to the output folder. Please get sure that you have installed DIAMOND v0.8.17 or higher, otherwise the program will not be able to read the ``nr.dmnd`` file. Moreover, please do not forget to ensure that the JVM gets enough memory. 

The output is automatically written into the output folder specified by the ``-a`` option. It consists of the two files ``reads.daa`` and ``reads.runs``.

- The binary file ``reads.daa`` containing all matches in the typical *DAA* format. Note that this file can be converted into other formats by using DIAMOND.
- The file ``reads.runs`` listing all runs of matches covering huge parts of protein sequences.

Options
=======

============== ======= ===========
Option         Default Description
============== ======= ===========
-exe                   Path to DIAMOND binary file.
-d                     Path to DIAMOND database file.
-q                     Path to query input file in FASTA or FASTQ format.
-o                     Path to output folder.
-p             max     Number of CPU threads.
--minBitScore  30      Minimum BitScore for a reported hit. 
--realign              Realings all alignments reported by DIAMOND.
-c             90      Minimum percentage of the reference that has to be covered by the query.
-m             10      Size of available main memory (in GB).
============== ======= ===========

Additionally, all DIAMOND-specific options can be set for configuring DIAMOND.

If no specific DIAMOND commands are defined, DIAMOND_XXL calls DIAMOND with its default parameters scoring an alignment with the *BLOSUM62* Matrix, a *gap open penalty* of 11 and a *gap extension penalty* of 1. 
