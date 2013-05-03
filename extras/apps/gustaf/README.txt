***  Gustaf - Generic mUlti-SpliT Alignment Finder ***
http://www.seqan.de/projects/gustaf.html

---------------------------------------------------------------------------
Table of Contents
---------------------------------------------------------------------------
  1.   Overview
  2.   Installation
  3.   Usage
  4.   Output Format
  5.   Example
  6.   Contact

---------------------------------------------------------------------------
1. Overview
---------------------------------------------------------------------------


Gustaf is a tool primarily designed for multi-split mapping of sequencing reads.

We present GUSTAF, a sound generic multi-split detection method implemented in
the C++ library SeqAn. GUSTAF uses SeqAn's exact local aligner Stellar to find
partial read alignments. Compatible partial alignments are identified, and a
split-read graph storing all compatibility information is constructed for each
read. Vertices in the graph represent partial alignments, edges represent
possible split positions. Using an exact dynamic programming approach, we
refine the alignments around possible split positions to determine precise
breakpoint locations at single-nucleotide level. We use a DAG shortest path
algorithm to determine the best combination of refined alignments, and report
those breakpoints supported by multiple reads.

---------------------------------------------------------------------------
2. Installation
---------------------------------------------------------------------------

Precompiled binaries (Linux 64-bit, Windows, Mac OS X) of Gustaf can be
downloaded from the SeqAn projects download page:
http://www.seqan.de/downloads/projects.html

Gustaf is distributed with SeqAn - The C++ Sequence Analysis Library (see
http://www.seqan.de). To build Gustaf yourself, you can check out the latest
SVN version of Gustaf and SeqAn with:

  1)  svn co http://svn.mi.fu-berlin.de/seqan/trunk/seqan
  2)  cd seqan/build/
  3)  cmake .. -DCMAKE_BUILD_TYPE=Release
  4)  make gustaf
  5)  ./bin/gustaf --help

If succesful, an executable file gustaf was built and a brief usage
description is dumped.

---------------------------------------------------------------------------
3. Usage
---------------------------------------------------------------------------

To get a short usage description of Gustaf, you can execute gustaf -h or
gustaf --help.

Usage: gustaf <GENOME FASTA FILE> <READ FASTA FILE> [Options]
       gustaf <GENOME FASTA FILE> <READ FASTA FILE> -m <GFF MATCH FILE> [Options]

  Using the first line, Gustaf will first run Stellar internally on the given
  input files. Please have a look at the Stellar options and default values for
  match length, error rate, etc.
  Using the --matchfile (-m) option, Gustaf skips running Stellar and
  directly uses the matches, preferable precalculated using Stellar, from
  the given GFF file.

---------------------------------------------------------------------------
3.1. Non-optional arguments
---------------------------------------------------------------------------

  Gustaf always expects both a database and a query file in Fasta format.
  All query sequences will be compared to all database sequences.
  Important: Following conventions, all Ids (line starting
  with '>') have to be unique already within the first part until the
  first whitespace! For example:
  >Gustaf|group=6|reads=9 readId=1
  >Gustaf|group=6|reads=9 readId=2
  are not unique.

  Without any additional parameters, Gustaf would call Stellar and then chain
  those matches of each read, that have either an overlap of at least 0.5 (50%
  of each match length) or a distance between the matches of at most 10 bp,
  and that miss at most 15 bp at the end or beginning of the read.
  The program calls Stellar with default options: Stellar would
  compare the query sequence(s) to both strands of the database sequence(s)
  with an error rate of 0.05 (i.e. 5% errors, an identity of 95%) and a
  minimal length of 100, and dump all local alignments in an output file named
  "stellar.gff" (see Stellar Readme).
  
  Two output files will be generated: breakpoint.gff includes all structural
  variant breakpoints in GFF format, indel.gff includes small indel (usually
  <10bp) in GFF format.

  The default behaviour can be modified by adding the following options to
  the command line.

---------------------------------------------------------------------------
3.2. Main Options
---------------------------------------------------------------------------

---------------------------------------------------------------------------
3.2.1. Gustaf-Specific Main Options
---------------------------------------------------------------------------

  [ -m FILE ],  [ --matchfile FILE ]

  Set the name of a file containing matches in GFF format.

  When setting option --matchfile (-m), Gustaf uses the matches from the
  given input file instead of calling Stellar (practical for running multiple
  tests on the same data). See also sample file stellar.gff.

  [ -oth REAL ],  [ --overlapThresh REAL ]

  Required overlap for matches of one read to be considered for chaining
  (main criterion, default 0.5).

  [ -gth NUM ],  [ --gapThresh NUM ]

  Maximal allowed distance between two matches to be considered for
  chaining. This criterion only applies if the matches do not overlap
  (default 10).

  [ -ith NUM], [ --initGapThresh NUM ]

  Maximal allowed length of leading or ending gap for the whole read(!)
  (default 15 for each).

  [ -st NUM], [ --support NUM ]

  Required number of supporting reads (or contigs) for a breakpoint to be
  written to the output file (see breakpoint output, default 2).

  [ -tp NUM], [ --transPen NUM ]

  Interchromosomal translocation penalty (default 5). Added to edge weight
  between matches that map to different database sequences.

  [ -ip NUM], [ --invPen NUM ]

  Inversion penalty (default 5). Added to edge weight between matches that
  map to different strands of the same database sequence. Note: Penalty is
  only added if there is no transPen already on the edge.

  [ -op NUM], [ --orderPen NUM ]

  Order penalty (default 5). Added to edge weight between matches that
  map to a database sequence in a different order than they are in the read
  sequence. Note: Penalty is  only added if there is no transPen  or invPen 
  already on the edge.

---------------------------------------------------------------------------
3.2.2. Stellar-Inherited Options
---------------------------------------------------------------------------

  This is only a tiny snapshot showing those Stellars options most crucial
  to get started using Gustaf. For more Stellar options, please see Stellar
  documentation and README.

  [ -e NUM ],  [ --epsilon NUM ]

  Set the maximal error rate for local alignments. NUM must be a floating
  point number between 0 and 0.25. The default value is 0.05. NUM is the
  number of edit operations needed to transform the aligned query substring
  to the aligned database substring divided by the length of the local
  alignment. For example specify '-e 0.1' for a maximal error rate of 10%
  (90% identity).

  [ -l NUM ],  [ --minLength NUM ]

  Set the minimal length of a local alignment. The default value is 100.
  NUM is the minimal length of an alignment, not the length of the
  aligned substrings.

  [ -f ],  [ --forward ]

  Only compare query sequence(s) to the positive/forward strand of the
  database sequence(s). By default, both strands are scanned.

  [ -r ],  [ --reverse ]

  Only compare query sequence(s) to the negative/reverse-complemented
  database sequence(s). By default, both strands are scanned.

---------------------------------------------------------------------------
3.3 Output Options
---------------------------------------------------------------------------

  [ -j ],  [ --jobName ]

  Optional jobname that will be added to all dot output files in the format
  read#_jobname.dot. # stands for the sequence number of the read in the
  query input file.

  [ -do ], [ --dots ]

  Switches output of DOT files on and off (default off). Each DOT file
  includes the graph representation of one read/query.

---------------------------------------------------------------------------
4. Output Formats
---------------------------------------------------------------------------

Gustaf currently supports the GFF output format for reporting breakpoints.

---------------------------------------------------------------------------
4.1. General Feature Format (GFF)
---------------------------------------------------------------------------

The General Feature Format is specified by the Sanger Institute as a tab-
delimited text format with the following columns:

<seqname> <src> <feat> <start> <end> <score> <strand> <frame> [attr] [cmts]

See also: http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml
Consistent with this specification Gustaf GFF output looks as follows:

GNAME Gustaf SVTYPE GBEGIN GEND . GSTRAND . ATTRIBUTES
StartSeqId      Label   SV type sPos    sPos    .       +/-     .
Match value description:

  GNAME        Name of the genome sequence (see --genome-naming)
  GUSTAF       Constant
  SVTYPE       SV type (deletion, insertion, inversion, translocation)
  GBEGIN       Beginning position in the genome sequence
               (positions are counted from 1)
  GEND         End position in the genome sequence (included!)
  .            Constant
  GSTRAND      '+'=forward strand or '-'=reverse strand
  .            Constant
  ATTRIBUTES   A list of attributes in the format <tag_name>[=<tag>]
               separated by ';'

Attributes are:

  ID=          Random number as identifier
  size=        Indel size (deletion and insertion only), "-" before the
               size indicates an insertion
  seq=         Sequence content of insertion (obviously insertion only).
  endChr=      Name of the second genome sequence (inversion and
               interchromosomal translocation only)
  endPos=      First position in second genome sequence (inversion and
               interchromosomal translocation only)
  endStrand=   '+'=forward strand or '-'=reverse strand of second genome
               sequence (inversion and interchromosomal translocation only)
  support=     Number of supporting reads for this breakpoint
  supportIds=  IDs of supporting reads
  breakpoint:  Breakpoint position within read

For matches on the reverse strand, GBEGIN and GEND are positions on the
related forward strand. It holds GBEGIN < GEND, regardless of GSTRAND.

---------------------------------------------------------------------------
5. Example
---------------------------------------------------------------------------

See the 'tests' folder for example runs of Gustaf. There is a genome file
'adeno.fa' and a modified genome file 'adeno_modified.fa', from the latter
we created the read file 'adeno_modified_reads.fa'. There is a pre-calculated
Stellar output file 'stellar.gff' computed by calling

./stellar adeno.fa adeno_modified_reads.fa -l 30 -o stellar.gff

The default calls for Gustaf would then be

./gustaf adeno.fa adeno_modified_reads.fa -st 1 -l 30 -bpo st1_l30.gff
or
./gustaf adeno.fa adeno_modified_reads.fa -st 1 -m stellar.gff -bpo st1_l30_m.gff

Both calls produce an output file containing the same breakpoints.
In the first run, Gustaf internally calls Stellar with parameter -l 30.
In the second run, Gustaf used the pre-calculated file with Stellar matches
(use this option for larger datasets where you want to run Stellar separately
and reuse the Stellar output for multiple Gustaf runs).

---------------------------------------------------------------------------
6. Contact
---------------------------------------------------------------------------

For questions or comments, contact:
  Kathrin Trappe  <kathrin.trappe@fu-berlin.de>
