This file is from:

    http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP17way/README.txt

This directory contains conservation scoring by phyloP (phylogenetic p-values)
from the PHAST package (http://compgen.cshl.edu/phast/) for multiple
alignments of 16 primate genomes to the human genome.

Files in this directory:

    hg38.phyloP17way.mod - phyloP tree model with branch lengths for
                            all 17 species
    hg38.phyloP17way.bw - scores in bigWig format
    hg38.phyloP17way.wigFix.gz - scores in fixed step wiggle format
    md5sum.txt - MD5 checksums to verify download files

See also: http://genome.ucsc.edu/FAQ/FAQformat.html#format6.1
   for bigWig format description

For a description of the fixed step wiggle data file format, see:
http://genome.ucsc.edu/goldenPath/help/phastCons.html

The phastCons data can be found at:
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons17way/

The multiple alignments and methods description are at:
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/multiz17way/

For more information about this data, see the track
description for the Conservation track:
    http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons17way

Assemblies used in these alignments:                                 (alignment
                                                                           type)

Human - Homo sapiens              Dec. 2013 (GRCh38/hg38)            reference

Baboon - Papio anubis             Mar. 2012 (Baylor Panu_2.0/papAnu2) syntenic
Bushbaby - Otolemur garnettii       Mar. 2011 (Broad/otoGar3) reciprocal best
Bonobo - Pan paniscus             May. 2012 (Max-Planck/panPan1) reciprocal best
Chimp - Pan troglodytes           Feb. 2011 (CSAC 2.1.4/panTro4)      syntenic
Crab-eating macaque - Macaca fascicularis
                          Jun 2013 (Macaca_fascicularis_5.0/macFas5)  syntenic
Gibbon - Nomascus leucogenys      Oct. 2012 (GGSC Nleu3.0/nomLeu3)    syntenic
Golden snub-nosed monkey - Rhinopithecus roxellana
                          Oct. 2014 (Rrox_v1/rhiRox1)           reciprocal best
Gorilla - Gorilla gorilla gorilla May 2011 (gorGor3.1/gorGor3)   reciprocal best
Green monkey - Chlorocebus sabaeus
                          Mar. 2014 (Chlorocebus_sabeus 1.1/chlSab2)  syntenic
Marmoset - Callithrix jacchus       Mar. 2009 (WUGSC 3.2/calJac3)     syntenic
Mouse lemur - Microcebus murinus    Jul. 2007 (Broad/micMur1) reciprocal best
Orangutan - Pongo pygmaeus abelii   July 2007 (WUGSC 2.0.2/ponAbe2)   syntenic
Proboscis monkey - Nasalis larvatus
                          Nov. 2014 (Charlie1.0/nasLar1)        reciprocal best
Rhesus - Macaca mulatta           Oct. 2010 (BGI CR_1.0/rheMac3)      syntenic
Squirrel monkey - Saimiri boliviensis Oct. 2011 (Broad/saiBol1) reciprocal best
Tarsier - Tarsius syrichta
                   Sep. 2013 (Tarsius_syrichta-2.0.1/tarSyr2) reciprocal best

---------------------------------------------------------------
To download a large file or multiple files from this directory, we recommend 
that you use rsync or ftp rather than downloading the files via our website:

   rsync -avz --progress \
        rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP17way/ ./

      OR

    ftp hgdownload.cse.ucsc.edu 
    user name: anonymous
    password: <your email address>
    go to the directory goldenPath/hg38/phyloP17way

To download multiple files from the UNIX command line, use the "mget" command. 
    mget <filename1> <filename2> ...
    - or -
    mget -a (to download all the files in the directory) 
Use the "prompt" command to toggle the interactive mode if you do not want 
to be prompted for each file that you download.

---------------------------------------------------------------
All the files in this directory are freely available for public use.
For data use restrictions regarding the genome assemblies used in this
annotation, see http://genome.ucsc.edu/goldenPath/credits.html.

---------------------------------------------------------------
References for phastCons and phyloP:

Pollard KS, Hubisz MJ, Siepel A. Detection of non-neutral substitution rates
on mammalian phylogenies. Genome Res. 2010 Jan;20(1):110-21.
(http://genome.cshlp.org/content/20/1/110.long)

Siepel A, Bejerano G, Pedersen JS, Hinrichs AS, Hou M, Rosenbloom K, Clawson
H, Spieth J, Hillier LW, Richards S, et al. Evolutionarily conserved elements
in vertebrate, insect, worm, and yeast genomes. Genome Res. 2005
Aug;15(8):1034-50.  (http://genome.cshlp.org/content/15/8/1034.full)

---------------------------------------------------------------
