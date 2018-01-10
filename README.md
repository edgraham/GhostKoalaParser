---
layout: post
title: "Importing GhostKoala annotations into anvi'o"
excerpt: "Kegg "
modified: 2018-01-10
tags: [tutorial, GhostKoala]
categories: [anvio]
comments: true
authors: [elaina]
---

{% capture images %}{{site.url}}/images/2018-01-10-working-with-GhostKoala{% endcapture %}

{% include _toc.html %}

This tutorial walks you through annotating your contigs with [GhostKOALA](http://www.kegg.jp/ghostkoala/), and compiling the results into a format easily importable to your anvi'o contigs database using `Kegg-to-Anvio.py`. You can download all of the assocaited files for this workflow [here](https://github.com/edgraham/GhostKoalaParser).

All you'll need to run this tutorial is an installation of anvi'o, python, and internet access. So first lets download all of files we'll need for this tutorial. They are all on github so we can clone the repository as shown below:

In python you'll need the pandas and BioPython modules. I was running the following versions:

[Pandas](https://pandas.pydata.org/): v0.22.0
[BioPython](http://biopython.org/): v1.70

``` bash
 $ git clone https://github.com/edgraham/GhostKoalaParser.git
```

This will produce a directory containing all of the files and scripts we'll use throughout this tutorial.

This tutorial assumes that you have assembled your reads into contigs, lets call your contigs file `assembly.fasta`, and have followed the [metagenomic workflow tutorial](http://merenlab.org/2016/06/22/anvio-tutorial-v2/) and generated an anvi'o contigs database, lets call it `CONTIGS.db`.

So lets start from the beginning. 

### Export Anvio Gene Calls

The first step is exporting amino acid sequences from your anvi'o contigs database. 

``` bash
 $ anvi-get-aa-sequences-for-gene-calls -c CONTIGS.db -o protein-sequences.fa
```


### Run GhostKOALA

{:.notice}
GhostKOALA does have an upload limit of 300MB. If your anvi'o database is rather large and you find that your `protein-sequences.fa` file is larger than that you can split the file and append the GhostKOALA outputs after this step.

Before we run GhostKOALA we need to make a few modifications to our 'protein-sequences.fa' file. The parser that GhostKOALA uses for fasta sequences appears to have an unfriendly relationship with sequence id's that begin with a number. An example is shown below:

```

>0
MAEYQNIFTQVQVQGPAEMGVDPAGTLSRERTNGTSFSKLAGLFGNAQLGPIYLGTFGLI
SLVTGFAWFFMVGLSFWDQVDYSPALFLRELFWLALEPPAEEYGLSIPPMAEGGYFLLAS
FFLLISVICWWVRTYLRAEELGMGKHVAWAFASAIWLFLVLGLFRPILMGSWSEMVPYGI
FPHLDWTNLFSLTYGNLFYNPFHALSIVFLYGSALLFAMHGATILAVSRYGGEREIEQIV
DRGTASERAALFWRWTM
>1
TYLRAEELGMGKHVAWAFASAIWLFLVLGLFRPILMGSWSEMVPYGI
FPHLDWTNLFSLTYGNLFYNPFHALSIVFLYGSALLFAMHGATILAVSRYGGEREIEQIV
DRGTASERAALFWRWTMGFNATMEGIHRWAWWFAVLTTLTGGIGILLTGTVVDNWFIWAQ
DHGYAPLN

```

To remedy this use the `fix_seq.py` script in the repository you clones from github earlier.

```

fix_seqs.py protein-sequences.fa > protein-sequences.renammed.fa

```

You will now have a new file `protein-sequences.renammed.fa` that contains all your gene calls with the id's ammended. This file is now ready to be submitted to GhostKOALA for annotation. To do this go to the [GhostKOALA webserver](http://www.kegg.jp/ghostkoala/). When you get to the page there will be a section that says **Upload query amino acid sequences in FASTA format**. You'll select `Choose File` and upload the file `protein-sequences.renammed.fa`. You'll be asked to input an email address (you can only run one instance of GhostKOALA at a time per email address!). Before the run starts you'll get an email that will look like this:

 ```
 Your GhostKOALA job request
  Query dataset: 1 entries
  KEGG database to be searched: c_family_euk+genus_prok


Please click on the link below to either submit or cancel your job.

https://urldefense.proofpoint.com/v2/url?u=http-3A__www.kegg.jp_kegg-2Dbin_blastkoala-5Fsubmit-3Fid-3D9ae0dd8ccf3d6e7559c1e4fc00d3d47f7fa36e63-26passwd-3DVApgKX-26type-3Dghostkoala&d=DwICAg&c=clK7kQUTWtAVEOVIgvi0NU5BOUHhpN0H8p7CSfnc_gI&r=ODewbdal_BGCBKDYMtgzwQ&m=JshUeXOCVSmUSAysX_ToLHXiWYwfGwb0fAvdsJnnenk&s=WxeiNpY2tIYugW4SI9DY5i5djiycNYoiqj0XV6ZsSaI&e= (Submit)

https://urldefense.proofpoint.com/v2/url?u=http-3A__www.kegg.jp_kegg-2Dbin_blastkoala-5Fcancel-3Fid-3D9ae0dd8ccf3d6e7559c1e4fc00d3d47f7fa36e63-26passwd-3DVApgKX-26mode-3Dcancel-26type-3Dghostkoala&d=DwICAg&c=clK7kQUTWtAVEOVIgvi0NU5BOUHhpN0H8p7CSfnc_gI&r=ODewbdal_BGCBKDYMtgzwQ&m=JshUeXOCVSmUSAysX_ToLHXiWYwfGwb0fAvdsJnnenk&s=KV1zTnS11QX2jG4iHnCYXmXdj599W2KCPciTbckyPOc&e= (Cancel)

If no action is taken within 24 hours, your request will be deleted.
```

 click the link saying `submit` and your run will begin processing.
 

{:.notice}
The KEGG parser I wrote also has the option to combine your GhostKOALA results with interproscan. If you want to incorporate both genecalls from interproscan and KEGG in the same table follow the tutorial [here](http://merenlab.org/2016/06/18/importing-functions/), **but make the following modification** and run interproscan with the flags `-f tsv --goterms --iprlookup --pathways`. 

### Generate the KEGG orthology table

**Note this step is optional and is there for anyone who is curious about how I generated this file :)**

In the repository you cloned earlier there is a file called `KO_Orthology_ko00001.txt`. If you take a peak at this file it looks like:

```
Metabolism      Overview        01200 Carbon metabolism [PATH:ko01200]  K00844  HK; hexokinase [EC:2.7.1.1]
Metabolism      Overview        01200 Carbon metabolism [PATH:ko01200]  K12407  GCK; glucokinase [EC:2.7.1.2]
Metabolism      Overview        01200 Carbon metabolism [PATH:ko01200]  K00845  glk; glucokinase [EC:2.7.1.2]
Metabolism      Overview        01200 Carbon metabolism [PATH:ko01200]  K00886  ppgK; polyphosphate glucokinase [EC:2.7.1.63]
Metabolism      Overview        01200 Carbon metabolism [PATH:ko01200]  K08074  ADPGK; ADP-dependent glucokinase [EC:2.7.1.147]
```
This is the information that we'll use to convert the KEGG Orthology assignments to function. Because the KEGG database is currently working under a subscription model I had to find a workaround to access the information to match the Orthologies with function. To do this you have to go to the [KEGG website](http://www.genome.jp/kegg-bin/get_htext?ko00000.keg) and download the htext file. It will download as something along the lines of `ko00001.keg`. Next you can use this very ugly code snippet to parse that file into a readble format.

```
kegfile="ko00001.keg"
while read -r prefix content; 
do case "$prefix" in A) col1="$content" ;;B) col2="$content" ;; C) col3="$content";; D) echo -e "$col1\t$col2\t$col3\t$content";; esac; 
done < <(sed '/^[#!+]/d;s/<[^>]*>//g;s/^./& /' < "$kegfile") > KO_Orthology_ko00001.txt

```
What this is doing is going through the hierarchical ‘.keg’ file you downloaded and extracting the different layers. The output should be a tab delimited file where column 1 corresponds to the broadest classification and column 5 corresponds to the gene itself. You may notice that some of the keg identifiers appear multiple times in this parsed folder. This is because many of the metabolism genes are constituents of multiple pathways.

Now that you know how I generated the `KO_Orthology_ko00001.txt` lets get back to the fun part and import our functions!

### Parsing the results from GhostKOALA

Once GhostKOALA has finished running it will send you an email with a link to the results. Download the annotation file, wthe file will be called `user_ko.txt`. The results should look like this:

``` bash
 $ head
 genecall_0       K01923
 genecall_1
 genecall_2       K03611
 genecall_3
 genecall_4       K01952
 genecall_5       K01693
 genecall_6       K00817
 genecall_7       K00013
 genecall_8       K00765

```

For further analysis you'll need to make sure that row 1 contains both a contig and accession id. if for example the file looks like this:

``` bash
 $ head
 genecall_0       
 genecall_1
 genecall_2       K03611
 genecall_3
 genecall_4       K01952
 genecall_5       K01693
 genecall_6       K00817
 genecall_7       K00013
 genecall_8       K00765

```
You will get an error in python when you run `Kegg-to-Anvio.py`. The error will say something along the lines of:

`pandas.parser.CParserError: Error tokenizing data. C error: Expected 1 fields in line 2, saw 2`

This error is because when Pandas reads the tab delimited file into a dataframe it sets a precedent using the first row, so if the first row has no accession data then pandas will assume all subsequent rows are only intended to have a single column, rather than two. 

Once you have checked and fixed this issue we can convert this table into a format that is easily importable to your anvi'o database. We'll use the `Kegg-to-Anvio.py` parser. The help menu is shown below:

```bash
 $ Kegg-to-Anvio.py -h
 
 usage: Kegg-to-Anvio.py [-h] [--KeggDB KEGGDB] [-i I]
                        [--interproscan INTERPROSCAN] [-o O]

Combines annotation Data for input to anvio

optional arguments:
  -h, --help            show this help message and exit
  --KeggDB KEGGDB       identify the Kegg Orthology file (modified from htext
                        using given bash script)
  -i I                  specify the file containing GhostKoala Results
  --interproscan INTERPROSCAN
                        interproscan results
  -o O                  Specify an output file
```

Now the next step it pretty simple. If you don't have interproscan results you want to include you'll run this as:

```bash
$ Kegg-to-Anvio.py --KeggDB KO_Orthology_ko00001.txt -i user_ko.txt -o KeggAnnotations-AnviImportable.txt
```

if you have interproscan results, lets say in a output file called `interproscan-results.txt` (ran with the flags `-f tsv --goterms --iprlookup --pathways`) then you will run this as:

```bash
$ Kegg-to-Anvio.py --KeggDB KO_Orthology_ko00001.txt -i user_ko.txt -o KeggAnnotations-AnviImportable.txt --interproscan interproscan-results.txt
```

Now you will have a file `KeggAnnotations-AnviImportable.txt` that can be imported into anvio using `anvi-import-functions`!! 

```bash
$ anvi-import-functions -c CONTIGS.db -i KeggAnnotations-AnviImportable.txt

```
{:.notice}
If you decide you want to knock two birds out with one stone you can also take the taxonomy data produced in ghost koala and convert that to an anvio importable format using the script `GhostKoalaTaxonomy-to-Anvio.py` found in the github repository we downloaded. The taxonomy file will download as a file called `user.out.top`. You can then run the parser like this: `GhostKoalaTaxonomy-to-Anvio.py user.out.top KeggTaxonomy.txt`. Then import into your anvi'o database `anvi-import-taxonomy -c CONTIGS.db -i KeggTaxonomy.txt`.


And we now have KEGG functions in our Anvi'o database! Hopefully (in the near future) we can find a better way to do this without us having to do all the convuluted steps with GhostKOALA. In the meantime 

<div style="margin:50px">&nbsp;</div>
