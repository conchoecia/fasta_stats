# fasta_stats

This is the `fasta_stats` program that is included with the Meraculous genome assembly software. I am not sure of who the original author of this software is.

## Installation

Install using a typical `git clone`, then `make` routine.

```
git clone https://github.com/conchoecia/fasta_stats.git
cd fasta_stats
make
```

Then, open `~/.bash_profile` or `~/.bashrc`, or however you add new programs to your path, 
then add the program to your path like so. Change the absolute path to match wherever your files are saved.

```
PATH=$PATH:/path/to/fasta_stats/bin
export PATH
```

Restart your terminal (or close your SSH, then login again), and you should be able to run the program.

## Usage

Simply run the following command to get stats about a fasta file, typically a genome assembly.

```
fasta_stats genome_assembly_of_interest.fa
```

The output will look like this. Keep in mind that N/L50 are accidentally swapped. This is a bug in the software that I have not yet fixed.

```
Main genome scaffold total: 45
Main genome contig total:   352
Main genome scaffold sequence total: 110.7 MB
Main genome contig sequence total:   110.7 MB (->  0.0% gap)
Main genome scaffold N/L50: 6/8.5 MB
Main genome contig N/L50:   61/580.4 KB
Number of scaffolds > 50 KB: 13
% main genome in scaffolds > 50 KB: 99.5%

 Minimum    Number    Number     Total        Total     Scaffold
Scaffold      of        of      Scaffold      Contig     Contig
 Length   Scaffolds  Contigs     Length       Length    Coverage
--------  ---------  -------  -----------  -----------  --------
    All        45        352  110,691,255  110,660,632    99.97%
   1 kb        45        352  110,691,255  110,660,632    99.97%
 2.5 kb        45        352  110,691,255  110,660,632    99.97%
   5 kb        44        351  110,686,898  110,656,275    99.97%
  10 kb        36        343  110,629,701  110,599,078    99.97%
  25 kb        23        330  110,434,894  110,404,271    99.97%
  50 kb        13        319  110,109,091  110,078,568    99.97%
 100 kb        13        319  110,109,091  110,078,568    99.97%
 250 kb        13        319  110,109,091  110,078,568    99.97%
 500 kb        13        319  110,109,091  110,078,568    99.97%
   1 mb        13        319  110,109,091  110,078,568    99.97%
 2.5 mb        13        319  110,109,091  110,078,568    99.97%
   5 mb        13        319  110,109,091  110,078,568    99.97%
```
