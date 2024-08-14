## **Introduction**
Hey! This is the pipeline I used to assemble the genome of _Plestiodon fasciatus_.

Unlike my other genome assembly repository, this long read sequencing went very well. I extracted from blood stored in EDTA using a long read MagAttract kit. 


## **Genome Assembly with hifiasm -- Adapted from [Amanda Markee](https://github.com/amandamarkee/actias-luna-genome.git)**

[hifiasm](https://hifiasm.readthedocs.io/en/latest/) is a fast and easy haplotype-resolved de novo assembly software for PacBio HiFi reads
 - hifiasm documentation explaining input parameters: https://hifiasm.readthedocs.io/en/latest/pa-assembly.html
 - hifiasm documentation explaining output files: https://hifiasm.readthedocs.io/en/latest/interpreting-output.html

With purging:
```
#!/bin/sh
#SBATCH --job-name hoff_hifiasm
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --tasks-per-node=1 # Number of cores per node
#SBATCH --time=30:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org
#SBATCH --output=slurm-%j-%x.outt
#conda init
source ~/.bash_profile
conda activate fasciatus_ass

hifiasm -o hoff_hifi_assembly.asm -l 2 -t 32 /home/jhoffman1/mendel-nas1/fasciatus_genome/6354_import-dataset/hifi_reads/m84082_240409_161410_s4.hifi_reads.bc2049.fastq
```

## **Genome Assembly Quality Assessment with assemblystats.py -- Again, Adapted From [Amanda Markee](https://github.com/amandamarkee/actias-luna-genome.git)**

- After assembly with hifiasm, we can assess assembly quality using the [assemblystats.py script](https://github.com/MikeTrizna/assembly_stats/tree/0.1.4) created by Mike Trizna.
- The version of assemblystats.py used here was modified by Paul Frandsen (Brigham Young University).

First, I copied this script into my working directory, and called it assemblystats.py

```
#!/usr/bin/env python

import numpy as np
from itertools import groupby
import json
import sys


def fasta_iter(fasta_file):
    """Takes a FASTA file, and produces a generator of Header and Sequences.
    This is a memory-efficient way of analyzing a FASTA files -- without
    reading the entire file into memory.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    header: str
        The string contained in the header portion of the sequence record
        (everything after the '>')
    seq: str
        The sequence portion of the sequence record
    """

    fh = open(fasta_file)
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.upper().strip() for s in next(fa_iter))
        yield header, seq


def read_genome(fasta_file):
    """Takes a FASTA file, and produces 2 lists of sequence lengths. It also
    calculates the GC Content, since this is the only statistic that is not
    calculated based on sequence lengths.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    contig_lens: list
        A list of lengths of all contigs in the genome.
    scaffold_lens: list
        A list of lengths of all scaffolds in the genome.
    gc_cont: float
        The percentage of total basepairs in the genome that are either G or C.
    """

    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []
    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        if "NN" in seq:
            contig_list = seq.split("NN")
        else:
            contig_list = [seq]
        for contig in contig_list:
            if len(contig):
                gc += contig.count('G') + contig.count('C')
                total_len += len(contig)
                contig_lens.append(len(contig))
    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    seq_array = np.array(seq_lens)
    stats['sequence_count'] = seq_array.size
    stats['gc_content'] = gc_cont
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['longest'] = int(sorted_lens[0])
    stats['shortest'] = int(sorted_lens[-1])
    stats['median'] = np.median(sorted_lens)
    stats['mean'] = np.mean(sorted_lens)
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)
    for level in [10, 20, 30, 40, 50]:
        nx = int(stats['total_bps'] * (level / 100))
        csumn = min(csum[csum >= nx])
        l_level = int(np.where(csum == csumn)[0])
        n_level = int(sorted_lens[l_level])

        stats['L' + str(level)] = l_level
        stats['N' + str(level)] = n_level
    return stats


if __name__ == "__main__":
    infilename = sys.argv[1]
    contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    stat_output = {'Contig Stats': contig_stats,
                   'Scaffold Stats': scaffold_stats}
    print(json.dumps(stat_output, indent=2, sort_keys=True))
```

Next, I changed permissions as follows to allow execution permissions.
```
chmod +x assemblystats.py
```

Then, I produced a FASTA file from the initial GFA output files from the hifiasm assembly output for haplotype 1, haplotype 2, and total. I used the primary contig file, as indicated with asm.bp.p_ctg.fa (p_ctg meaning primary contig)
```
awk '/^S/{print ">"$2;print $3}' hoff_hifi_assembly.asm.bp.hap1.p_ctg.gfa > hoff_hifi_assembly.hap1-l0.ctg.fa
awk '/^S/{print ">"$2;print $3}' hoff_hifi_assembly.asm.bp.hap2.p_ctg.gfa > hoff_hifi_assembly.hap2-l0.ctg.fa
awk '/^S/{print ">"$2;print $3}' hoff_hifi_assembly.asm.bp.p_ctg.gfa > hoff_hifi_assembly.total-l0.ctg.fa
```

Lastly, I ran the assemblystats.py script on the newly generated fasta file of the fga in the format of scriptfilepath/scirptname.py nameofassembly.fa and save as a txt file 
```
./assemblystats.py ./hoff_hifi_assembly.hap1-l0.ctg.fa >> hoff_hifi_ass.stats-h1l0.txt
./assemblystats.py ./hoff_hifi_assembly.hap2-l0.ctg.fa >> hoff_hifi_ass.stats-h2l0.txt
./assemblystats.py ./hoff_hifi_assembly.total-l0.ctg.fa >> hoff_hifi_ass.stats-Tl0.txt
```

Let's check the stats on this assembly:

Contig Stats:
    L10: 0,
    L20: 1,
    L30: 1,
    L40: 2,
    L50: 3,
    N10: 244253184,
    N20: 227050938,
    N30: 227050938,
    N40: 212007450,
    N50: 166278287,
    gc_content: 45.59562809329659,
    longest: 244253184,
    mean: 45226658.02941176,
    median: 704375.0,
    sequence_count: 34,
    shortest: 8262,
    total_bps: 1537706373

Wow!! That looks really great. Since _Plestiodon fasciatus_ has 26 chromosomes and this assembly has 34 scaffolds, it is very close to chromosome level.
 

## **Checking Assembly Completeness with BUSCO**

[BUSCO](https://busco.ezlab.org/busco_userguide.html) is a program that estimates genome completeness based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs.

Since the worker nodes of the AMNH's computational clusters don't have access to the internet, it is necesarry to install BUSCO locally:
```
#clone repository
git clone https://gitlab.com/ezlab/busco.git
cd busco/

#pip install
python -m pip install .
```

Make sure you have all the [dependencies](https://busco.ezlab.org/busco_userguide.html#editing-busco-run-configuration) installed for the type of BUSCO run you are planning on running.

Then, I ran this shell file:
```
#!/bin/sh
#SBATCH --job-name Busco_Genomepfas
#SBATCH --nodes=20
#SBATCH --mem=100gb
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org

source ~/.bash_profile
conda activate fasciatus_ass
pfas="/home/jhoffman1/mendel-nas1/fasciatus_genome/6354_import-dataset/hifi_reads/hoff_hifi_assembly.total-l0.ctg.fa"
busco -m genome -i $pfas -o pfasBUSCO_saur -l sauropsida_odb10 -f --metaeuk --offline --download_path /home/jhoffman1/mendel-nas1/fasciatus_genome
```

BUSCO  results:

        ***** Results: *****

        C:95.1%[S:93.9%,D:1.2%],F:1.2%,M:3.7%,n:7480
        7110    Complete BUSCOs (C)
        7023    Complete and single-copy BUSCOs (S)
        87      Complete and duplicated BUSCOs (D)
        90      Fragmented BUSCOs (F)
        280     Missing BUSCOs (M)
        76      Total BUSCO groups searched

These BUSCO results are about as good as you can hope for!


