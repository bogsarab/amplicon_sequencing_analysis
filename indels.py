from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os
import re
from configs import RunConfig, FileConfig
import sys

def indel_analysis(runconfig: RunConfig, fileconfig: FileConfig):
    samples = runconfig.barcodes.dict().values()
    barcodes = runconfig.barcodes.dict()
    ref_seq = runconfig.ref_seq    
    odir = fileconfig.out_dir
    idir = fileconfig.out_dir
    fdir = fileconfig.figures_dir

    if not os.path.isdir(idir):
        sys.exit("No input directory found.")

    for bc in barcodes.keys():

        print(f"Processing {barcodes[bc]} indels...")
        idx1, idx2 = 0,0
        indels = []

        file = os.path.join(idir, (barcodes[bc] + "_final.fastq"))

        fregex = bc
        rregex = str(Seq.reverse_complement(Seq(bc)))

        for rec in SeqIO.parse(file, "fastq"):
            idx1 += 1

            f = re.search(fregex, str(rec.seq[0:9]))
            r = re.search(rregex, str(rec.seq[-10:-1]))
            if f and r:

                start = f.end()
                end = -(r.end() + 2)


                f = re.search(fregex, str(rec.seq[0:9]))
                r = re.search(rregex, str(rec.seq[-10:-1]))

                start = f.end()
                end = -(r.end()+2)

                rec = rec[start:end]

                if len(rec.seq[25:-17]) is not 82:
                    idx2 += 1
                    indels.append(rec[25:-17])

                    alignment = pairwise2.align.globalms(rec.seq[25:-17], Seq(ref_seq), 2, -1, -1, -.1)

                    print(rec)
                    print(format_alignment(*alignment[0]))
                    print("\n")
        print(f"In sample {barcodes[bc]}, ratio of indels ia {idx2/idx1}.")
        print(f"Saving indel records...")
        SeqIO.write(indels, os.path.join(odir, (barcodes[bc] + "_indels.fastq")), "fastq")
