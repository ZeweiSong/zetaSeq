#! /bin/bash

python zetaseq/test/spawn_test.py
gzip tst_ctg_*.fa
mkdir contigs/
mv tst_ctg_1.fa.gz contigs/
mv tst_ctg_2.fq.gz contigs/
ls -hl contigs/

echo Create a working folder for DereCo
zetaseq/dereco create -i contigs/ -o tst/
ls -hl tst/

echo Calculate the NR set
zetaseq/dereco calculate -i tst/ -t 4
ls -hl tst/
ls -hl tst/dereco_final.fa
more tst/options.json

echo Add a new contig and update the working folder
mv tst_ctg_3.fa.gz contigs/
zetaseq/dereco update -i tst/

echo Calculate again
zetaseq/dereco calculate -i tst/ -t 4
ls -hl tst/
