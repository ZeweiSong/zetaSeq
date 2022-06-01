#! /bin/bash
clear

echo ===== STEP 1 =====
echo Generate random sequences and three contigs
rm -r tst/
rm -r contigs/
python zetaseq/test/spawn_test.py
gzip tst_ctg_*.fa
mkdir contigs/
mv tst_ctg_1.fa.gz contigs/
mv tst_ctg_2.fa.gz contigs/
ls -hl contigs/

read -n 1 -s -r -p "Press any key to continue"
echo ""
echo ===== STEP 2 =====
echo Create a working folder for DereCo
zetaseq/dereco create -i contigs/ -o tst/ -s
ls -hl tst/

read -n 1 -s -r -p "Press any key to continue"
echo ""
echo ===== STEP 3 =====
echo Calculate the none-redundant set
zetaseq/dereco devour -i tst/ -t 4
ls -hl tst/
ls -hl tst/dereco_final.fa

read -n 1 -s -r -p "Press any key to continue"
echo ""
echo ===== STEP 4 =====
echo Add a new contig and update the working folder
mv tst_ctg_3.fa.gz contigs/
zetaseq/dereco update -i tst/

read -n 1 -s -r -p "Press any key to continue"
echo ""
echo ===== STEP 5 =====
echo Calculate the NR set again
zetaseq/dereco devour -i tst/ -t 4
ls -hl tst/

read -n 1 -s -r -p "Press any key to continue"
echo ""
echo ===== STEP 6 =====
echo Check the status of a working folder
zetaseq/dereco status -i tst/
