#!/bin/bash

echo "Generating sequences for simulating layout module"

zetaseq/test/spawn_test_layout.py
ls -hl sim.fa
echo ""

echo "====== Dereco layout ======"
zetaseq/dereco layout -i sim.fa -o tst_layout/ -t 4 -l 1000 -d 0.99 -m 1000
echo ""

echo "====== FINISHED ======"
ls -hl tst_layout
