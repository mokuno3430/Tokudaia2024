#!/bin/sh

gff=input.gff3
prefix=test
python ../traincars_main.py  -i input.tsv --out $prefix --gff3 $gff \
       --scaffold_font_size 12 --width 34000 --scale both --gene_font_size 4 --figure_size 8 4 1> log.$prefix
