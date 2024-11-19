#!/bin/sh

db=db.txt
query=query.txt
blastn=blastn_samples_outfmt6.txt
highlight=highlight.txt

prefix=test1
python ../blastn2dotplot.py -i1 $db -i2 $query \
       --blastn $blastn --min_identity 95 --min_alignlen 1000 --show_grid \
       --ytitle_rotate 90 --out $prefix --share --wspace 0.08 --hspace 0.08 --font_size 8  1> log.${prefix}

prefix=test2
python ../blastn2dotplot.py -i1 $db -i2 $query \
       --blastn $blastn --min_identity 95 --min_alignlen 1000 \
       --highlight $highlight \
       --ytitle_rotate 90 --out $prefix --wspace 0.16 --hspace 0.16 --font_size 8  1> log.${prefix}
