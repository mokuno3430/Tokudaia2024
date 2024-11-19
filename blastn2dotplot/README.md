
## Usage

```
blastn2dotplot.py [-h] -i1 str -i2 str --blastn str  
                  [--highlight str] [--highlight_crossed str]
                  [--out str] [--min_identity int] [--min_alignlen int]
                  [--line_width float] [--share] [--show_grid] 
                  [--xtitle_rotate int] [--ytitle_rotate int]
                  [--font_size float] [--figure_size [float ...]] 
                  [--hspace float] [--wspace float] [--h_alpha float]
                  [--tick_label_size float] [--tick_width int] 
                  [--colormap {0,1,2,3,4,5}] [-v]
```


## Options
```
  -h, --help            show this help message and exit
  -i1, --input1 str     sequence IDs of database at blastn (=row)
  -i2, --input2 str     sequence IDs of query at blastn (=column)
  --blastn str          blastn.tsv (-outfmt '6 std qlen slen' at blastn)
  --highlight str       optional: positions.tsv (scaffold start end color)
  --highlight_crossed str
                        optional: positions.tsv (scaffold start end color)
  --out str             optional: prefix of pdf file (default out)
  --min_identity int    optional: minimum sequence identity (default 90)
  --min_alignlen int    optional: minimum alignment length (default 100)
  --line_width float    optional: line width of dotplots (default 1.0)
  --share               optional: sharing axis scales among subplots.
  --show_grid           optional: show grid-line
  --xtitle_rotate int   optional: (default 0)
  --ytitle_rotate int   optional: (default 90)
  --font_size float     optional: (default 8)
  --figure_size [float ...]
                        optional: (default [8,8])
  --hspace float        optional: (default -1 means auto without --share, 0.04 with --share)
  --wspace float        optional: (default -1 means auto without --share, 0.04 with --share)
  --h_alpha float       optional: transparency ratio of highlights (default 0.3)
  --tick_label_size float
                        optional: (default 6)
  --tick_width int      optional: scale width of axis (bp) (default -1 means auto)
  --colormap {0,1,2,3,4,5}
                        optional: colormap for identity of alignments
                         ( 0:binary, 1:bone_r, 2:inferno_r, 3:hot_r, 4:YlGnBu, 5:original) (default 5)
  -v, --version         show program's version number and exit
```