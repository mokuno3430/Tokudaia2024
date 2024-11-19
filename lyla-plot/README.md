## Directory Structure
- `lyla/`: Contains library modules for the main script.
- `lyla_main_for_local.py`: The main script to run the analysis.

## Command line options
```
-h, --help            show this help message and exit
  -i, --input [str ...]
                        scaffold_info.tsv (arranged from bottom to top)
  --xlsx str            scaffold_info.xlsx
  -a, --alignment [str ...]
                        alignment.tsv
  --blastn [str ...]    blastn.tsv (-outfmt 6)
  --lastz [str ...]     lastz.tsv (--format=general)
  --mummer [str ...]    show-coords.tsv (--format=show-coords -H)
  --min_identity int    optional: (default -1, -1 means auto)
  --max_identity int    optional: (default 100)
  --min_alignment_len int
                        optional: (default 0)
  --alignment_height float
                        optional: (default 1.5)
  --alignment_alpha float
                        optional: (default 0.5)
  --colormap {0,1,2,3,4,5}
                        optional: colormap for identity of alignments ( 0:binary, 1:bone_r, 2:inferno_r, 3:hot_r, 4:YlGnBu, 5:original) (default 0)
  --gff3 [str ...]      gff3 files
  --gff_xlsx [str ...]  excel files
  --gene_thickness float
                        optional: (default 3)
  --gene_font_size float
                        optional: (default 0)
  --gene_font_rotation float
                        optional: (default 0)
  --gene_color str      optional: (default black)
  --gene_edge_color str
                        optional: same means same as fill color (default black)
  --hist [str ...]      histogram file (scaffold_ID position value)
  --mark_v str          fileformat in the air
  --scaffold_layout {left,center,right}
                        optional: (default center)
  --scaffold_color str  optional: (default grey)
  --margin_bw_scaffolds float
                        optional: (default -1, -1 means auto)
  --xlim_max float      optional: (default -1, -1 means auto)
  --scaffold_font_size float
                        optional: (default 0, 0 means not shown)
  --figure_size [float ...]
                        optional: (default [8,6])
  --seq_thickness float
                        optional: (default 0.1)
  --out str             optional: prefix of pdf file (default out)
  ```