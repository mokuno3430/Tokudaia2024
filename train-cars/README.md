## Directory Structure
- `traincars_main.py`: The main script to run the analysis.
- `traincars/`: Contains library modules for the main script.
- `sample_data/`: Includes sample data and a script for demonstration.

## Command line options
```
  -h, --help            show this help message and exit
  -i, --input str       scaffold_info.tsv
  --gff3 str            gff3 file
  --gff_excel str       gff3 file
  --highlight str       fileformat in the air
  --out str             optional: prefix of pdf file (default out)
  --margin_bw_scaffolds float
                        optional: (default -1, -1 means auto)
  --width float         optional: (default 0, 0 means auto)
  --row_count float     optional: (default 4)
  --gene_ratio float    optional: (default 8)
  --gene_color str      optional: (default grey)
  --figure_size [float ...]
                        optional: (default [6,4])
  --scaffold_font_size float
                        optional: (default 0, 0 means not shown)
  --gene_font_size float
                        optional: (default 0, 0 means not shown)
  --scale {legend,tick,both}
                        optional: [legend], [tick] or [both]
  -v, --version         show program's version number and exit
```