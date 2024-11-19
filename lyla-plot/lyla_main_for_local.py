import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os.path

from lyla import alignment, common, gff, histogram, mark_v, scaffold
from lyla import scaffold as scf


argv = sys.argv

def func_set_axes( ax, size ):
    ##set data range
    ax.set_xlim(0, size.xlim_max)
    ax.set_ylim(0, size.ylim_max)

    ##non-display axis
    ax.spines["right"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.tick_params( length=0 )

    ax.patch.set_alpha(0.0)


def func_print_messages():
    print( 'start' ) 
    print( ' '.join( argv ))


def main():

    args = common.get_args()
    func_print_messages( )

    if args.input:
        seqs = scf.input_scaffold_tsv( args.input, args.scaffold_color, args.seq_thickness )
    elif args.xlsx:
        seqs = scf.input_scaffold_xlsx( args.xlsx, args.scaffold_color, args.seq_thickness )


    size=common.Size( seqs, args.margin_bw_scaffolds, args.xlim_max, args.alignment_height, args.figure_size, args.seq_thickness, args.gene_thickness )
    #histograms = histogram.set_space( args.hist, seqs, size.histogram_height )
    size.set_histogram_space( seqs, args.hist )
    size.set_scaffold_layout( seqs, args.scaffold_layout )
    size.output_parameters()

    fig = plt.figure( figsize=size.figsize_inch )
    ax = fig.add_subplot(111)
    # fig.patch.set_alpha( 0.0 )
    func_set_axes( ax, size )

    scf.plot_scaffolds( ax, seqs, args.scaffold_font_size )

    ##plot scale bar
    scalebar = scf.Scalebar( size )
    scalebar.plot( ax )
    scalebar.output_parameters()

    ##plot alignment
    max_identity = args.max_identity
    input_formats = [ args.alignment, args.blastn, args.lastz, args.mummer ]
    func_plot_alignmment = [ alignment.plot_alignment4original, alignment.plot_alignment4blastn, alignment.plot_alignment4lastz, alignment.plot_alignment4mummer ]
    valid_files = alignment.count_alignment_files( args )
    if valid_files == 0:
        pass
    else:
        min_identity = alignment.set_min_identity( args )
        ##set colormap
        heatmap = alignment.Colormap( min_identity, max_identity, args.colormap, args.alignment_alpha )
        heatmap.output_parameters()
        ##set and plot colormap legend
        heatmap_legend = alignment.Colorbox( size )
        heatmap_legend.plot( ax, heatmap )
        heatmap_legend.output_parameters()

        for files, func_plot in zip( input_formats, func_plot_alignmment ):
            if files is None:
                continue
            for fn in files:
                if not os.path.isfile( fn ):
                    continue
                func_plot( seqs, ax, heatmap, size, fn, args )

    ##plot mark_v
    if args.mark_v is not None:
        if os.path.isfile( args.mark_v ):
            mark_v.plot_mark_v( seqs, ax, size, args.mark_v )

    ##plot gene
    if args.gff3 is not None:
        for fn in args.gff3:
            if not os.path.isfile( fn ):
                continue
            gff.plot_genes( seqs, ax, size, fn, args.gene_thickness, args.gene_font_size, args.gene_font_rotation, args.gene_color, args.gene_edge_color )

    if args.gff_xlsx is not None:
        for fn in args.gff_xlsx:
            if not os.path.isfile( fn ):
                continue
            gff.plot_genes_xlsx( seqs, ax, size, fn, args.gene_thickness, args.gene_font_size, args.gene_font_rotation, args.gene_color, args.gene_edge_color )



    ##plot histogram
    histogram.plot_background( seqs, ax, size )
    if args.hist is not None:
        for fn in args.hist:
            if not os.path.isfile( fn ):
                continue
            histogram.plot_histogram( seqs, ax, size, fn )

    pdf_file = args.out + '.pdf'
    pp = PdfPages( pdf_file )
    pp.savefig( fig, bbox_inches='tight' )
    pp.close()


if __name__ == '__main__':
    main()
