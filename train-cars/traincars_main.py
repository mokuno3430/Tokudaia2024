import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import os.path

from traincars import common
from traincars import scaffold as scf
from traincars import highlight
from traincars import gff

argv = sys.argv

def func_set_axes( ax, size, seq ):
    ##set data range
    width = 0.05
    ax.set_xlim( -1 * width * seq.width, seq.width * ( 1 + width ) )
    ax.set_ylim( -1 * (size.bottom_margin + (len( seq.pos ) - 1 ) * seq.MARGIN ), size.top_margin )
    
    ##non-display axis
    ax.spines["right"].set_color("none")  
    ax.spines["left"].set_color("none")   
    ax.spines["top"].set_color("none")    
    ax.spines["bottom"].set_color("none") 
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.tick_params( length=0 )
    #ax.patch.set_alpha(0.0)

    
def func_print_messages():
    print( 'start' ) 
    print( ' '.join( argv ))

    
def main():

    args = common.get_args()
    func_print_messages( )
    seq = scf.input_scaffoldtsv( args.input, args.width, args.row_count )
    if( args.gff3 is not None and os.path.isfile( args.gff3 ) ):
        intergenic = gff.cal_intergenic( args.gff3, seq )
        seq.correct_width_by_genes( intergenic )
    elif( args.gff_excel is not None and os.path.isfile( args.gff_excel ) ):
        intergenic = gff.cal_intergenic_excel( args.gff_excel, seq )
        seq.correct_width_by_genes( intergenic )

    size=common.Size( seq, args.figure_size, args.gene_ratio )
    size.output_parameters()
    seq.output_parameters()
        
    fig = plt.figure( figsize=size.figsize_inch )
    ax = fig.add_subplot(111)
    fig.patch.set_alpha( 0.0 )
    func_set_axes( ax, size, seq )
    scf.plot_scaffolds( ax, seq, args.scaffold_font_size )

    #plot scale bar
    scalebar = scf.Scalebar( size, seq )
    if args.scale == 'legend' or args.scale == 'both' :
        scalebar.plot_legend( ax )
    if args.scale == 'tick' or args.scale == 'both' :
        scalebar.plot_tick( ax, seq, size )
    scalebar.output_parameters()

    ##plot highlight
    if args.highlight is not None:
        highlight.plot_highlight( seq, ax, size, args.highlight )

    ##plot gene
    if( args.gff3 is not None and os.path.isfile( args.gff3 ) ):
        gff.plot_genes( seq, ax, size, args.gff3, common.Color( args.gene_color, 1 ), args.gene_font_size )

    if( args.gff_excel is not None and os.path.isfile( args.gff_excel ) ):
        gff.plot_genes_excel( seq, ax, size, args.gff_excel, common.Color( args.gene_color, 1 ), args.gene_font_size )


    pdf_file = args.out + '.pdf'
    pp = PdfPages( pdf_file )
    pp.savefig( fig, bbox_inches='tight' )
    pp.close()
    

if __name__ == '__main__':
    main()
