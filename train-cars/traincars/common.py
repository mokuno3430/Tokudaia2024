import argparse
from matplotlib import colors as mcolors
import matplotlib.colors as clr
import os.path

def get_args():
    parser = argparse.ArgumentParser( formatter_class=argparse.MetavarTypeHelpFormatter )
    
    parser.add_argument( '-i', '--input', help='scaffold_info.tsv', type=str, required=True )
    parser.add_argument( '--gff3', help='gff3 file', type=str )
    parser.add_argument( '--gff_excel', help='gff3 file', type=str )
    parser.add_argument( '--highlight', help='fileformat in the air', type=str )
    parser.add_argument( '--out', help='optional: prefix of pdf file (default out)', type=str, default='out' )
    parser.add_argument( '--margin_bw_scaffolds', help='optional: (default -1, -1 means auto)', type=float, default=-1 )
    parser.add_argument( '--width', help='optional: (default 0, 0 means auto)', type=float, default=0 )
    parser.add_argument( '--row_count', help='optional: (default 4)', type=float, default=4 )
    parser.add_argument( '--gene_ratio', help='optional: (default 8)', type=float, default=8 )
    parser.add_argument( '--gene_color', help='optional: (default grey)', type=str, default='grey' )
    parser.add_argument( '--figure_size', help='optional: (default [6,4])', type=float, nargs='*', default=[6,4] )
    parser.add_argument( '--scaffold_font_size', help='optional: (default 0, 0 means not shown)', type=float, default=0 )
    parser.add_argument( '--gene_font_size', help='optional: (default 0, 0 means not shown)', type=float, default=0 )
    parser.add_argument('--scale', help='optional: [legend], [tick] or [both]', type=str, choices=['legend', 'tick', 'both'], default='legend')
    parser.add_argument( '-v', '--version', action='version', version='%(prog)s 1.0', default=False)
    args = parser.parse_args()

    return ( args )


class Size:
    def __init__( self, seq, size, ratio ):
        self.bottom_margin = 1.25
        self.top_margin = 1.0
        self.figsize_inch = [ size[0], size[1] ]
        self.gene_ratio = ratio
    def output_parameters( self ):
        print( '\n##Size parameters:' )
        print( '  figure size inch: %.1f,%.1f' % ( self.figsize_inch[0], self.figsize_inch[1] ))
        print( '  top_margin: %.2f' % ( self.top_margin ))
        print( '  bottom_margin: %.2f' % ( self.bottom_margin ))
        print( '  gene ratio: %d' % ( self.gene_ratio ))
        print( '' )


class Text:
    def __init__( self, label, size, x, y, ha, va, color='black' ):
        self.label = label
        self.size = size
        self.origin_x = x
        self.origin_y = y
        self.color = color
        self.ha = ha
        self.va = va

    def output( self, ax ):
        if self.label == 'BLANK':
            return
        if self.size == 0:
            return
        ax.text( self.origin_x, self.origin_y, self.label, fontsize = self.size, color = self.color, ha=self.ha, va=self.va )


class Polygon:
    def __init__( self, x, y, color, alpha ):
        self.x = x
        self.y = y
        self.color = Color( color, alpha )

    def plot( self, ax ):
        ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0 )


class Color:
    replaced_color = 'grey'
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

    def __init__( self, color, alpha ):
        self.color = self.__evaluate_input( color )
        self.alpha = alpha

    def __evaluate_input( self, color ):
        if clr.is_color_like( color ):
            return color
        else:
            print( "invalid color :{c}".format( c = color ) )
            return Color.replaced_color

"""
class Color:
    replaced_color = 'grey'
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

    def __init__( self, color, alpha ):
        self.color = self.__evaluate_input( color )
        self.alpha = alpha

    def __evaluate_input( self, color ):
        if color in Color.colors:
            return color
        if color[0:1] == '#' and len( color ) == 7:
            RGB = (int(color[1:3],16),int(color[3:5],16),int(color[5:7],16))
            if( RGB[0] <= 255 and RGB[1] <= 255 and RGB[2] <= 255 ):
                return color
        print( "invalid color :{c}".format( c = color ) )
        return Color.replaced_color
"""

def cal_total_length( seqs ):
    length = []
    for sample in seqs:
        total = 0
        for scaf in sample:
            total += sample[scaf].length
        length.append( total )
    return length


def detect_index( ID, start, end, seq ):
    if( seq.name == 'BLANK' ):
        return -1
    if( ID == seq.ID and seq.start <= start and end <= seq.end ):
        return 0
    else:
        return -1
