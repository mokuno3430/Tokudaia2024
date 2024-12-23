import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clr
from matplotlib.backends.backend_pdf import PdfPages
import csv
import math
from lyla import common
import pandas as pd

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

class Scaffold:
    def __init__( self, order, ID, start, end, strand, name, color, thickness ):
        self.order = order
        self.ID = ID
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        self.length = self.end - self.start
        self.origin_x = 0
        self.origin_y = 0
        self.height = thickness
        if not clr.is_color_like( color ):
            print(f"Error: {color} is invalid color. Use default color for {self.ID}.")
            color = 'grey'
        else:
            color = color
        #color = '#95949a'
        alpha = 0 if self.name == 'BLANK' else 1
        self.color = common.Color( color, alpha )

    def set_origins( self, x, y ):
        self.origin_x = x
        self.origin_y = y

    def convert_position2xcoord( self, pos ):
        if( self.strand == '+' ):
            x = self.origin_x + pos - self.start
        else:
            x = self.origin_x + self.length - pos + self.start
        return x

    def convert_position2ycoord( self, pos ):
        return self.origin_y + pos

    def plot_line( self, ax ):
        x = [self.origin_x, self.origin_x + self.length, self.origin_x + self.length, self.origin_x ]
        y = [self.origin_y, self.origin_y, self.origin_y + self.height , self.origin_y + self.height ]
        ax.fill( x, y, color=self.color.color, lw=0, alpha=self.color.alpha )

    def plot_line_on_histogram( self, ax, h_height ):
        x = [self.origin_x, self.origin_x + self.length, self.origin_x + self.length, self.origin_x ]
        y = [self.origin_y + h_height, self.origin_y + h_height, self.origin_y + self.height + h_height, self.origin_y + self.height + h_height ]
        ax.fill( x, y, color=self.color.color, lw=0, alpha=self.color.alpha )


def input_scaffold_tsv( scf_files, color, thickness ):
    seqs = []
    for fn in scf_files:
        with open( fn , 'r' ) as file:
            buf = csv.reader( file, delimiter = '\t' )
            n = 0
            seqlist = {}
            for scaffold, scaf_start, scaf_end, scaf_strand, scaf_name in buf:
                aninstance = Scaffold( n, scaffold, int( scaf_start ), int( scaf_end ), scaf_strand, scaf_name, color, thickness )
                seqlist[ scaffold ] = aninstance
                n += 1
            seqs.append( seqlist )
    return seqs


def input_scaffold_xlsx( xlsx, color, thickness ):
    seqs = []
    df = pd.read_excel( xlsx )
    max_value = df['n'].max()
    seqs = [None] * ( max_value + 1 )

    for row in df.itertuples(index=False):
        seqlist = {}
        aninstance = Scaffold( 0, row.ID, int( row.start ), int( row.end ), row.strand, row.name, color, thickness )
        seqlist[ row.ID ] = aninstance
        seqs[ row.n ] = seqlist

    """
    for fn in scf_files:
        with open( fn , 'r' ) as file:
            buf = csv.reader( file, delimiter = '\t' )
            n = 0
            seqlist = {}
            for scaffold, scaf_start, scaf_end, scaf_strand, scaf_name in buf:
                aninstance = Scaffold( n, scaffold, int( scaf_start ), int( scaf_end ), scaf_strand, scaf_name )
                seqlist[ scaffold ] = aninstance
                n += 1
            seqs.append( seqlist )
    """
    return seqs


def plot_scaffolds( ax, seqs, s_font_size ):
    for i in range( len( seqs )):
        for scaf in seqs[i]:
            ##plot scaffold
            seqs[i][scaf].plot_line( ax )
            y = seqs[i][scaf].origin_y + seqs[i][scaf].height /2
            aninstance = common.Text( seqs[i][scaf].name, s_font_size, 0, y, 'right', 'center' )

            aninstance.output( ax )


class Scalebar:
    def __init__( self, size ):
        self.width = self.set_width( size )
        self.bar_color = 'black'
        self.legend = self.set_label()
        self.origin_x, self.origin_y = self.set_origin( size )
        self.label = common.Text( self.legend, 8, self.origin_x + self.width/2, self.origin_y - 0.05, 'center', 'top' )

    def set_label( self ):
        if 6 <= math.log10( self.width ):
            return str( int( self.width/ math.pow( 10, 6 )) ) + ' Mbp'
        if 3 <= math.log10( self.width ):
            return str( int( self.width/ math.pow( 10, 3 )) ) + ' kbp'
        return str( int( self.width )) + ' bp'

    def set_origin( self, size ):
        x = size.xlim_max - self.width
        y = size.bottom_margin / 2
        return x, y

    def set_width( self, size ):
        if  math.pow( 10, int( math.log10( size.xlim_max )))/ size.xlim_max <= 0.2 :
            return math.pow( 10, int( math.log10( size.xlim_max )))
        if math.pow( 10, int( math.log10( size.xlim_max )))/ size.xlim_max <= 0.4 :
            return math.pow( 10, int( math.log10( size.xlim_max )))/2
        return math.pow( 10, int( math.log10( size.xlim_max )))/5

    def plot( self, ax ):
        ax.plot( [ self.origin_x, self.origin_x + self.width ], [ self.origin_y, self.origin_y ], color=self.bar_color, lw=1 )
        self.label.output( ax )

    def output_parameters( self ):
        print( '##Scalebar paramenters:' )
        print( '  width: %d' % ( self.width ))
        print( '  origin_x %.2f' % ( self.origin_x ))
        print( '  origin_y: %.2f' % ( self.origin_y ))
        print( '' )

