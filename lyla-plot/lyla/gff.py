from lyla import common
import sys
import csv
import re
import matplotlib.pyplot as plt
import pandas as pd

class Gene:
    def __init__( self, scaffold, buf, height, tail_length, color, edge_color ):
        self.start, self.end, self.strand = self.__convert_position2coord( scaffold, buf )
        self.y_mid_coord = scaffold.convert_position2ycoord( scaffold.height /2 )
        self.width = self.end - self.start
        self.height = height
        self.color = common.Color( color, 1 )
        if edge_color == 'same' :
            self.edge_color = color
        else:
            self.edge_color = edge_color
        self.x, self.y  = self.__set_gene_coord( tail_length )

    def __set_gene_coord( self, tail_length ):
        x1 = self.start
        x3 = self.end
        y1 = self.y_mid_coord - self.height / 2
        #y1 = self.y_mid_coord - self.height * 0.8
        y2 = self.y_mid_coord
        y3 = self.y_mid_coord + self.height / 2
        #y3 = self.y_mid_coord + self.height * 0.8
        if( self.strand == '+' ): # => or >
            x2 = self.end - tail_length if self.width > tail_length else self.start
            x = [ x1, x2, x3, x2, x1 ]
            y = [ y1, y1, y2, y3, y3 ]
        elif( self.strand == '-' ): # <= or <
            x2 = self.start + tail_length if self.width > tail_length else self.end
            x = [ x1, x2, x3, x3, x2 ]
            y = [ y2, y1, y1, y3, y3 ]
        else:
            x = [ x1, x3, x3, x1 ]
            y = [ y1, y1, y3, y3 ]
        return x, y

    def __convert_position2coord( self, scaffold, buf ):
        if( scaffold.strand == '+' ):
            start = scaffold.convert_position2xcoord( int( buf[3] ) )
            end = scaffold.convert_position2xcoord( int( buf[4] ) )
            return start, end, buf[6]
        else:
            start = scaffold.convert_position2xcoord( int( buf[4] ) )
            end = scaffold.convert_position2xcoord( int( buf[3] ) )
            if( buf[6] == '+' ): #逆にする
                return  start, end, '-'
            elif( buf[6] == '-' ): #逆にする
                return  start, end, '+'
            else:
                return start, end, buf[6]

    def plot( self, ax ):
        if self.edge_color == 'none':
            ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0 )
        else:
            ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0.15, edgecolor=self.edge_color )
        """
        if self.color.color == '#ffffff' :
            ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0.15, edgecolor='black' )
        else:
            ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0.15 )
        """


def func_plot_gene_text( i, seqs, ax, buf, gene_font_size, height, rotation ):
    #遺伝子名の表示
    if gene_font_size == 0:
        return
    if re.search( r'Name=(\S+?);', buf[8] ):
        gene_name = re.search( r'Name=(\S+?);', buf[8] ).group(1)
        text_x = seqs[i][buf[0]].convert_position2xcoord( (int( buf[3] ) + int( buf[4]))/2 )
        text_y = seqs[i][buf[0]].convert_position2ycoord( -1 * height/2 )
        ax.text( text_x, text_y, gene_name, fontsize = gene_font_size, color = 'black', ha='center', va='top', rotation=rotation, style='italic' )


def plot_genes( seqs, ax, size, fn, thickness, gene_font_size, rotation, color, edge_color ):
    #color = '#4d4c61'
    RATIO = thickness
    color = color
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            if( buf[2] != 'gene' ):
                continue
            i = common.detect_index( buf[0], int( buf[3] ), int( buf[4] ), seqs )
            if( i == -1 ):
                continue
            height = seqs[i][buf[0]].height * RATIO
            tail_length = ( height / 2 ) * ( size.xlim_max / size.ylim_max ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
            if( len( buf ) == 10 ):
                color1 = common.Color( buf[9], 1 )
                aninstance = Gene( seqs[i][buf[0]], buf, height, tail_length, color1.color, edge_color )
            else:
                aninstance = Gene( seqs[i][buf[0]], buf, height, tail_length, color, edge_color )
            aninstance.plot( ax )

            func_plot_gene_text( i, seqs, ax, buf, gene_font_size, height, rotation )


def plot_genes_xlsx( seqs, ax, size, fn, thickness, gene_font_size, rotation, color, edge_color ):
    #color = '#4d4c61'
    RATIO = thickness
    color = color
    df = pd.read_excel( fn, header=None )
    for index, row in df.iterrows():
        buf = row.tolist()
        if( buf[2] != 'gene' ):
            continue

        i = common.detect_index( buf[0], int( buf[3] ), int( buf[4] ), seqs )
        if( i == -1 ):
            continue

        height = seqs[i][buf[0]].height * RATIO
        tail_length = ( height / 2 ) * ( size.xlim_max / size.ylim_max ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
        if( len( buf ) == 10 ):
            color1 = common.Color( buf[9], 1 )
            aninstance = Gene( seqs[i][buf[0]], buf, height, tail_length, color1.color, edge_color )
        else:
            aninstance = Gene( seqs[i][buf[0]], buf, height, tail_length, color, edge_color )
        aninstance.plot( ax )

        func_plot_gene_text( i, seqs, ax, buf, gene_font_size, height, rotation )
