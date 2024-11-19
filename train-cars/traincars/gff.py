import matplotlib.pyplot as plt
import sys
import csv
from traincars import common
import re
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill

class Gene:

    def __init__( self, scaffold, buf, height, tail_length, color ):
        self.start, self.end, self.strand = self.__convert_position2coord( scaffold, buf )
        self.y_mid_coord = scaffold.convert_position2ycoord( int( buf[3] ), scaffold.height /2 )
        self.y_mid_coord_end = scaffold.convert_position2ycoord( int( buf[4] ), scaffold.height /2 )
        self.width = self.end - self.start
        self.height = height
        self.color = common.Color( color, 1 )
        self.x, self.y  = self.__set_gene_coord( tail_length )
        
    def __set_gene_coord( self, tail_length ):
        decimal_place = 3 #小数点3位
        x1 = self.start
        x3 = self.end
        y1 = round( self.y_mid_coord - self.height / 2, decimal_place )
        #y1 = self.y_mid_coord - self.height * 0.8
        y2 = round( self.y_mid_coord, decimal_place )
        y3 = round( self.y_mid_coord + self.height / 2, decimal_place )
        #y3 = self.y_mid_coord + self.height * 0.8
        if( self.strand == '+' ): # => or >
            if self.width > tail_length :
                x2 = self.end - tail_length
                x = [ x1, x2, x3, x2, x1 ]
                y = [ y1, y1, y2, y3, y3 ]
            else:
                x = [ x1, x3, x1 ]
                y = [ y1, y2, y3 ]
        elif( self.strand == '-' ): # <= or <
            if self.width > tail_length:
                x2 = self.start + tail_length
                x = [ x1, x2, x3, x3, x2 ]
                y = [ y2, y1, y1, y3, y3 ]
            else:
                x = [ x1, x3, x3 ]
                y = [ y2, y1, y3 ]
        else: 
            x = [ x1, x3, x3, x1 ]
            y = [ y1, y1, y3, y3 ]
        return x, y

    def __convert_position2coord( self, scaffold, buf ):
        start = scaffold.convert_position2xcoord( int( buf[3] ) )
        end = scaffold.convert_position2xcoord( int( buf[4] ) )
        return start, end, buf[6]

    def plot( self, ax ):
        ax.fill( self.x, self.y, color=self.color.color, alpha=self.color.alpha, lw=0.25, edgecolor='black' )


def plot_genes( seq, ax, size, fn, color, gene_font_size ):
    RATIO = size.gene_ratio
    color = color.color

    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            if( buf[2] != 'gene' ):
                continue
            i = common.detect_index( buf[0], int( buf[3] ), int( buf[4] ), seq )
            if( i == -1 ):
                continue

            height = seq.height * RATIO
            y_total = (( len( seq.pos ) -1 ) * seq.MARGIN ) + size.top_margin + size.bottom_margin 
            tail_length = ( height / 3 ) * ( seq.width / y_total ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
            if( len( buf ) == 10 ):
                color1 = common.Color( buf[9], 1 )
                aninstance = Gene( seq, buf, height, tail_length, color1.color )
            else:
                aninstance = Gene( seq, buf, height, tail_length, color )
            aninstance.plot( ax )

            #遺伝子名の表示
            if gene_font_size == 0:
                continue
            if re.search( r'Name=(\S+?);', buf[8] ):
                gene_name = re.search( r'Name=(\S+?);', buf[8] ).group(1)
                text_x = seq.convert_position2xcoord( ( int( buf[3] )*3/7 + int( buf[4] )*4/7 ) )
                text_y = seq.convert_position2ycoord( ( int( buf[3] )*3/7 + int( buf[4] )*4/7 ), -1 * height/2 - 0.04 )
                ax.text( text_x, text_y, gene_name, fontsize = gene_font_size, color = 'black', ha='center', va='top', rotation=75, style='italic' )


def plot_genes_excel( seq, ax, size, fn, color, gene_font_size ):
    RATIO = size.gene_ratio
    color = color.color

    df = pd.read_excel( fn, header=None )


    for index, row in df.iterrows():
        buf = row.tolist()
        if( buf[2] != 'gene' ):
            continue
        i = common.detect_index( buf[0], int( buf[3] ), int( buf[4] ), seq )
        if( i == -1 ):
            continue

        height = seq.height * RATIO
        y_total = (( len( seq.pos ) -1 ) * seq.MARGIN ) + size.top_margin + size.bottom_margin
        tail_length = ( height / 3 ) * ( seq.width / y_total ) * ( size.figsize_inch[1] / size.figsize_inch[0] )
        if( len( buf ) == 10 ):
            color1 = common.Color( buf[9], 1 )
            aninstance = Gene( seq, buf, height, tail_length, color1.color )
        else:
            aninstance = Gene( seq, buf, height, tail_length, color )
        aninstance.plot( ax )

        #遺伝子名の表示
        if gene_font_size == 0:
            continue
        if re.search( r'Name=(\S+?);', buf[8] ):
            gene_name = re.search( r'Name=(\S+?);', buf[8] ).group(1)
            text_x = seq.convert_position2xcoord( ( int( buf[3] )*3/7 + int( buf[4] )*4/7 ) )
            text_y = seq.convert_position2ycoord( ( int( buf[3] )*3/7 + int( buf[4] )*4/7 ), -1 * height/2 - 0.04 )
            ax.text( text_x, text_y, gene_name, fontsize = gene_font_size, color = 'black', ha='center', va='top', rotation=75, style='italic' )


def cal_intergenic( fn, seq ):
    intergenic = []
    end = seq.start
    with open( fn , 'r' ) as file:
        for line in file:
            buf = line.rstrip( '\n' ).split( '\t' )
            if( len( buf ) < 4 ):
                continue
            if( buf[2] != 'gene' and buf[2] != 'pseudogene' ):
                continue

            i = common.detect_index( buf[0], int( buf[3] ), int( buf[4] ), seq )
            if( i == -1 ):
                continue
            start = int( buf[3] )
            if( start - end >= 10 ):
                # intergenic.append( int(( start + end )/2) )
                intergenic.append( end + 5 ) # 2024/10/16 
            end = int( buf[4] )
    return intergenic


def cal_intergenic_excel( fn, seq ):
    intergenic = []
    end = seq.start

    df = pd.read_excel( fn )

    for index, row in df.iterrows():
        buf = row.tolist()
        if( len( buf ) < 4 ):
            continue
        if( buf[2] != 'gene' and buf[2] != 'pseudogene' ):
            continue

        i = common.detect_index( buf[0], int( buf[3] ), int( buf[4] ), seq )
        if( i == -1 ):
            continue
        start = int( buf[3] )
        if( start - end >= 10 ):
            intergenic.append( int(( start + end )/2) )
        end = int( buf[4] )

    return intergenic

