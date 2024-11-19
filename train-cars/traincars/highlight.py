import matplotlib.pyplot as plt
import os
import sys
import csv
from traincars import common

class Highlight:

    def __init__( self, scaffold, buf, height, color ):
        self.height = height
        self.color = color
        self.pos = self.__convert_position2coord( scaffold, buf )


    def __convert_position2coord( self, scaffold, buf ):
        start = int( buf[1] )
        pos = []
        for i in range( len( scaffold.pos )):
            if scaffold.pos[i]['pos'] + scaffold.pos[i]['length'] < start:
                continue
            if int( buf[2] ) <= scaffold.pos[i]['pos'] + scaffold.pos[i]['length']: #1段で済んだ場合
                end = int( buf[2] )
                x1 = scaffold.convert_position2xcoord( start )
                x2 = scaffold.convert_position2xcoord( end )
                y1 = scaffold.convert_position2ycoord( end, scaffold.height /2 + self.height / 2 )
                y2 = scaffold.convert_position2ycoord( end, scaffold.height /2 - self.height / 2 )
                x = [ x1, x2, x2, x1 ]
                y = [ y1, y1, y2, y2 ]
                pos.append( { 'x': x, 'y': y } )
                break
            else : #またがる場合
                end = scaffold.pos[i]['pos'] + scaffold.pos[i]['length'] -1
                x1 = scaffold.convert_position2xcoord( start )
                x2 = scaffold.convert_position2xcoord( end )
                y1 = scaffold.convert_position2ycoord( end, scaffold.height /2 + self.height / 2 )
                y2 = scaffold.convert_position2ycoord( end, scaffold.height /2 - self.height / 2 )
                x = [ x1, x2, x2, x1 ]
                y = [ y1, y1, y2, y2 ]
                pos.append( { 'x': x, 'y': y } )
                start = scaffold.pos[i+1]['pos']
        return pos

    def plot( self, ax ):
        for i in range( len( self.pos ) ):
            ax.fill( self.pos[i]['x'], self.pos[i]['y'], color=self.color.color, alpha=self.color.alpha, lw=0, zorder=0 )

def plot_highlight( seq, ax, size, fn ):
    RATIO = size.gene_ratio + 2.4
    try:
        # ファイルの存在を確認
        if not os.path.exists(fn):
            raise FileNotFoundError(f"Error: The file '{fn}' does not exist. Please check the file path.")
        
        # ファイルが空かどうかを確認
        if os.path.getsize(fn) == 0:
            raise ValueError(f"Error: The file '{fn}' is empty. Please provide a valid input file.")
        
        with open( fn , 'r' ) as file:
            for line in file:
                buf = line.rstrip( '\n' ).split( '\t' )
                i = common.detect_index( buf[0], int( buf[1] ), int( buf[2] ), seq )
                if( i == -1 ):
                    continue
                color1 = common.Color( buf[3], 0.3 )
                height = seq.height * RATIO
                aninstance = Highlight( seq, buf, height, color1 )
                aninstance.plot( ax )

    except (FileNotFoundError, ValueError) as e:
        # エラーメッセージを標準エラー出力に表示
        print(e, file=sys.stderr)
