import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import csv
import os
import sys
import math
from traincars import common

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

class Scaffold:
    
    def __init__( self, ID, start, end, name, width, row_count ):
        self.ID = ID
        self.start = start
        self.end = end
        self.name = name
        self.length = self.end - self.start
        self.MARGIN = 1
        ( self.width, self.row_count ) = self.__set_width_and_row( width, row_count )
        self.pos = self.__set_positions()
        self.height = 0.03
        color = '#111111'
        #color = '#95949a'
        alpha = 1
        self.color = common.Color( color, alpha )

    def __set_positions( self ):
        positions = [];
        remaining_length = self.length
        x = 0
        y = 0
        pos = self.start
        while( remaining_length >= 0 ):
            if( remaining_length - self.width >= 0 ):
                positions.append( {'x': x, 'y' : y, 'pos': pos, 'length': self.width } )
                remaining_length -= self.width
                pos += self.width
            else:
                positions.append( {'x': x, 'y' : y, 'pos': pos, 'length': remaining_length } )
                break
            y -= self.MARGIN

        #print( positions )
        return positions

    
    def __set_width_and_row( self, width, row_count ):
        if( width == 0 ):
            tmp = ( self.length / row_count )
            if 6 <= math.log10( tmp ):
                width = math.ceil( tmp/math.pow( 10, 6 ) ) * math.pow( 10, 6 );
            elif 3 <= math.log10( tmp ):
                width = math.ceil( tmp/math.pow( 10, 3 ) ) * math.pow( 10, 3 );
        else:
            row_count = math.ceil( self.length / width )

        return width, row_count

    def correct_width_by_genes( self, intergenic ):
        l_count = 0
        s = 0
        for i in range( s, len( intergenic )):
            if l_count == len( self.pos ):
                break
            if intergenic[i] < self.pos[ l_count ][ 'pos' ] + self.pos[ l_count ][ 'length' ] :
                continue
            s = i - 1
            self.pos[ l_count ][ 'length' ] =  intergenic[i-1] - self.pos[ l_count ]['pos']
            diff = self.width - self.pos[ l_count ][ 'length' ]
            l_count += 1
            for j in range( l_count, len( self.pos ) ):
                self.pos[ j ][ 'pos' ] -= diff
            self.pos[ -1 ][ 'length' ] += diff
            if self.pos[-1]['length'] > self.width : #はみ出たから新しい行を作成する
                x = 0
                y = self.pos[-1]['y'] - 1
                self.pos[-1]['length'] = self.width
                pos = self.pos[-1]['pos'] + self.pos[-1]['length']
                length = self.length - ( self.pos[-1]['pos'] + self.pos[-1]['length'] - self.start )
                self.pos.append( {'x': x, 'y' : y, 'pos': pos, 'length': length } )
            self.row_count = len( self.pos )
                
    def convert_position2xcoord( self, pos ):
        x = pos - self.start
        for i in range( len( self.pos )):
            if self.pos[i]['length'] <= x :
                x -= self.pos[i]['length']
            else:
                break
        return x

    
    def convert_position2ycoord( self, pos, h ):
        l = pos - self.start
        for i in range( len( self.pos )):
            if self.pos[i]['length'] <= l :
                l -= self.pos[i]['length']
            else:
                y = self.pos[i]['y']
                break
        return y + h

    
    def plot_line( self, ax ):
        slen = self.length
        for i in range( len( self.pos )):
            y = [self.pos[i]['y'], self.pos[i]['y'], self.pos[i]['y'] + self.height , self.pos[i]['y'] + self.height ]            
            x = [self.pos[i]['x'], self.pos[i]['x'] + self.pos[i]['length'], self.pos[i]['x'] + self.pos[i]['length'], self.pos[i]['x'] ]
            ax.fill( x, y, color=self.color.color, lw=0, alpha=self.color.alpha )
            slen -= self.width

            
    def output_parameters( self ):
        print( '##Scaffold paramenters:' )
        print( '  seq ID: %s' % ( self.ID ))
        print( '  start position: %d' % ( self.start ))
        print( '  end position: %d' % ( self.end ))
        print( '  seq length: %d' % ( self.length ))
        print( '  set width: %d' % ( self.width ))
        for i in range( len( self.pos )):
            print( '  actial width: %d' % ( self.pos[i]['length'] ))
        print( '  row count: %d' % ( self.row_count ))
        print( '' )
            

def input_scaffoldtsv( fn, width, row_count ):        
    try:
        # ファイルの存在を確認
        if not os.path.exists(fn):
            raise FileNotFoundError(f"Error: The file '{fn}' does not exist. Please check the file path.")
    
        # ファイルが空かどうかを確認
        if os.path.getsize(fn) == 0:
            raise ValueError(f"Error: The file '{fn}' is empty. Please provide a valid input file.")

        with open( fn , 'r' ) as file:
            buf = csv.reader( file, delimiter = '\t' )
            for scaffold, scaf_start, scaf_end, scaf_name in buf:
                seq = Scaffold( scaffold, int( scaf_start ), int( scaf_end ), scaf_name, width, row_count )
        return seq

    except (FileNotFoundError, ValueError) as e:
        print( e, file=sys.stderr )
        sys.exit(1)


def plot_scaffolds( ax, seq, s_font_size ):
    ##plot scaffold
    seq.plot_line( ax )
    ##plot scaffold name
    x = seq.width / 2 
    y = 1
    aninstance = common.Text( seq.name, s_font_size, x, y, 'center', 'top' )
    aninstance.output( ax )

            
class Scalebar:
    def __init__( self, size, seq ):
        self.width = self.__set_width( seq.width )
        self.bar_color = 'black'
        self.legend = self.__set_label()
        self.origin_x, self.origin_y = self.__set_origin( seq )
        self.label = common.Text( self.legend, 8, self.origin_x + self.width/2, self.origin_y - 0.4/size.figsize_inch[1], 'center', 'top' )
        
    def __set_label( self ):
        if 6 <= math.log10( self.width ):
            return str( int( self.width/ math.pow( 10, 6 )) ) + ' Mbp'
        if 3 <= math.log10( self.width ):
            return str( int( self.width/ math.pow( 10, 3 )) ) + ' kbp'
        return str( int( self.width )) + ' bp'

    def __set_origin( self, seq ):
        x = 0
        y = len( seq.pos ) * seq.MARGIN * -1 + 0.25
        return x, y
        
    def __set_width( self, width ):
        if  math.pow( 10, int( math.log10( width )))/ width <= 0.2 :
            return math.pow( 10, int( math.log10( width )))
        if math.pow( 10, int( math.log10( width )))/ width <= 0.4 :
            return math.pow( 10, int( math.log10( width )))/2
        return math.pow( 10, int( math.log10( width )))/5

    def plot_legend( self, ax ):
        ax.plot( [ self.origin_x, self.origin_x + self.width ], [ self.origin_y, self.origin_y ], color=self.bar_color, lw=0.5 )
        self.label.output( ax )

    def plot_tick( self, ax, seq, size ):
        if 6 <= math.log10( seq.end ):
            mr = 10 ** 6
            unit = 'Mbp'
        elif 3 <= math.log10( seq.end ):
            mr = 10 ** 3
            unit = 'kbp'
        else:
            mr = 1
            unit = 'bp'
        start = int( math.ceil( seq.start / self.width ) * self.width )

        for i in range( start, seq.end, int( self.width ) ):
            x = seq.convert_position2xcoord( i )
            y1 = seq.convert_position2ycoord( i, ( 1 + size.gene_ratio/2 ) * seq.height )
            y2 = seq.convert_position2ycoord( i, ( 2 + size.gene_ratio/2 ) * seq.height )
            y3 = seq.convert_position2ycoord( i, ( 2.25 + size.gene_ratio/2 ) * seq.height )
            ax.plot( [x, x], [y1, y2], color=self.bar_color, lw=0.25 )
            if mr / self.width <= 1 :
                aninstance = common.Text( "{0} {1}".format( int( i/mr ), unit ), 5, x, y3, 'center', 'bottom' )
            elif mr / self.width <= 10 :
                aninstance = common.Text( "{0:0.1f} {1}".format( i/mr, unit ), 5, x, y3, 'center', 'bottom' )
            elif mr / self.width <= 100 :
                aninstance = common.Text( "{0:0.2f} {1}".format( i/mr, unit ), 5, x, y3, 'center', 'bottom' )
            else :
                #aninstance = common.Text( "{0:0.3f} {1}".format( i/mr, unit ), 5, x, y3, 'center', 'bottom' )
                if unit == 'Mbp' :
                    aninstance = common.Text( "{0:,} {1}".format( int(i*1000/mr), 'kbp' ), 5, x, y3, 'center', 'bottom' )
                else :
                    aninstance = common.Text( "{0:,} {1}".format( int(i*1000/mr), 'bp' ), 5, x, y3, 'center', 'bottom' )

            aninstance.output( ax )

            
    def output_parameters( self ):
        print( '##Scalebar paramenters:' )
        print( '  width: %d' % ( self.width ))
        print( '  origin_x %.2f' % ( self.origin_x ))
        print( '  origin_y: %.2f' % ( self.origin_y ))
        print( '' )

