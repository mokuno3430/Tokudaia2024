import sys
import argparse
import os.path
import math
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from matplotlib.backends.backend_pdf import PdfPages


def main():
    args = func_get_args()
    func_print_messages( )
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    #plt.rcParams['font.family'] = 'Arial'

    #blastのアライメント結果の読み込み
    blastn = func_input_blastn6( args.blastn, args.min_identity, args.min_alignlen )

    #ラベル名と配列長の情報の読み込み
    table_db = func_input_list( args.input1, blastn, 'slen' )
    table_qry = func_input_list( args.input2, blastn, 'qlen' )

    #subplotのサイズ、軸の設定など
    fig, ax = func_set_subplots_parameters( table_db, table_qry, args )

    #blastnの結果を可視化
    heatmap = Colormap( args.min_identity, 100, args.colormap )
    func_plot_alignments( blastn, table_db, table_qry, heatmap, ax, args.line_width )
    colormap_legend = Colorbox( heatmap, ax, args.figure_size )
    colormap_legend.plot( ax, heatmap, fig, args.tick_label_size )

    #highlightの表示
    if args.highlight != None:
        highlight = func_input_tsv_for_highlight( args.highlight )
        if highlight is not None:
            func_plot_highlights( highlight, table_db, table_qry, ax, args.h_alpha )

    #highlight crossedの表示
    if args.highlight_crossed != None:
        highlight_c = func_input_tsv_for_highlight( args.highlight_crossed )
        if highlight_c is not None:
            func_plot_crossed( highlight_c, table_db, table_qry, ax, args.h_alpha )

    #結果の保存
    pdf_file = args.out + '.pdf'
    pp = PdfPages( pdf_file )
    pp.savefig( fig )
    pp.close()
    plt.clf()
    print(f"{pdf_file} has been generated.")


class Colormap:
    pallet = {
        'red'  :((0.00, 0.10, 0.10), (0.25, 0.10, 0.10), (0.40, 0.30, 0.30), (0.60, 1.00, 1.00), (0.80, 0.90, 0.90), (1.00, 0.70, 0.70)),
        'green':((0.00, 0.10, 0.10), (0.25, 0.60, 0.60), (0.40, 0.80, 0.80), (0.60, 0.75, 0.75), (0.80, 0.30, 0.30), (1.00, 0.15, 0.15)),
        'blue' :((0.00, 0.40, 0.40), (0.25, 1.00, 1.00), (0.40, 0.25, 0.25), (0.60, 0.00, 0.00), (0.80, 0.05, 0.05), (1.00, 0.20, 0.20))
    }
    mycm = clr.LinearSegmentedColormap('original', pallet )
    cmaps = [ cm.binary, cm.bone_r, cm.inferno_r, cm.hot_r, cm.YlGnBu, mycm ]
    cmap_list = [ 'binary', 'bone_r', 'inferno_r', 'hot_r', 'YlGnBu', 'original' ]

    def __init__( self, min_identity, max_identity, cm, alpha=1 ):
        self.min_identity = min_identity
        self.max_identity = max_identity
        self.alpha = alpha
        self.cm = cm

    def convert_identity2color( self, identity ):
        color_value=( identity - self.min_identity )/( self.max_identity - self.min_identity )
        align_color=Colormap.cmaps[ self.cm ]( color_value )
        return align_color

    def output_parameters( self ):
        print( '##Colormap paramenters:' )
        print( '  min_identity: %d' % ( self.min_identity))
        print( '  max_identity: %d' % ( self.max_identity ))
        print( '  alpha: %.2f' % ( self.alpha ))
        print( '  colormap: %d (%s)' % ( self.cm, Colormap.cmap_list[self.cm]  ))
        print( '' )


class Colorbox:
    def __init__( self, heatmap, ax, figure_size ):
        if heatmap.max_identity - heatmap.min_identity <= 2:
            self.BIN = 0.025
        elif heatmap.max_identity - heatmap.min_identity <= 5:
            self.BIN = 0.1
        elif heatmap.max_identity - heatmap.min_identity <= 10:
            self.BIN = 0.2
        else:
            self.BIN = 0.5

        self.height = 0.015 if figure_size[0] <= figure_size[1] else 0.015 * figure_size[0] / figure_size[1]
        pos = ax[-1][-1].get_position()
        self.width = pos.x1 - pos.x0 if len( ax[0] ) >= 3 else 0.20
        self.cell_width = self.width * self.BIN /( heatmap.max_identity - heatmap.min_identity + self.BIN )
        self.origin_x = pos.x1 - self.width
        self.origin_y = 0.025

    def plot( self, ax, heatmap, fig, font_size ):
        ax2 = ax[-1][0]
        cell_originx = self.origin_x
        SCALE = 5
        if heatmap.max_identity - heatmap.min_identity <= 1 : #identity 99%以上を表示する時は0.2刻み
            SCALE = 0.25
        elif heatmap.max_identity - heatmap.min_identity <= 2 : #identity 98%以上を表示する時は0.5刻み
            SCALE = 0.5
        elif heatmap.max_identity - heatmap.min_identity <= 5: #identity 95%以上を表示する時は１刻み
            SCALE = 1
        elif heatmap.max_identity - heatmap.min_identity <= 10: #identity 90%以上の時は2刻み
            SCALE = 2
        #colormapの凡例を表示する
        ax2.text(self.origin_x, self.origin_y + self.height/2, 'identity (%)  ', va='center', ha='right', transform=fig.transFigure, fontsize=font_size *1.5 )
        for i in range( int(heatmap.min_identity/self.BIN), int(heatmap.max_identity/self.BIN) + 1 ):
            align_color=heatmap.convert_identity2color( i*self.BIN )
            colorbox_x = [ cell_originx, cell_originx + self.cell_width, cell_originx + self.cell_width, cell_originx ]
            colorbox_y = [ self.origin_y, self.origin_y, self.origin_y + self.height, self.origin_y + self.height ]
            cell_originx += self.cell_width
            # color cellのプロット
            ax2.fill( colorbox_x, colorbox_y, color=align_color, alpha=heatmap.alpha, linewidth=0, transform=fig.transFigure, clip_on=False )
            # 数字のプロット
            if( i*self.BIN % SCALE == 0 ):
                if SCALE < 1:
                    ax2.text( cell_originx - self.cell_width/2, self.origin_y-0.004, i*self.BIN, fontsize = font_size, color = 'black', ha='center', va='top', transform=fig.transFigure )
                else:
                    ax2.text( cell_originx - self.cell_width/2, self.origin_y-0.004, int( i*self.BIN ), fontsize = font_size, color = 'black', ha='center', va='top', transform=fig.transFigure )

    def output_parameters( self ):
        print( '##Colorbox paramenters:' )
        print( '  width: %.2f' % ( self.width ))
        print( '  height: %.2f' % ( self.height ))
        print( '  origin_x %.2f' % ( self.origin_x ))
        print( '  origin_y: %.2f' % ( self.origin_y ))
        print( '' )


def func_get_args():
    parser = argparse.ArgumentParser( formatter_class=argparse.MetavarTypeHelpFormatter )
    parser.add_argument( '-i1', '--input1', help='sequence IDs of database at blastn (=row)', type=str, required=True )
    parser.add_argument( '-i2', '--input2', help='sequence IDs of query at blastn (=column)', type=str, required=True )
    parser.add_argument( '--blastn', help='blastn.tsv (-outfmt \'6 std qlen slen\' at blastn)', type=str, required=True )
    parser.add_argument( '--highlight', help='optional: positions.tsv (scaffold start end color)', type=str )
    parser.add_argument( '--highlight_crossed', help='optional: positions.tsv (scaffold start end color)', type=str )
    parser.add_argument( '--out', help='optional: prefix of pdf file (default out)', type=str, default='out' )
    parser.add_argument( '--min_identity', help='optional: minimum sequence identity (default 90)', type=int, default=90 )
    parser.add_argument( '--min_alignlen', help='optional: minimum alignment length (default 100)', type=int, default=100 )
    parser.add_argument( '--line_width', help='optional: line width of dotplots (default 1.0)', type=float, default=1.0 )
    parser.add_argument( '--share', help='optional: sharing axis scales among subplots.', action='store_true' )
    parser.add_argument( '--show_grid', help='optional: show grid-line', action='store_true' )
    parser.add_argument( '--xtitle_rotate', help='optional: (default 0)', type=int, default=0 )
    parser.add_argument( '--ytitle_rotate', help='optional: (default 90)', type=int, default=90 )
    parser.add_argument( '--font_size', help='optional: (default 8)', type=float, default=8 )
    parser.add_argument( '--figure_size', help='optional: (default [8,8])', type=float, nargs='*', default=[8,8] )
    parser.add_argument( '--hspace', help='optional: (default -1 means auto without --share, 0.04 with --share)', type=float, default=-1 )
    parser.add_argument( '--wspace', help='optional: (default -1 means auto without --share, 0.04 with --share)', type=float, default=-1 )
    parser.add_argument( '--h_alpha', help='optional: transparency ratio of highlights (default 0.3)', type=float, default=0.3 )
    parser.add_argument( '--tick_label_size', help='optional: (default 6)', type=float, default=6 )
    parser.add_argument( '--tick_width', help='optional: scale width of axis (bp) (default -1 means auto)', type=int, default=-1 )
    parser.add_argument( '--colormap', help='optional: colormap for identity of alignments ( 0:binary, 1:bone_r, 2:inferno_r, 3:hot_r, 4:YlGnBu, 5:original) (default 5)', choices=[ 0, 1, 2, 3, 4, 5 ], default=5, type=int )
    parser.add_argument( '-v', '--version', action='version', version='%(prog)s 4.0', default=False)
    args = parser.parse_args()
    return ( args )


def func_set_subplots_parameters( table_db, table_qry, args ):
    share = args.share
    figure_size = args.figure_size
    #それぞれのsubplotのサイズを配列長に調整
    #--share y軸目盛を一番左のplotのみ, x軸目盛を一番下のplotのみにする場合に使用
    if not share:
        ws = 0.15 if args.wspace == -1 else args.wspace
        hs = 0.15 if args.hspace == -1 else args.hspace
        fig, ax  = plt.subplots( len(table_db), len(table_qry), figsize=(figure_size[0], figure_size[1]),
                                gridspec_kw=dict(width_ratios=table_qry['len'].tolist(),
                                    height_ratios=table_db['len'].tolist() ))
        fig.subplots_adjust(hspace=hs, wspace=ws, left=0.15, right=0.98, top=0.90, bottom=0.07)
    else:
        ws = 0.04 if args.wspace == -1 else args.wspace
        hs = 0.04 if args.hspace == -1 else args.hspace
        fig, ax  = plt.subplots( len(table_db), len(table_qry), figsize=(figure_size[0], figure_size[1]),sharex='col', sharey='row',
                                gridspec_kw=dict(width_ratios=table_qry['len'].tolist(),
                                    height_ratios=table_db['len'].tolist() ))
        # fig.subplots_adjust(hspace=hs, wspace=ws, left=0.15, right=0.98, top=0.90, bottom=0.07)
        fig.subplots_adjust(hspace=hs, wspace=ws )
    #どんなsubplotsも擬似的に2次元配列として扱うための操作
    if isinstance(ax, plt.Axes):
        ax = np.array([[ax]])
    else:
        ax = ax.reshape( len(table_db), len( table_qry) )
    #subplotの大きさ、目盛, grid線などの設定
    fig, ax = func_for_multi_plots( fig, ax, table_db, table_qry, args )
    return fig, ax


def func_for_multi_plots( fig, ax, table_db, table_qry, args ):
    x_rotate = args.xtitle_rotate
    y_rotate = args.ytitle_rotate
    grid = args.show_grid
    title_font_size = args.font_size
    tick_label_size = args.tick_label_size
    max_width = table_db['len'].max() if table_db['len'].max() > table_qry['len'].max() else table_qry['len'].max()
    sum_width = table_db['len'].sum() if table_db['len'].sum() > table_qry['len'].sum() else table_qry['len'].sum()
    mr, unit, width = func_set_tick_label( max_width, sum_width )
    width = width if args.tick_width == -1 else args.tick_width
    # 軸の範囲の指定など
    for i in range( 0, len( table_db )):
        strain_y=table_db.at[ i, 'label']
        y_max = table_db.at[ i, 'len']
        for j in range( 0, len( table_qry ) ):
            strain_x=table_qry.at[ j, 'label']
            x_max = table_qry.at[j, 'len']
            if table_qry.at[j, 'strand'] == '+':
                ax[i][j].set_xlim( 0, x_max )
            else:
                ax[i][j].set_xlim( x_max, 0 )
            if table_db.at[i, 'strand'] == '+':
                ax[i][j].set_ylim( 0, y_max )
            else:
                ax[i][j].set_ylim( y_max, 0 )
            ax[i][j].spines["right"].set_linewidth(0.5) #線の太さは0.5
            ax[i][j].spines["left"].set_linewidth(0.5) #線の太さは0.5
            ax[i][j].spines["bottom"].set_linewidth(0.5) #線の太さは0.5
            ax[i][j].spines["top"].set_linewidth(0.5) #線の太さは0.5

            ax[i][j].xaxis.set_tick_params( width=0.5 )
            ax[i][j].yaxis.set_tick_params( width=0.5 )
            if grid: #grid線を引く
                ax[i][j].grid( axis='both', color='grey', lw=0.3, linestyle='dotted' )
            ax[i][j].set_aspect('equal', adjustable='box') #縦横比が1:1になるようにする

            # 一番左のdotplotsにtitleを表示
            if( j == 0 ):
                ypad = 10 if y_rotate == 90 else 30 #y titleの位置を設定するための小手先の操作
                ax[i][j].set_ylabel( strain_y, fontsize=title_font_size, labelpad=ypad, rotation=y_rotate, va='center' )

            y_spos = table_db.at[i, 'start_pos']
            s_ytick = (( int( y_spos / width ) + 1) * width) - y_spos
            x_spos = table_qry.at[j, 'start_pos']
            s_xtick = (( int( x_spos / width ) + 1) * width) - x_spos

            # y軸目盛の表示
            ax[i][j].tick_params( axis='y', which='both', length=1, pad=1 )
            ax[i][j].set_yticks( np.arange( s_ytick, y_max, width ) )
            # 一番上のdotplotsにtitleを表示
            if( i == 0 ):
                ax[i][j].set_title( strain_x, fontsize=title_font_size, rotation=x_rotate ) #rotateで角度調整
            #x軸の目盛を表示
            ax[i][j].tick_params( axis='x', which='both', length=1, pad=1 )
            ax[i][j].set_xticks( np.arange( s_xtick, x_max, width ) )

            if width % mr == 0:
                ax[i][j].set_yticklabels( (np.arange( ( y_spos + s_ytick), (y_spos + y_max), width)/mr).astype(int), fontsize=tick_label_size )
                ax[i][j].set_xticklabels( (np.arange( ( x_spos + s_xtick), (x_spos + x_max), width)/mr).astype(int), fontsize=tick_label_size )
            else:
                ax[i][j].set_yticklabels( np.arange( ( y_spos + s_ytick), (y_spos + y_max), width )/mr , fontsize=tick_label_size )
                ax[i][j].set_xticklabels( np.arange( ( x_spos + s_xtick), (x_spos + x_max), width )/mr, fontsize=tick_label_size )

    #左上に単位を表示する
    ax[0][0].text( 0, 1.02, unit, fontsize = tick_label_size, color = 'black', ha='right', va='bottom', transform=ax[0][0].transAxes )
    return fig, ax


def func_set_tick_label( max_len, sum_len ):
    mr = 1 #3桁区切り
    unit = '(bp)' #単位
    if 9 <= math.log10( max_len ):
        mr = 10 ** 9
        unit = '(Gbp)'
    elif 6 <= math.log10( max_len ):
        mr = 10 ** 6
        unit = '(Mbp)'
    elif 3 <= math.log10( max_len ):
        mr = 10 ** 3
        unit = '(kbp)'
    #要はいい感じの目盛幅にしたい 目盛の数字は1、２, 5, 10, 20, 50, 100, 200, 500のどれか
    width=int( math.pow( 10, int( math.log10( max_len )))/5 )
    if math.pow( 10, int( math.log10( max_len )))/ max_len <= 0.2:
        width=int( math.pow( 10, int( math.log10( max_len ))) )
    elif math.pow( 10, int( math.log10( max_len )))/ max_len <= 0.4:
        width=int( math.pow( 10, int( math.log10( max_len )))/2)
    return( mr, unit, width )


def func_plot_alignments( blastn, table_row, table_column, heatmap, ax, line_width ):
    #アライメント結果のplot
    for b_out in blastn.itertuples():
        if b_out.sseqid not in table_row['scaffold'].values:
            continue
        if b_out.qseqid not in table_column['scaffold'].values:
            continue
        c_value = heatmap.convert_identity2color( float( b_out.identity ))
        i = table_row.index[ table_row['scaffold'] == b_out.sseqid ][0]
        j = table_column.index[ table_column['scaffold'] == b_out.qseqid ][0]
        ax[ i ][ j ].plot( [ b_out.qstart, b_out.qend ], [ b_out.sstart, b_out.send ],color=c_value, lw=line_width )


def func_input_blastn6( fn, min_identity, min_alignlength ):
    required_columns = [0, 1, 2, 6, 7, 8, 9, 12, 13]
    try:
        if not os.path.exists(fn): #ファイルがない場合
            raise FileNotFoundError(f"Error: {fn} does not exist.")

        if os.path.getsize(fn) == 0 : #ファイルが空の場合
            raise ValueError(f"Error: {fn} is empty.")

        df = pd.read_csv( fn, sep='\t', header=None )
        if df.shape[1] < 14:  # 要求された列数に足りない場合
            raise ValueError(
                "Error: The input file does not contain the required columns.\n"
                "Please run blastn with the option \"-outfmt '6 std qlen slen'\" to generate the correct output format."
            )

        df = df[required_columns]  # 要求された列だけを抽出
        df.columns = ['qseqid', 'sseqid', 'identity', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen' ]

        df['qseqid'] = df['qseqid'].astype(str)
        df['sseqid'] = df['sseqid'].astype(str)

        # identityとalignment lengthでフィルタリング
        blastn = df[ (df['identity'] >= min_identity) & (df['qend'] - df['qstart'] + 1 >= min_alignlength ) ]

        if blastn.empty:
            raise ValueError(
                f"Error: No alignments available for visualization in {fn}.\n"
                "Please check the input file and consider adjusting the values for --min_identity or --min_alignlen."
            )
        return blastn

    except (FileNotFoundError, ValueError) as e:
        # エラーメッセージを表示してプログラムを中断
        print( e, file=sys.stderr )
        sys.exit( 1 )
        # raise  # エラーを再度発生させ、プログラム全体の中断を促す


def func_input_list( fn, blastn, len_ID ):
    try:
        if not os.path.exists(fn):
            raise FileNotFoundError(f"Error: {fn} does not exist.")

        if os.path.getsize(fn) == 0:
            raise ValueError(f"Error: {fn} is empty.")

        table = pd.read_csv( fn, sep='\t', header=None )
        num_columns = table.shape[1]
        scaffold_col = 'scaffold'
        label_col = 'label'
        strand_col = 'strand'
        start_pos_col = 'start_pos'
        if( num_columns == 1 ):
            table.columns = [scaffold_col]
            table[label_col] = table[ scaffold_col ]
            table[strand_col] = '+'
            table[start_pos_col] = 0
        elif( num_columns == 2 ):
            table.columns = [scaffold_col, label_col]
            table[strand_col] = '+'
            table[start_pos_col] = 0
        else:
            unique_strand_values = set(table[2].unique())
            if unique_strand_values.issubset({'+', '-'}):
                table.columns = [scaffold_col, label_col, strand_col]
                table[start_pos_col] = 0
            else:
                # 一時的な対応: strand情報がない場合
                table.columns = [scaffold_col, label_col, start_pos_col]
                table[strand_col] = '+'

        table['scaffold'] = table['scaffold'].astype(str)
        #blastnの結果から配列長を取り出す
        if len_ID == 'qlen':
            lens = blastn.set_index('qseqid')['qlen'].to_dict()
        else:
            lens = blastn.set_index('sseqid')['slen'].to_dict()
        #配列長のデータを追加
        table["len"] = table["scaffold"].map( lens )
        #blastにhitがない配列はlengthの情報が取り出せないので除外
        table = table.dropna(subset=['len']).reset_index(drop=True)
        if table.empty:
            raise ValueError( "The DataFrame is empty. Please check your sequence IDs in {}.".format( fn ) )
        return table
    
    except (FileNotFoundError, ValueError) as e:
        # エラーメッセージを表示してプログラムを中断
        print(e, file=sys.stderr)
        sys.exit(1)


def func_input_tsv_for_highlight( fn ):
    try:
        if not os.path.exists(fn):
            raise FileNotFoundError(f"Error: {fn} does not exist. Calculation will continue.")
        if os.path.getsize(fn) == 0:
            raise ValueError(f"Error: {fn} is empty. Calculation will continue.")

        df = pd.read_csv( fn, sep='\t', header=None )
        if df.shape[1] != 4:
            raise ValueError("Error: The file should have at four columns: contig_ID, start_pos, end_pos and colorcode. Calculation will continue.")
        df.columns = ['seqid', 'start', 'end', 'color' ]
        df['seqid'] = df['seqid'].astype(str)
        return( df )
    except (FileNotFoundError, ValueError) as e:
        # エラーメッセージを表示
        print( e, file=sys.stderr )


def func_plot_highlights( highlight, table_row, table_column, ax, alpha ):
    alpha = alpha
    for buf in highlight.itertuples():
        if not clr.is_color_like( buf.color ):
            print( f"Error: {buf.color} is invalid color. Calculation will continue.", file=sys.stderr )
            continue
        #y軸scaffoldへのhighlight
        if buf.seqid in table_row['scaffold'].values:
            y = [ buf.start, buf.start, buf.end, buf.end ]
            i = table_row.index[ table_row['scaffold'] == buf.seqid ][0]
            for j in range( 0, len( table_column ) ):
                x_max = table_column.at[j, 'len']
                x = [ 0, x_max, x_max, 0 ]
                ax[i][j].fill( x, y, color=buf.color, linewidth=0, alpha=alpha )
        #x軸scaffoldへのhighlight
        if buf.seqid in table_column['scaffold'].values:
            x = [ buf.start, buf.end, buf.end, buf.start ]
            j = table_column.index[ table_column['scaffold'] == buf.seqid ][0]
            for i in range( 0, len( table_row)):
                y_max = table_row.at[i, 'len']
                y = [0, 0, y_max, y_max ]
                #print( x, y )
                ax[i][j].fill( x, y, color=buf.color, linewidth=0, alpha=alpha )


def func_plot_crossed( highlight, table_row, table_column, ax, alpha ):
    alpha = alpha
    for buf in highlight.itertuples():
        if not clr.is_color_like( buf.color ):
            print( f"Error: {buf.color} is invalid color. Calculation will continue.", file=sys.stderr )
            continue
        if buf.seqid in table_row['scaffold'].values and buf.seqid in table_column['scaffold'].values :
            y = [ buf.start, buf.start, buf.end, buf.end ]
            i = table_row.index[ table_row['scaffold'] == buf.seqid ][0]
            x = [ buf.start, buf.end, buf.end, buf.start ]
            j = table_column.index[ table_column['scaffold'] == buf.seqid ][0]
            ax[i][j].fill( x, y, color=buf.color, linewidth=0, alpha=alpha )

        for buf2 in highlight.itertuples():
            if buf.Index >= buf2.Index:
                continue
            if not clr.is_color_like( buf2.color ):
                print( f"Error: {buf2.color} is invalid color. Calculation will continue.", file=sys.stderr )
                continue
            if buf.color != buf2.color :
                continue
            #y軸scaffoldへのhighlight
            if buf.seqid in table_row['scaffold'].values and buf2.seqid in table_column['scaffold'].values :
                y = [ buf.start, buf.start, buf.end, buf.end ]
                i = table_row.index[ table_row['scaffold'] == buf.seqid ][0]
                x = [ buf2.start, buf2.end, buf2.end, buf2.start ]
                j = table_column.index[ table_column['scaffold'] == buf2.seqid ][0]
                ax[i][j].fill( x, y, color=buf.color, linewidth=0, alpha=alpha )

            #x軸scaffoldへのhighlight
            if buf.seqid in table_column['scaffold'].values and buf2.seqid in table_row['scaffold'].values :
                x = [ buf.start, buf.end, buf.end, buf.start ]
                j = table_column.index[ table_column['scaffold'] == buf.seqid ][0]
                y = [ buf2.start, buf2.start, buf2.end, buf2.end ]
                i = table_row.index[ table_row['scaffold'] == buf2.seqid ][0]
                ax[i][j].fill( x, y, color=buf.color, linewidth=0, alpha=alpha )


def func_print_messages():
    print( 'start' )
    print( ' '.join( sys.argv ))


if __name__ == '__main__':
    main()
