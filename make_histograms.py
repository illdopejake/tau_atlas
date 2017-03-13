from matplotlib import pyplot as plt
from matplotlib import mlab
import pandas

def main(sheet, out, cols=[],skp=False,style='single',lof=False,frt = 'pdf', xlab= 'Variable X',title='I am a histogram'):

# right now, LOF isn't working...

    if style != 'single' and style != 'overlay' and style != 'interleave':
        raise NameError('style must be set to single, overlay or interleave.Name %s not recognized'%(style))

    if type(sheet) == str:
        if 'xls' in sheet[-4:]:
            try:
                df = pandas.ExcelFile(sheet).parse('Sheet1')
            except:
                df = pandas.read_csv(sheet)
        else:
            df = pandas.read_csv(sheet)
    else:
        df = sheet

    if len(cols) > 0:
        ncols = []
        for col in cols:
            if type(col) == str:
                if col in df.columns:
                    ncols.append(col)
                else:
                    print('%s not found in spreadsheet, skipping'%(col))
            else:
                try:
                    ncols.append(df.columns[col])
                except:
                    raise ValueError('cols must be empty, or must contain only string, integers and slices')
    else:
        ncols = df.columns.tolist()

    if skp:
        ncols.remove(ncols[0])

    # make histogram
    if style == 'overlay' or style == 'interleave':
        hist = plt.figure()
        plt.title(title)
        plt.xlabel(xlab)
        plt.ylabel("Frequency")

    if style == 'interleave':
        histyz = []
        if lof:
            yz = []

    for col in ncols:
        distr = df[:][col].tolist()
        if style == 'single':
            hist = plt.figure()
            n,bins,patches = plt.hist(distr)
            plt.title(title)
            plt.xlabel(xlab)
            plt.ylabel("Frequency")
            if lof:
                y = line_of_fit(n,bins,patches)
                plt.plot(bins, y, '--')
            hist.savefig('%s_%s.pdf'%(out,col),bbox_inches = 'tight')
            print('created histogram for column %s'%(col))
        elif style == 'overlay':
            n,bins,patches = plt.hist(distr)
            if lof:
                y = line_of_fit(n,bins,patches)
                plt.plot(bins, y, '--')
        else:
            histyz.append(distr)

    if style == 'overlay':
        hist.savefig('%s_hist_overlay.pdf'%(out),bbox_inches = 'tight')
        print('histogram written')
    if style == 'interleave':
        plt.hist(histyz)
        hist.savefig('%s_hist_interleave.pdf'%(out),bbox_inches = 'tight')
        print('histogram written')

def line_of_fit(n,bins,patches):
    y = mlab.normpdf(bins, mu, sigma)
    return y
