import os,sys
from glob import glob
import pandas
import numpy as np
import nibabel as ni
import statsmodels.formula.api as smf
sys.path.insert(0,'/Users/jakevogel/git/tau_atlas/')
from image_transforms import slices_to_skip, compile_slice, reorient_mtx

def main(scans,vardict,outnm,outpth=None,mask=None,resids=False):

    if not outpth:
        outpth = os.getcwd()

    print 'creating dataframe'
    df = create_data(vardict)
    if 'x' not in df.columns.tolist():
        raise IOError('vardict must contain a key called x')

    print 'retriving information'
    inf = ni.load(scans[0])
    x,y,z = inf.get_shape()
    affine = inf.get_affine()

    if mask:
        print 'masking'
        skp,msk_d = slices_to_skip(mask)
        tst_block = compile_slice(scans,0)
        bx,by,bz = tst_block.shape

    for i in range(x):
        print 'working on block %s'%(i)
        if mask:
            if i in skp:
                print 'skipping %s'%(i)
                emt_slc = np.full((bx,by,1),np.nan)
                if i == 0:
                    res_t = emt_slc
                    res_p = emt_slc
#                    if resids:
#                        blaaaa
                else:
                    tcont = [res_t]
                    pcont = [res_p]
                    tcont.append(emt_slc)
                    pcont.append(emt_slc)
                    res_t = np.concatenate(tcont,axis=2)
                    res_p = np.concatenate(pcont,axis=2)
            else:
                block = compile_slice(scans,i)
                if i == 0:
                    if resids:
                        res_t,res_p,resid = run_stats(block,df,resids)
                    else:
                        res_t,res_p = run_stats(block,df)
                else:
                    tcont = [res_t]
                    pcont = [res_p]
                    if resids:
                        t,p,resid = run_stats(block,df,resids)
                    else:
                        t,p = run_stats(block,df)
                        tcont.append(t)
                        pcont.append(p)
                        res_t = np.concatenate(tcont,axis=2)
                        res_p = np.concatenate(pcont,axis=2)

        else:
            block = compile_slice(scans,i)
            if i == 0:
                if resids:
                    res_t,res_p,resid = run_stats(block,df,resids)
                else:
                    res_t,res_p = run_stats(block,df)
            else:
                tcont = [res_t]
                pcont = [res_p]
                if resids:
                    t,p,resid = run_stats(block,df,resids)
                else:
                    t,p = run_stats(block,df)
                    tcont.append(t)
                    pcont.append(p)
                    res_t = np.concatenate(tcont,axis=2)
                    res_p = np.concatenate(pcont,axis=2)

    print 'reorienting images...'
    res_t = reorient_mtx(res_t)
    res_p = reorient_mtx(res_p)

    print 'writing output'
    n_timg = ni.Nifti1Image(res_t, affine)
    n_timg.to_filename(os.path.join(outpth,'%s_t'%(outnm)))
    n_pimg = ni.Nifti1Image(res_p, affine)
    n_pimg.to_filename(os.path.join(outpth,'%s_p'%(outnm)))

    print 'results written to %s'%(outpth)

def create_data(vardict):
    df = pandas.DataFrame(np.full((len(vardict['x']),len(vardict)),np.nan),columns=vardict.keys())
    for k,v in vardict.iteritems():
        for i in range(len(v)):
            df.ix[i,k] = v[i]

    return df

def run_stats(block, df, resids=False):

    x,y,z = block.shape
    res_t = np.full((x,y),np.nan)
    res_p = np.full((x,y),np.nan)
    if resids:
        resid = np.full((x,y,z),np.nan)
    for i in range(x):
        for j in range(y):
#            print i,j
            ny = block[i][j][:]
            check = check_relevance(ny)
            if not check:
                continue
            else:
                for ind in df.index.tolist():
                    df.ix[ind,'y'] = ny[ind]
                stmnt = build_statement(df)
                lm = smf.ols(stmnt,data=df).fit()
                res_t[i][j] = lm.tvalues[1]
                res_p[i][j] = lm.pvalues[1]
                if resids:
                    for k in range(z):
                        resid[i][j][k] = lm.resid[k]

    res_t = res_t.reshape(x,y,1)
    res_p = res_p.reshape(x,y,1)

    if resids:
        return res_t,res_p,resid
    else:
        return res_t,res_p

def build_statement(df):
    stmnt = 'y ~ x'
    varlist = df.columns.tolist()
    for col in varlist:
        if col != 'x':
            stmnt = stmnt+'+ %s'%(col)

    return stmnt

def check_relevance(y):
    check = []
    for x in y:
        if not abs(x) < 0.000000000001:
            check.append(x)

    return check
#def handle_resids(sub_n,iter_n,slc,aff):

#    x.y,z = slc.shape
#    if iter_n == 0:
#        for i in range(z):
#            

