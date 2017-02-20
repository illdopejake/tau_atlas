import os
from glob import glob
import pandas
#import random, math
import numpy as np
import nibabel as ni
import math

def z_transform(scans,outpth=None,int_flz =[],keep_ints=False,mask=None):
    ''' Given a list of images, will z-transform them

        Scans is a list of paths to images to be transformed

        outpth is a path where the output should be written to. If None,
        current working directory is used

        int_flz is a list containing mean image and SD image, in that order. If
        None, voxelwise means and SDs will be computed.

        keep_ints -- if True, mean and SD images will be written

        mask -- working on it...
        '''

    if not outpth:
        outpth = os.getcwd()

    if int_flz:
        keep_ints=False

#	try:
#        skp,msk_d = slices_to_skip(mask)
#    except:
#    	pass
    inf = ni.load(scans[0])
    x,y,z = inf.get_shape()
    affine = inf.get_affine()

    if int_flz:
        mnz = ni.load(int_flz[0]).get_data()
        sdz = ni.load(int_flz[1]).get_data()

    else:
        for i in range(z):
#    	    if mask:
#    		    if i in skp:
#    			    continue
            print 'preparing slice %s of %s'%(i,z)
            block = compile_slice(scans,i)
            if i == 0:
                mnz,sdz = get_means_and_sdz(block)
            else:
                mcont = [mnz]
                scont = [sdz]
                mn,sd = get_means_and_sdz(block)
                mcont.append(mn)
                scont.append(sd)
                mnz = np.concatenate(mcont,axis=2)
                sdz = np.concatenate(scont,axis=2)

        print 'reorienting images...'
        mnz = reorient_mtx(mnz)
        sdz = reorient_mtx(sdz)

    if keep_ints:
        i_mnz = ni.Nifti1Image(mnz,affine)
        i_sdz = ni.Nifti1Image(sdz,affine)
        i_mnz.to_filename(os.path.join(outpth,'means'))
        i_sdz.to_filename(os.path.join(outpth,'SDs'))

    for scan in scans:
        aff,dat = compute_zscore(scan,mnz,sdz)
        if len(os.path.split(scan)) > 2:
            pth,fl = os.path.split(scan)[-1]
            newfl = os.path.join(outpth,'zscored_%s'%(fl))
        else:
            newfl = os.path.join(outpth,'zscored_%s'%(scan))
        n_img = ni.Nifti1Image(dat, aff)
        n_img.to_filename(newfl)

def compile_slice(scans, slc_no):

    jnk = ni.load(scans[0])
    block = jnk.get_data()[slc_no][:][:]
    x,y,z = jnk.get_shape()
    block = block.reshape(y,x,1)
    for scan in scans[1:]:
        cont = [block]
        slc = ni.load(scan).get_data()[slc_no][:][:]
        slc = slc.reshape(y,x,1)
        cont.append(slc)
        block = np.concatenate(cont,axis=2)

    return block


def slices_to_skip(mask):

    skp = []
    jnk = ni.load(mask)
    msk_d = jnk.get_data()
    x,y,z = jnk.get_shape()
    for itr in range(x):
        print '%s'%(itr)
        nonz=[]
        lst = msk_d[itr][:][:].flatten()
        for i in lst:
            if i != 0:
                nonz.append(i)
        if len(nonz) == 0:
            skp.append(itr)

    return skp,msk_d


def get_means_and_sdz(block):
    x,y,z = block.shape

    mtx_rsl = block.reshape(x*y,z)
    mnz = np.array([i/float(z) for i in map(np.sum, mtx_rsl)]).reshape(x,y,1)
    var = [i/float(z) for i in map(np.sum,np.subtract(mtx_rsl,np.array(mnz).reshape(x*y,1))**2)]
    stdz = np.array([math.sqrt(k) for k in var]).reshape(x,y,1)

    # the above code should double the speed of the code below, simply by using
    # map and pure python instead of np. Strange.

#    mnz = np.full((x,y),np.nan)
#    stdz = np.full((x,y),np.nan)

#    for i in range(x):
#        for j in range(y):
#            if block[i][j][0] != 0:
#                mnz[i][j] = np.mean(block[i][j][:])
#                stdz[i][j] = np.std(block[i][j][:])

#    mnz = mnz.reshape(x,y,1)
#    stdz = stdz.reshape(x,y,1)

    return mnz,stdz

def compute_zscore(scan,mnz,sdz):
    print 'z-scoring subject %s'%(scan)
    jnk = ni.load(scan)
    x,y,z = jnk.get_shape()
    aff = jnk.get_affine()
    dat = jnk.get_data()
    for i,j,k in np.nditer([dat,mnz,sdz],flags=['external_loop'],op_flags=['readwrite']):
        i[...] = (i-j)/k

### THE ABOVE LOOP SPEEDS UP THE LOOP BELOW CODE MORE THAN TENFOLD ### 
#    for i in range(x):
#        for j in range(y):
#            for k in range(z):
#                val = dat[i][j][k]
#                dat[i][j][k] = (val-mnz[i][j][k])/sdz[i][j][k]

    return aff,dat

def reorient_mtx(in_mtx):
    x,y,z = in_mtx.shape
    o_mtx = np.full((z,x,y),np.nan)

    for i in range(x):
        for j in range(y):
            for k in range(z):
                o_mtx[k][i][j] = in_mtx[i][j][k]

    return o_mtx
