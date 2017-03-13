import os,sys
from glob import glob
import pandas
import subprocess
import nibabel as ni
import numpy as np
sys.path.insert(0,'/Users/jakevogel/git/tPSO_scripts/')
#from propagation_correlations import codegen

def codegen(N):
    cde = ''.join(random.choice(string.ascii_lowercase) for _ in range(N)) 
    return cde

def extract_roi_values(img_dict, atlas, otpt, img_names, roi_key = None,scale=None, comp = 'mean', method='nibabel',weight = False):
    '''
    extracts values from several images given an atlas

    img_dict is a dict where the key is the subject ID and the value is a
    list of images to extract values from

    atlas is a nifti image with several ROIs such that each voxel corresponds
    with an integer label that references the label of that ROI

    otpt is the file name of the output spreadsheet. DONT PUT AN EXTENSION

    img_names -- since the keys in img_dict are lists, the script needs labels
    that correspond to each "type" of image in the list. img_names is therefore
    a dict where the value is a search string that is unique to a type of image
    in the img_dict value list, and the key is a string label that corresponds
    to that image type

    NOTE: You do not need to pass roi_key OR scale. It looks like the current
    version requires you to pass one or the other, but it also looks like I
    wrote in something that will work even if none are passed... so maybe pass
    one or the other for now and I'll get back to this later.

    roi_key is a dictionary where key is an integer reflecting the ROI label
    corresponding to ROI in atlas, and value is a string label for the atlas.

    scale is an integer that corresponds to the number of distinct ROIs in the
    atlas.

    method can be set to either nibabel or fsl. fsl uses command line fsl
    commands. nibabel actually opens the image and does the extraction
    manually. nibabel is like 10-20x faster, probably more. So, probably leave
    that as the method...

    comp can be set to either mean or count. mean returns the mean of voxels
    within each ROI. count returns the number of non-zero voxels in the roi
    (for the case of binarized images)

    weight can be set to True if you wish the extracted values to be weighted
    by ROI size. This will only work if method is set to nibabel
    '''

    if method != 'nibabel' and method != 'fsl':
        raise IOError('please set method to either fsl or nibabel')


    if method != 'nibabel' and weight == True:
        raise IOError('if weight is True, method must be set to nibabel')
    print('using method %s'%(method))

    wdir = os.getcwd()
    if roi_key == None and scale == None:
        raise IOError('scale or roi_key  must be passed')

    df = pandas.DataFrame(index = img_dict.keys())

    for sid,imgz in img_dict.items():
        for img in imgz:
            for k,v in img_names.items():
                if v in img:
                    itp = k
            if method == 'fsl':
                df = extract_values_from_img(img,atlas,df,sid,itp,roi_key=roi_key,scale=scale,comp=comp)
            elif method == 'nibabel':
                df = extract_values_nibabel(img,atlas,df,sid,itp,roi_key=roi_key,scale=scale,comp=comp)
            df.to_csv(os.path.join(wdir,'%s.csv'%(otpt)))
            print('file written to %s'%(os.path.join(wdir,'%s.csv'%(otpt))))

    return df

def extract_values_from_img(img,atlas,df,sid,itp,roi_key = None, scale = None,comp = 'mean'):

    print('extracting values from %s'%(img))
    if not roi_key:
        for i in range(1,(scale+1)):
            print('working on roi %s'%(i))
            cde = codegen(6)
            os.system('fslmaths %s -thr %s -uthr %s %s_msk'%(atlas,i,i,cde))
            os.system('fslmaths %s -mas %s_msk.nii.gz %s_valz'%(img,cde,cde))
            if comp == 'count':
                val = subprocess.check_output('fslstats %s_valz.nii.gz -V'%(cde),shell = True)
            else:
                val = subprocess.check_output('fslstats %s_valz.nii.gz -M'%(cde),shell = True)
            os.system('rm %s_*'%(cde))
            df.ix[sid, '%s_%s'%(itp,i)] = float(val)

    else:
        for i,roi in roi_key.items():
            print('working on roi %s'%(roi))
            cde = codegen(6)
            os.system('fslmaths %s -thr %s -uthr %s %s_msk'%(atlas,i,i,cde))
            os.system('fslmaths %s -mas %s_msk.nii.gz %s_valz'%(img,cde,cde))
            val = subprocess.check_output('fslstats %s_valz.nii.gz -M'%(cde),shell=True)
            os.system('rm %s_*'%(cde))
            df.ix[sid, '%s_%s'%(itp,roi)] = float(val)

    return df

def extract_values_nibabel(img,atlas,df,sid,itp,roi_key = None, scale = None, comp = 'mean', weight = False):

    adata = ni.load(atlas).get_data()
    idata = ni.load(img).get_data()
    if not roi_key and not scale:
        print('determining atlas characteristics')
        a,b,c = adata.shape
        unique = list(set(adata.flat))
#        unique = []
#        for i in range(a):
#            for j in range(b):
#                for k in range(c):
#                    if adata[i][j][k] not in unique:
#                        unique.append(adata[i][j][k])
        print('%s unique atlas elements found in %s'%(len(unique), atlas))
        for i in unique:
            print('working on roi %s'%(i))
            msk = np.ma.masked_array(adata,adata==i)
            if weight:
                sz = len(adata[adata==i])
                df.ix[sid,'%s_%s'%(itp,i)] = (np.mean(idata[msk.mask]) * sz) / (np.mean(idata[msk.mask]) + sz)
            else:
                df.ix[sid, '%s_%s'%(itp,i)] = np.mean(idata[msk.mask])

    if not roi_key and scale:
        print('extracting values from %s using scale %s'%(img,scale))
        for i in range(1,(scale+1)):
            print('working on roi %s'%(i))
            msk = np.ma.masked_array(adata,adata==i)
            if comp == 'count':
                df.ix[sid, '%s_%s'%(itp,i)] = np.mean(idata[msk.mask]) * len(idata[msk.mask])
            else:
                if weight:
                    sz = len(adata[adata==i])
                    df.ix[sid,'%s_%s'%(itp,i)] = (np.mean(idata[msk.mask]) * sz) / (np.mean(idata[msk.mask]) + sz)
                else:
                    df.ix[sid, '%s_%s'%(itp,i)] = np.mean(idata[msk.mask])

    if roi_key:
        print('extracting values from %s using roi_key'%(img))
        for i,roi in roi_key.items():
            print('working on roi %s'%(roi))
            msk = np.ma.masked_array(adata,adata==i)
            df.ix[sid, '%s_%s'%(itp,roi)] = np.mean(idata[msk.mask])
            if comp == 'count':
                df.ix[sid, '%s_%s'%(itp,i)] = np.mean(idata[msk.mask]) * len(idata[msk.mask])
            else:
                if weight:
                    sz = len(adata[adata==i])
                    df.ix[sid,'%s_%s'%(itp,i)] = (np.mean(idata[msk.mask]) * sz) / (np.mean(idata[msk.mask]) + sz)
                else:
                    df.ix[sid, '%s_%s'%(itp,i)] = np.mean(idata[msk.mask])
    return df
