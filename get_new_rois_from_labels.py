import os
import nibabel as ni
import numpy as np
from copy import deepcopy
import pandas

def main(in_atlas,new_rois,outfl,staging='Off',sizes='True'):
    '''Given a dictionary (new_rois), where key=index(new label) and value = a list of
    existing ROI labels, will output a new atlas or new atlases corresponding
    to the labels given in new_rois. Use staging if you want cumulative ROIs
    instead or in addition.

    in_atlas = a nifti image whose voxels are integers corresponding to ROI
    labels

    new_rois = a dictionary where key = a new label and value = a list of
    existing ROI labels. All ROIs in value will be merged under the label of
    key. Keys should be consecutive integers starting from 1, particularly if
    staging is not set to 'Off'.

    outfl = a general filename structure. Should include path to output
    directory. Should not include any extension.

    staging = if not set to 'Off', will create cumulative ROIs, such that each
    consecutive ROI label will also include all ROIs from the previous label.
    Separate images will therefore be created for each key in new_rois.
        'Off' -> only generate re-atlased image based on new_rois, do not
        create stage ROIs
        'only' -> only generate stage ROI images. Do not create re-atlased
        image
        'also' -> generate both stage ROI images and re-atlased image.

    sizes = If True, will generate spreadsheets containing information about
    the size of each ROI (i.e. number of voxels). Will work no matter what
    staging is set to.
    '''

    valid = validate(new_rois,staging)
    if not valid:
        raise ValueError('in_atlas must be a dictionary where each value is a list (and each key is an integer)')

    jnk = ni.load(in_atlas)
    atl = jnk.get_data()
    aff = jnk.affine

    if staging != 'only':
        natl = create_basic_atlas(atl,new_rois)
        if sizes:
            sdf = get_roi_sizes(natl,new_rois)
            sdf.to_excel('%s_new_atlas_sizes.xls'%(outfl))
        print('writing image')
        nimg = ni.Nifti1Image(natl, aff)
        nimg.to_filename('%s_new_atlas'%(outfl))

    if staging != 'Off':
        satlz = create_stage_rois(atl,new_rois)
        for i,satl in satlz.items():
            if sizes:
                sdf = get_roi_sizes(satl,{i:[]})
                sdf.to_excel('%s_stage_%s_size.xls'%(outfl,i))
            print('writing stage image %s'%(i))
            nimg = ni.Nifti1Image(satl,aff)
            nimg.to_filename('%s_stage_%s'%(outfl,i))

def validate(roi_labs,staging):

    print('validating')
    valid = True
    if staging != 'Off':
        if any(type(x) != int for x in roi_labs.keys()):
            valid = False
    for k,v in roi_labs.items():
        if type(v) != list:
            valid = False

    return valid

def create_basic_atlas(atl,new_rois):

    print('reatlasing...')
    natl,uni = clean_atlas(atl,new_rois)
    lab1 = sorted(list(new_rois.keys()))[0]
    imgz = {}
    for k,v in new_rois.items():
        iatl = deepcopy(natl)
        if k not in v:
            iatl[iatl==k] = 0 # a consequence of relabeling
        for lab in v:
            iatl[iatl==lab] = k
        iatl[iatl!=k] = 0
        imgz.update({k:iatl})

    natl = imgz[lab1]
    for k,v in imgz.items():
        if k != lab1:
            natl = natl + v

    return natl

def create_stage_rois(atl,new_rois):

    print('creating stage ROIs')
    natl,uni = clean_atlas(atl,new_rois)
    n_dict = {}
    for i in range(1,len(new_rois.keys())+1):
        nstage = []
        for j in range(1,i+1):
            nstage = nstage+new_rois[j]
        n_dict.update({i: nstage})

    imgz = {}
    for k,v in n_dict.items():
        print('working on stage %s'%(k))
        iatl = deepcopy(natl)
        if k not in v:
            iatl[iatl==k] = 0 # a consequence of relabeling
        for lab in v:
            iatl[iatl==lab] = k
        iatl[iatl!=k] = 0
        imgz.update({k: iatl})

    return imgz

def clean_atlas(atl,new_rois):

    print('cleaning atlas')
    natl = deepcopy(atl)
    newlabz = []
    for k,v in new_rois.items():
        newlabz = newlabz+v
    uni = set(natl.flat)
    n_labz = list(uni - set(newlabz))
    for lab in n_labz:
        natl[natl==lab] = 0

    for k,v in new_rois.items():
        for lab in v:
            if lab not in uni:
                print('**WARNING*** passed label %s not included in in_atlas'%(lab))

    return natl,uni

def get_roi_sizes(atl,rois=None):

    print('generating ROI size information...')
    if rois == None:
        unique = list(set(atl.flat))
    else:
        unique = rois.keys()
    sdf = pandas.DataFrame(index=unique,columns=['size'])
    for k in unique:
        if k == 0 or k == 0.0:
            continue
        else:
            sdf.ix[k,'size'] = len(atl[atl==k].flat)

    return sdf
