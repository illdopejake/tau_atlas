import numpy as np
import nibabel as ni
import pandas
import itertools

def main(in_img, outnm, thr, thr_tp, k_thr=None,out_tp='df',stat2_img = None):
    '''given a t-map or p-map and an intensity threshold (and cluster threshold
    if desired), will return a df of independent clusters surpassing that
    threshold. Additionally, if desired, will output those clusters to an nifti
    image.

    in_img = an t-map or p-map image, nifti format. Image should be masked to
    make things quicker.

    outnm = the string indicating the name associated with output files. Please
    no extensions, but do include path or else will output to cwd.

    thr = The intensity threshold. Only voxes above this theshold will be
    considered in the cluster extraction.

    thr_tp = accepts 'tpos', 'tneg' or 'p'. Tpos and tneg indicate input image
    is a tmap, and pos/neg indicates the direction of the threshold (i.e. a t
    thr of -4.65 would be tneg, and only tvalues below that value would be
    assessed). Set to p for a p-map image.

    k_thr = Optional. Input a threshold for cluster size. Clusters below this
    cluster size will not be included in the output. Inputting even a low k_thr
    may speed up computation.

    out_tp = Choose from df, img_b, img_c. df will only return a dataframe with
    peak stats and coordinates. img_b will return an image of all voxels in
    clusters surviving given thresholds. img_c will do the same, but will have
    different values for each cluster so that they can be differentiated.

    stat2_img = Optional. A nifti image in the same space as in_img. Output
    dataframe will include values from this image at each peak. For example, if
    in_img is a tmap, you can input the associated pmap for stat2_img, and your
    output df will include both t and p values for each peak.
    '''

    # make sure inputs are set correctly
    if thr_tp != 'tpos' and thr_tp != 'tneg' and thr_tp != 'p':
        raise NameError('thr_tp must be set as either tpos, tneg or p')

    if thr_tp == 'tneg' and thr>0:
        thr = thr * -1

    if thr_tp == 'tpos' and thr<0:
        thr = thr * -1

    if thr_tp == 'p' and thr<0 or thr>1:
        raise ValueError('if thr_tp is set to p, thr must be a value between 0 and 1')

    if out_tp != 'df' and out_tp != 'img_b' and out_tp !='img_c':
        raise NameError('out_tp must be seat as either df, img_b or img_c')

    # load data
    img = ni.load(in_img)
    dat = img.get_data()
    aff = img.get_affine()

    if stat2_img:
        stat2 = ni.load(stat2_img).get_data()
    else:
        stat2 = None

    # create dataframe
    df = pandas.DataFrame(columns=['x','y','z','stat','stat2','k'])

    # intialize containers
    peaks = []
    count = 1000
    keep_searching = True

    #threshold data
    dat = threshold_data(dat,thr,thr_tp)

    while keep_searching == True:
        # find peak
        peak,coord,keep_searching = find_peak(dat,thr,thr_tp,k_thr)
        peaks.append(peak)
        # find associated cluster
        if coord:
            clust = find_cluster(dat,coord,thr,thr_tp)
            k = len(clust)
        # update dataframe
            if k_thr:
                if k > k_thr:
                    df = update_df(df,coord,peak,count,len(clust),stat2)
            else:
                df = update_df(df,coord,peak,count,len(clust),stat2)
        # mark voxels in cluster
            dat = mark_cluster(dat,clust,count)
            count = count+1

    print 'writing output, outfile prefix is %s'%(outnm)
    if out_tp != 'df':
        dat[dat<1000] = 0
        if out_tp == 'img_b':
            dat[dat>1000] = 1

        nimg = ni.Nifti1Image(dat, aff)
        nimg.to_filename('%s_clusters'%(outnm))

    df.to_excel('%s.xls'%(outnm))

## Use a while loop to recursively find peaks, explore clusters around peak,
## and save coordinates of within-cluster voxels. Then, label those clusters
## and repeat

def threshold_data(dat,thr,thr_tp):

    # denan
    nanlox = np.isnan(dat)
    dat[nanlox] = 0

    if thr_tp == 'tpos':
        dat[dat<thr] = 0
    else:
        dat[dat>thr] = 0

    return dat

def find_peak(ndat,thr,thr_tp,k_thr = None):

    keep_searching = True

    sdat = ndat[ndat<1000]
    if len(sdat) < 1:
        keep_searching = False
    else:
        print 'identifying peak'
        if thr_tp == 'tpos':
            dist = sorted(sdat[sdat>thr])
        else:
            dist = sorted(sdat[sdat<thr])
        if k_thr:
            if len(dist) < k_thr:
                keep_searching = False
        if len(dist) < 1:
            keep_searching = False
        else:
            if thr_tp == 'tpos':
                peak = dist[-1]
            else:
                peak = dist[0]

    if keep_searching == True:
        idx = np.where(ndat == peak)
        coord = (idx[0][0],idx[1][0],idx[2][0])
    else:
        peak = None
        coord = None

    return peak,coord,keep_searching

def update_df(df,coord,peak,count,k,stat2=None):
    nm = 'peak_%s'%(count-1000)
    df.ix[nm,'x'] = coord[0]
    df.ix[nm,'y'] = coord[1]
    df.ix[nm,'z'] = coord[2]
    df.ix[nm,'stat'] = peak
    if stat2 == np.ndarray:
        df.ix[nm,'stat2'] = stat2[x][y][z]
    df.ix[nm,'k'] = k

    return df

def find_cluster(dat,coord,thr,thr_tp):

    todo = [coord]
    clust = []

    while len(set(todo) - set(clust)) > 0:
        print 'new iteration'
        for vox in todo:
            neighbors,neighbors2 = get_neighbors(vox)
            clean = assess_neighborhood(dat,vox,clust,neighbors,neighbors2,thr,thr_tp)
            for cl in clean:
                if cl not in clust and cl not in todo:
                    todo.append(cl)
            clust.append(vox)
            #todo.remove(vox)
            print len(clust),len(todo)

    return clust

def get_neighbors(coord):
    neighbors = []
    nx = list(coord)
    for i in itertools.product(*[[0,1,-1]]*len(nx)):
        jnk = tuple(np.add(list(i),nx))
        neighbors.append(jnk)

    neighbors2 = []
    for i in itertools.product(*[[0,2,-2]]*len(nx)):
        jnk = tuple(np.add(list(i),nx))
        if jnk not in neighbors:
            neighbors2.append(jnk)

    return neighbors, neighbors2

def assess_neighborhood(dat,coord,clust,neighbors,neighbors2,thr,thr_tp):

    clean = []

    # make sure neighbors are real image coordinates
    neighbors = validate_neighbors(dat,neighbors)
    #neighbors,neighbors2 = validate_neighbors(dat,neighbors)

    for i,n in enumerate(neighbors):
        cval = dat[n[0]][n[1]][n[2]]
        if thr_tp == 'tpos':
            if cval > thr and cval < 1000 and n not in clust:
                clean.append(n)
#            else:
#                n2 = neighbors2[i]
#                cval2 = dat[n2[0]][n2[1]][n2[2]]
#                if cval2 > thr and cval2 < 1000 and n2 not in clust:
#                    clean.append(n)
        else:
            if cval < thr and cval < 1000 and n not in clust:
                clean.append(n)
#            else:
#                n2 = neighbors2[i]
#                cval2 = dat[n2[0]][n2[1]][n2[2]]
#                if cval2 < thr and cval2 < 1000 and n2 not in clust:
#                    clean.append(n)

    return clean

def validate_neighbors(dat,neighbors):
    invalid = []
    for n in neighbors:
        if n[0] < 0 or n[1] < 0 or n[2] < 0:
            invalid.append(n)
        else:
            try:
                dat[n[0]][n[1]][n[2]]
            except:
                invalid.append(n)

    # IF WE END UP INCLUDING N2
#    ndict = {}
#    for i,n in enumerate(neighbors):
#        ndict.update({n: i}) 

    for inv in invalid:
        neighbors.remove(inv)

    # IF WE END UP INCLUDING N2
#    nn2 = []
#    for n in neighbors:
#        nind = ndict[n]
#        nn2.append(neighbors2[nind])

#    return neighbors,nn2
    return neighbors

def mark_cluster(dat,clust,count):

    for vox in clust:
        dat[vox[0]][vox[1]][vox[2]] = count

    return dat
