'''
statistical analysis functions, fitting, optimization and the like
'''

import numpy as np
from numpy import cos, sin, tan, sqrt, log, exp, sum
from copy import copy, deepcopy#HMCC

def peaks(x, y, n=0):
    '''

    '''
    maxs = {}

    if len((np.where(y == y.max()))[0]) == 1:
        for i in np.arange(1, len(x)-1):
            if (y[i] - y[i-1]) > 0 and (y[i+1] - y[i]) < 0:
                maxs[y[i]] = x[i]
    else:
        for i in np.arange(2, len(x)-1):
            if bool((y[i-1] - y[i-2]) > 0 and (y[i] - y[i-1]) == 0 and (y[i+1] - y[i]) < 0):
                maxs[y[i]] = x[i]

    vals  = sorted(maxs.keys())
    f1 = []
    f2 = []
    for v in reversed(vals):
        f1.append(maxs[v])
        f2.append(v)

    if n>0:
        return np.array(f1[0:n]), np.array(f2[0:n])
    else:
        return np.array(f1), np.array(f2)

def gs_search(f, bracket, tol=1.e-6, nmax=50):
    '''
    golden section search
    '''
    n = 0
    a,b,c = bracket
    cur_tol = c-a

    w = 0.38197

    while True:
        if n > nmax: break
        if cur_tol < tol: break

        if c-b > b-a:
            x = b + (c-b)*w
            if f(b)<f(x):
                c = x
            else:
                a = b
                b = x
        else:
            x = b - (b-a)*w
            if f(b)<f(x):
                a = x
            else:
                c = b
                b = x

        n += 1
        cur_tol = c-a

    return np.array([a,b,c]), cur_tol

def fit_gauss_2d(x,y,F):

    n1 = len(x)
    n2 = len(y)

    mu1 = 0
    mu2 = 0

    for i in range( n1 ): mu1 += x[i]*sum(F[i,:])
    for i in range( n2 ): mu2 += y[i]*sum(F[:,i])
    mu2 /= sum(F[:,:])
    mu1 /= sum(F[:,:])

    sig1 = 0

    for i in range( n1 ): sig1 += (x[i]-mu1)**2*sum(F[i,:])
    sig1 /= sum(F[:,:])

    sig2 = 0

    for i in range( n2 ): sig2 += (y[i]-mu2)**2*sum(F[:,i])
    sig2 /= sum(F[:,:])

    sig12 = 0

    for i in range( n1 ):
        for j in range( n2 ):
            sig12 += (x[i]-mu1)*(y[j]-mu2)*F[i,j]

    sig12 /= sum(F[:,:])

    rho = sig12 / sqrt(sig1*sig2)

    return mu1, mu2, sqrt(sig1), sqrt(sig2), rho

def fit_gauss_1d(x,F):

    n1 = len(x)

    mu1 = 0

    for i in range( n1 ): mu1 += x[i]*F[i]
    mu1 /= sum(F[:])

    sig1 = 0

    for i in range( n1 ): sig1 += (x[i]-mu1)**2*F[i]
    sig1 /= sum(F[:])

    return mu1, sqrt(sig1)

def fwhm(x,F):
    m = np.max(F) / 2.0
    ups = []
    downs = []
    for i in range(len(x)-1):
        if F[i] <  m and F[i+1] > m:
            ups.append(i)
        if F[i] >=  m and F[i+1] < m:
            downs.append(i)

    #print ups, downs
    return x[downs[-1]] - x[ups[0]]

def stats(outputs):

    ''' return mean, std, median and extreme (farthest from mean) of a time series '''

    omean = np.zeros_like(outputs[0])

    for o in outputs:
        omean += np.array(o)
    omean = omean / len(outputs)

    difm = np.linalg.norm(omean, 2)
    difw = 0.0
    imed = 0
    iworst = 0
    std = np.zeros_like(outputs[0])

    for i in np.arange( len(outputs) ):

        difn = np.linalg.norm(omean - np.array(outputs[i]), 2)

        if difn < difm:
            imed = i
            difm = difn

        if difn > difw:
            iworst = i
            difw = difn

        std += (omean - np.array(outputs[i]))**2

    std = np.sqrt(std / len(outputs))

    #print 'median id=', imed
    return omean , std, outputs[imed], outputs[iworst], imed, iworst


def find_saturation(power, z, z_max, n_smooth=5):
    p = np.diff(np.log10(power))

    u = np.convolve(p, np.ones(n_smooth) / float(n_smooth), mode='same')
    um = np.max(u)

    ii = 0

    for i in range(len(u)):
        if u[i] < 0.0 * um and z[i] > z_max:
            ii = i
            break

    #plt.plot(g.z[1:], u, lw=3)
    #plt.plot(g.z[ii+1], p[ii], 'rd')

    #plt.plot(g.z, power, lw=3)
    #plt.plot(z[ii+1], np.log10(power[ii]), 'rd')

    return z[ii+1], ii+1

def find_nearest_idx(array, value):
    if value == -np.inf:
        value = np.amin(array)
    if value == np.inf:
        value = np.amax(array)
    return (np.abs(array-value)).argmin()

def second_derivative_saturation(x,y,zmin,zmax): #HMCC  find_saturation
    from scipy.optimize import brentq
    from scipy.interpolate import UnivariateSpline

    z_min = [i_t for i_t,t in enumerate(x) if t>zmin][0]
    z_max = [i_t for i_t,t in enumerate(x) if t>zmax][0]

    x_i = deepcopy(x[z_min:z_max])
    y_i=deepcopy(y[z_min:z_max])


    y_spl = UnivariateSpline(x_i,y_i,k=5)
    y_spl_2d = y_spl.derivative(n=2)

    root2 = brentq(lambda xii: y_spl_2d(xii),zmin,zmax)
    return np.where(x>=root2)[0][0],root2,y[np.where(x>=root2)[0][0]]


def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def n_moment(x, counts, c, n):
    x = np.squeeze(x)
    if x.ndim is not 1:
        raise ValueError("scale of x should be 1-dimensional")
    if x.size not in counts.shape:
        raise ValueError("operands could not be broadcast together with shapes %s %s" %(str(x.shape), str(counts.shape)))
    
    if np.sum(counts)==0:
        return 0
    else:
        if x.ndim == 1 and counts.ndim == 1:
            return (np.sum((x-c)**n*counts) / np.sum(counts))**(1./n)
        else:
            
            if x.size in counts.shape:
                dim_ = [i for i, v in enumerate(counts.shape) if v == x.size]
                counts = np.moveaxis(counts, dim_, -1)
                return (np.sum((x-c)**n*counts, axis=-1) / np.sum(counts, axis=-1))**(1./n)
#def n_moment(x, counts, c, n):
#    if np.sum(counts)==0:
#        return 0
#    else:
#        return (np.sum((x-c)**n*counts) / np.sum(counts))**(1./n)

def std_moment(x, counts):
    mean=n_moment(x, counts, 0, 1)
    return n_moment(x, counts, mean, 2)

def bin_array(array,bin_size):
    #bining the array by averaging values within bins with size bin_size pixels
    if bin_size > len(array):
        return np.mean(array)
    elif bin_size == 1:
        return array
    else:
        new_shape = (array.shape[0] // bin_size) * bin_size
        array_av = array[:new_shape]
        array_av = array_av.reshape(int(new_shape/bin_size), bin_size)
        array_av = np.mean(array_av, axis=1)
        return array_av

def bin_scale(scale,bin_size):
    if bin_size > len(scale):
        return np.array([0])
    elif bin_size == 1:
        return scale
    else:
        hbin = np.int(bin_size/2) #halfbin (to pick bin centers)
        new_shape = (scale.shape[0] // bin_size) * bin_size
        new_scale = scale[hbin : new_shape+hbin] [ :: bin_size]
        return new_scale

def index_of(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

#HMCC migrating this function from gui/genesis_plot.py
def fwhm3(valuelist, height=0.5, peakpos=-1, total=1):
    """calculates the full width at half maximum (fwhm) of the array.
    the function will return the fwhm with sub-pixel interpolation. 
    It will start at the maximum position and 'walk' left and right until it approaches the half values.
    if total==1, it will start at the edges and 'walk' towards peak until it approaches the half values.
    INPUT:
    - valuelist: e.g. the list containing the temporal shape of a pulse
    OPTIONAL INPUT:
    -peakpos: position of the peak to examine (list index)
    the global maximum will be used if omitted.
    if total = 1 - 
    OUTPUT:
    - peakpos(index), interpolated_width(npoints), [index_l, index_r]
    """
    if peakpos == -1:  # no peakpos given -> take maximum
        peak = np.max(valuelist)
        peakpos = np.min(np.nonzero(valuelist == peak))

    peakvalue = valuelist[peakpos]
    phalf = peakvalue * height

    if total == 0:
        # go left and right, starting from peakpos
        ind1 = peakpos
        ind2 = peakpos
        while ind1 > 2 and valuelist[ind1] > phalf:
            ind1 = ind1 - 1
        while ind2 < len(valuelist) - 1 and valuelist[ind2] > phalf:
            ind2 = ind2 + 1
        grad1 = valuelist[ind1 + 1] - valuelist[ind1]
        grad2 = valuelist[ind2] - valuelist[ind2 - 1]
        if grad1 == 0 or grad2 == 0:
            width = None
        else:
            # calculate the linear interpolations
            # print(ind1,ind2)
            p1interp = ind1 + (phalf - valuelist[ind1]) / grad1
            p2interp = ind2 + (phalf - valuelist[ind2]) / grad2
            # calculate the width
            width = p2interp - p1interp
    else:
        # go to center from edges
        ind1 = 1
        ind2 = valuelist.size-2
        # print(peakvalue,phalf)
        # print(ind1,ind2,valuelist[ind1],valuelist[ind2])
        while ind1 < peakpos and valuelist[ind1] < phalf:
            ind1 = ind1 + 1
        while ind2 > peakpos and valuelist[ind2] < phalf:
            ind2 = ind2 - 1
        # print(ind1,ind2)
        # ind1 and 2 are now just above phalf
        grad1 = valuelist[ind1] - valuelist[ind1 - 1]
        grad2 = valuelist[ind2 + 1] - valuelist[ind2]
        if grad1 == 0 or grad2 == 0:
            width = None
        else:
            # calculate the linear interpolations
            p1interp = ind1 + (phalf - valuelist[ind1]) / grad1
            p2interp = ind2 + (phalf - valuelist[ind2]) / grad2
            # calculate the width
            width = p2interp - p1interp
        # print(p1interp, p2interp)
            
    return (peakpos, width, np.array([ind1, ind2]))

    ############ HMCC new version OCELOT ####
    
def corr_f_np(corr, val, n_skip=1, norm=1, count=0):
    n_val = corr.shape[0]
    for i in range(n_val):
        if count:
            sys.stdout.write('\r')
            sys.stdout.write('slice %i of %i' %(i, n_val-1))
            sys.stdout.flush()
        for j in range(n_val):
            means = np.mean(val[i*n_skip,:] * val[j*n_skip,:])
            meanl = np.mean(val[i*n_skip,:])
            meanr = np.mean(val[j*n_skip,:])
            if norm:
                corr[i,j] = means / meanl / meanr
            else:
                corr[i,j] = means - meanl * meanr
    
    if norm:
        corr[np.isnan(corr)] = 1
    else:
        corr[np.isnan(corr)] = 0

             
def correlation2d(val, norm=0, n_skip=1):
    N = int(val.shape[0] / n_skip)
    corr = np.zeros([N,N])
    corr_f_np(corr, val, n_skip, norm)
    return corr
    
def corr_c_np(corr, n_corr, val, norm):
    n_val = len(val) - n_corr*2
    for i in range(n_val):
        for j in range(n_corr):
            if not j%2:
                ind_l = int(i - j/2 + n_corr)
                ind_r = int(i + j/2 + n_corr)
            else:
                ind_l = int(i - (j-1)/2 + n_corr)
                ind_r = int(i + (j-1)/2 + 1 + n_corr)
            means = np.mean(val[ind_l,:] * val[ind_r,:])
            meanl = np.mean(val[ind_l,:])
            meanr = np.mean(val[ind_r,:])

            if meanl == 0 or meanr == 0:
                corr[i,j] = 0
            else:
                if norm:
                    corr[i,j] = means / meanl / meanr
                else:
                    corr[i,j] = means - meanl * meanr


def correlation2d_center(n_corr, val, norm=0):
    n_val, n_event = val.shape
    zeros = np.zeros((n_corr, n_event))
    val = np.r_[zeros, val, zeros]
    corr = np.zeros([n_val, n_corr])
    corr_c_np(corr, n_corr, val, norm)
    return corr
    
def mut_coh_func(J, fld, norm=1):
        """
        Mutual Coherence function
        """
        n_x = len(fld[0,0,:])
        n_y = len(fld[0,:,0])
        n_z = len(fld[:,0,0])
        
        for i_x1 in range(n_x):
            for i_y1 in range(n_y):                      
                    for i_x2 in range(n_x):
                        for i_y2 in range(n_y):
                            j = 0
                            for k in range(n_z):
                                j += (fld[k, i_y1, i_x1] * fld[k, i_y2, i_x2].conjugate())
                            if norm:
                                AbsE1 = 0
                                AbsE2 = 0
                                for k in range(n_z):
                                    AbsE1 += abs(fld[k, i_y1, i_x1])
                                    AbsE2 += abs(fld[k, i_y2, i_x2])
                                J[i_y1, i_x1, i_y2, i_x2] = j / (AbsE1 * AbsE2 / n_z**2) / n_z
                            else:
                                J[i_y1, i_x1, i_y2, i_x2] = j / n_z
    
def gauss_fit(X, Y):
    import numpy as np
    import scipy.optimize as opt

    def gauss(x, p):  # p[0]==mean, p[1]==stdev p[2]==peak
        return p[2] / (p[1] * np.sqrt(2 * np.pi)) * np.exp(-(x - p[0])**2 / (2 * p[1]**2))

    p0 = [0, np.max(X) / 2, np.max(Y)]
    errfunc = lambda p, x, y: gauss(x, p) - y
    p1, success = opt.leastsq(errfunc, p0[:], args=(X, Y))
    fit_mu, fit_stdev, ampl = p1
    Y1 = gauss(X, p1)
    RMS = fit_stdev
    return (Y1, RMS)
