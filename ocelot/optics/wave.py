'''
wave optics
'''

from numpy import sin, cos, pi, sqrt, log, exp, array, random, sign
from numpy.linalg import norm
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
import scipy.integrate as integrate
#import matplotlib.animation as animation

from ocelot.optics.elements import *
from ocelot.common.globals import *


class TransferFunction(object):
    '''
    data container for Fourier Optics transfer functions
    '''
    def __init__(self):
        self.k = None # wave vector - 2*pi/wavelength
        self.tr = None # complex value of transmission - modulus*exp(-i*phase)
        self.ref = None # .. of reflection
        self.xlamds = None # carrier wavelength
        self.mid_k = None # center of feature in spectrum
        self.dk = None # width of feature in spectrum
    
    def ev(self):
        return self.k* h_eV_s/2/pi * speed_of_light
    
    def __mul__(self, f):
        if f.__class__ == TransferFunction:
            f2 = TransferFunction()
            f2.k = f.k
            f2.ev = f.ev
            # TODO check data grid alignment
            
            f2.tr = self.tr * f.tr
            f2.ref = self.ref * f.ref
            
            return f2
        return None


class WaveFront:
    def __init__(self):
        pass

class Scene:
    def __init__(self):
        pass


def init():
    global scene
    scene.line_wf.set_data([], [])
    scene.time_text.set_text('')
    #scene.profile_im.set_data(np.ones([51,51])*10)
    res = []
    if 'geometry' in scene.views:
        res.append(scene.line_wf)
        res.append(scene.time_text)
    if 'detectors' in scene.views:
        res.append(scene.profile_im)
    return res


def normal(x1,x2):
    d = (x2[0] - x1[0]) / (x2[1] - x1[1])
    #print 'x1,x2,d=', x1, x2, d 
    n1, n2 = 1.0 / np.sqrt(1+d**2), -d / np.sqrt(1+d**2)
    return np.array([n1, n2])


def rotate_pi(v,n):
    vrot = -v + 2*n*np.dot(n,v)
    return vrot 



'''
data structures for optical field propagation
'''
class Mesh:
    def __init__(self, nx, ny, dtype=np.float):
        self.nx = nx
        self.ny = ny
        
        self.points = np.zeros([nx,ny], dtype=dtype)
                
        self.x = 0 
        self.y = 0 
         
        self.dx = 1.0 
        self.dy = 1.0 
        
    def __getitem__(self, idx):
        return self.points[idx]
    def __setitem__(self, idx, val):
        self.points[idx] = val
        
    def __str__(self):
        s = "Mesh " + str(self.nx) + 'x' + str(self.ny) + ' '
        s += 'xmin='+ str(self.x) + ' xmax=' + str(self.x + self.dx * (self.nx - 1) ) + ' '  
        s += 'ymin='+ str(self.y) + ' ymax=' + str(self.y + self.dy * (self.ny - 1) ) + ' '
        s+= '\n' + str(self.points)
        return s
    def idx(self, x):
        '''
        return mesh point idx of left-bottom corner of the cell a coordinate belongs to
        if coordinate outside mesh return -1
        '''
        ix = (x[0] - self.x) / self.dx
        if ix < 0 or ix > self.nx - 1:
            ix = -1

        iy = (x[1] - self.y) / self.dy
        if iy < 0 or iy > self.ny - 1:
            iy = -1

        return ix,iy
    
    def init(self, f = lambda x, y : 0):
        
        x=0; y=0
        
        for i1 in range(self.nx):
            for i2 in range(self.ny):
                x = self.x + self.dx * i1
                y = self.y + self.dy * i2
                self[i1,i2] = f(x,y)
                


class ParaxialFieldSlice():
    '''
    complex transverse electric field E_x + i E_y
    '''
    def __init__(self, lam=1.0, nx=31, ny=31, size_x=1.0, size_y = 1.0):
        '''
        lam -- wavelength (m)
        '''
        c = 1.0
        self.lam = lam     
        self.k = 2*pi / self.lam 
        self.w = self.k / c 
        
        self.nx = nx
        self.ny = ny

        self.size_x = size_x
        self.size_y = size_y
        
        self.x = np.zeros(nx)
        self.y = np.zeros(ny)
            
    def init_field(self, f = lambda x, y : 0):
        self.mesh = Mesh(nx=self.nx, ny=self.ny, dtype = np.complex)
        self.mesh.x = -self.size_x
        self.mesh.y = -self.size_y
        
        self.mesh.dx = 2*(-self.mesh.x) / ( self.mesh.nx -1)
        self.mesh.dy = 2*(-self.mesh.y) / ( self.mesh.ny -1)
        
        x=0; y=0
        
        for i1 in range(self.mesh.nx):
            for i2 in range(self.mesh.ny):
                x = self.mesh.x + self.mesh.dx * i1
                y = self.mesh.y + self.mesh.dy * i2
                self.mesh[i1,i2] = f(x,y)
                
                self.x[i1] = x
                self.y[i2] = y
                
                #print i1, i2, x, y, self.mesh[i1,i2] 

    def __getitem__(self, idx):
        return self.mesh[idx]
    def __setitem__(self, idx, val):
        self.mesh[idx] = val

def rescale(of, scale=2.0):
    of.size_x /= scale
    of.size_y /= scale
        
    of.x /= scale
    of.y /= scale

    of.mesh.x /= scale
    of.mesh.y /= scale

    of_old = np.copy(of.mesh.points)
    for i in range(of.nx):
        for j in range(of.ny):
            i_new = int(i*scale)
            j_new = int(j*scale)
            try:
                of[i,j] = of_old[ int(i*scale). int(j*scale)]
            except:
                of[i,j] = 0
             
            
            
def propagate_fourier(of, dz, obj=None, scale=1.0):
    '''
    wave propagator
    '''
    
    if obj == None or obj.__class__ == OptDrift:
        debug('wave propagator: drift')
        spec = fft.fft2(of[:,:])
        
        kx = np.fft.fftfreq(of.nx, d=2*of.size_x/of.nx)
        ky = np.fft.fftfreq(of.ny, d=2*of.size_y/of.ny)
        
        for i in range(of.nx):
            for j in range(of.ny):
                k = 2*pi / of.lam #of.w / c
                #print (kx[i]/k), (ky[j]/k)
                #phi = k * sqrt(1 - (kx[i]/k)**2 - (ky[j]/k)**2)
                phi = -pi * of.lam * ( (kx[i])**2 + (ky[j])**2 )
                #print phi*dz
                spec[i,j] *= exp(1j * phi*dz + 1j *k*dz)
            
        of[:,:] = fft.ifft2(spec)
        
    if obj.__class__ == Aperture:
        debug('wave propagator: aperture', obj.d)
        for i in range(of.nx):
            for j in range(of.ny):
                #print of.x[i], obj.d[0]
                if (of.x[i]/obj.d[0])**2 + (of.y[j]/obj.d[1])**2 >1:
                    of[i,j] = 0 

    if obj.__class__ == Lense:
        debug('wave propagator: lense, f=', obj.f, " [m]")
        for i in range(of.nx):
            for j in range(of.ny):
                phi = pi*( (of.x[i]/sqrt(of.lam*obj.f))**2 + (of.y[j]/sqrt(of.lam*obj.f))**2 )
                of[i,j] *= np.exp(-1j*phi) 
                #of[i,j] *= i
    

def propagate_fresnel(of, dz, scale=1.0):
    '''
    Propagate paraxial field slice in free space, Fresnel 
    '''
    k = 2*pi/ of.lam

    #of_old = np.copy(of.mesh.points)

    for i in range(of.nx):
        for j in range(of.ny):
            tmp = 0.0 + 0.0j
            print(i,j)
            for i1 in range(of.nx):
                for j1 in range(of.ny):
                    phi = 1j * k * ( (of.x[i1] - of.x[i])**2 + (of.y[j1] - of.y[j])**2 ) / (2.0 * dz)
                    #print phi
                    tmp = tmp +  of[i1,j1] * exp(phi) / (of.nx * of.ny)
            of[i,j] = tmp * exp(1j*k*dz) / (1j * of.lam * dz)
            print(of[i,j])
        
#---------------- GENESIS v4 and other functions ----------------- HMCC
def calc_ph_sp_dens(spec, freq_ev, n_photons, spec_squared=1):
    '''
    calculates number of photons per electronvolt
    '''
    if spec.ndim == 1:
        axis=0
    else:
        if spec.shape[0] == freq_ev.shape[0]:
            spec = spec.T
        axis=1
            # axis=0
        # elif spec.shape[1] == freq_ev.shape[0]:
            # axis=1
        # else:
            # raise ValueError('operands could not be broadcast together with shapes ', spec.shape, ' and ', freq_ev.shape)
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    if spec_squared:
        spec_sum = np.trapz(spec, x=freq_ev, axis=axis)
    else:
        spec_sum = np.trapz(abs(spec)**2, x=freq_ev, axis=axis)
        
    if np.size(spec_sum) == 1:
        if spec_sum == 0: #avoid division by zero
            spec_sum = np.inf 
    else:
        spec_sum[spec_sum == 0] = np.inf #avoid division by zero
        
    if spec_squared:
        norm_factor = n_photons / spec_sum
    else:
        norm_factor = np.sqrt(n_photons / spec_sum)
    
    if spec.ndim == 2:
        norm_factor = norm_factor[:,np.newaxis]
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    # _logger.debug('norm_factor.shape = {}'.format(norm_factor.shape))
    spec = spec * norm_factor
    if axis==1:
        spec = spec.T
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    return spec
    



def imitate_1d_sase_like(td_scale, td_env, fd_scale, fd_env, td_phase = None, fd_phase = None, phen0 = None, en_pulse = None, fit_scale = 'td', n_events = 1):
    '''
    Models FEL pulse(s) based on Gaussian statistics
    td_scale - scale of the pulse on time domain [m]
    td_env - expected pulse envelope in time domain [W] fd_scale - scale of the pulse in frequency domain [eV]
    fd_env - expected pulse envelope in frequency domain [a.u.]
    td_phase - additional phase chirp to be added in time domain
    fd_phase - additional phase chirp to be added in frequency domain
    phen0 - sampling wavelength expressed in photon energy [eV]
    en_pulse - expected average energy of the pulses [J]
    fit_scale - defines the scale in which outputs should be returned:
        'td' - time domain scale td_scale is used for the outputs, frequency domain phase and envelope will be re-interpolated
        'fd' - frequency domain scale fd_scale is used for the outputs, time domain phase and envelope will be re-interpolated
    n_events - number of spectra to be generated
        
    returns tuple of 4 arguments: (ph_en, fd, s, td)
    fd_scale - colunm of photon energies in eV
    fd - matrix of radiation in frequency domain with shape, normalized such that np.sum(abs(fd)**2) is photon spectral density, i.e: np.sum(abs(fd)**2)*fd_scale = N_photons
    td - matrix of radiation in time domain, normalized such that abs(td)**2 = radiation_power in [w]
    '''
    
    _logger.info('generating 1d radiation field imitating SASE')
    
    if fit_scale == 'td':
        
        n_points = len(td_scale)
        s = td_scale
        Ds = (td_scale[-1] - td_scale[0])
        ds = Ds / n_points
        
        td = np.random.randn(n_points,n_events) + 1j * np.random.randn(n_points,n_events)
        td *= np.sqrt(td_env[:, np.newaxis])
        fd = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(td, axes=0), axis=0), axes=0) 
        # fd = np.fft.ifft(td, axis=0)
        # fd = np.fft.fftshift(fd, axes=0)
        
        if phen0 is not None:
            e_0 = phen0
        else:
            e_0 = np.mean(fd_scale)
        
        # internal interpolated values
        fd_scale_i = h_eV_s * np.fft.fftfreq(n_points, d=(ds / speed_of_light)) + e_0 # internal freq.domain scale based on td_scale
        fd_scale_i = np.fft.fftshift(fd_scale_i, axes=0)
        fd_env_i = np.interp(fd_scale_i,fd_scale,fd_env, right=0, left=0)
        
        if fd_phase is None:
            fd_phase_i = np.zeros_like(fd_env_i)
        else:
            fd_phase_i = np.interp(fd_scale_i,fd_scale,fd_phase, right=0, left=0)
        
        fd *= np.sqrt(fd_env_i[:, np.newaxis]) * np.exp(1j * fd_phase_i[:, np.newaxis])
        
        # td = np.fft.ifftshift(fd, axes=0) 
        # td = np.fft.fft(td, axis=0)
        td = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(fd, axes=0), axis=0), axes=0)
        
        td_scale_i = td_scale
        
    elif fit_scale == 'fd':
        
        n_points = len(fd_scale)
        Df = abs(fd_scale[-1]-fd_scale[0]) / h_eV_s
        df = Df / n_points
             
        fd = np.random.randn(n_points,n_events) + 1j * np.random.randn(n_points,n_events)
        fd *= np.sqrt(fd_env[:, np.newaxis])
        td = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(fd, axes=0), axis=0), axes=0)
        
        td_scale_i = np.fft.fftfreq(n_points, d=df) * speed_of_light
        td_scale_i = np.fft.fftshift(td_scale_i, axes=0)
        td_scale_i -= np.amin(td_scale_i)
        td_env_i = np.interp(td_scale_i, td_scale, td_env, right=0, left=0)
        
        if td_phase is None:
            td_phase_i = np.zeros_like(td_env_i)
        else:
            td_phase_i = np.interp(td_scale_i, td_scale, td_phase, right=0, left=0)
        
        td *= np.sqrt(td_env_i[:, np.newaxis]) * np.exp(1j * td_phase_i[:, np.newaxis])

        fd = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(td, axes=0), axis=0), axes=0) 

        fd_scale_i = fd_scale
        
    else:
        raise ValueError('fit_scale should be either "td" of "fd"')
    
    #normalization for pulse energy
    if en_pulse == None:
        en_pulse = np.trapz(td_env_i, td_scale_i / speed_of_light) # CALCULATE FOR fit_scale == 'td' !!!!!!!!!!
    pulse_energies = np.trapz(abs(td)**2, td_scale_i / speed_of_light, axis=0)
    scale_coeff = en_pulse / np.mean(pulse_energies)
    td *= np.sqrt(scale_coeff)

    #normalization for photon spectral density
    spec = np.mean(np.abs(fd)**2, axis=1)
    spec_center = np.sum(spec * fd_scale_i) / np.sum(spec)
    
    n_photons = pulse_energies * scale_coeff / q_e / spec_center
    fd = calc_ph_sp_dens(fd, fd_scale_i, n_photons, spec_squared=0)
    td_scale, fd_scale = td_scale_i, fd_scale_i
    
    return (td_scale, td, fd_scale, fd)
    

def imitate_1d_sase(spec_center = 500, spec_res = 0.01, spec_width = 2.5, spec_range = (None,None), pulse_length = 6, en_pulse = 1e-3, flattop = 0, n_events = 1, spec_extend = 5):
    '''
    Models FEL pulse(s) based on Gaussian statistics
    spec_center - central photon energy in eV
    spec_res - spectral resolution in eV
    spec_width - width of spectrum in eV (fwhm of E**2)
    spec_range = (E1, E2) - energy range of the spectrum. If not defined, spec_range = (spec_center - spec_width * spec_extend, spec_center + spec_width * spec_extend)
    pulse_length - longitudinal size of the pulse in um (fwhm of E**2)
    en_pulse - expected average energy of the pulses in Joules
    flattop - if true, flat-top pulse in time domain is generated with length 'pulse_length' in um
    n_events - number of spectra to be generated

    return tuple of 4 arguments: (s, td, ph_en, fd)
    ph_en - colunm of photon energies in eV with size (spec_range[2]-spec_range[1])/spec_res
    fd - matrix of radiation in frequency domain with shape ((spec_range[2]-spec_range[1])/spec_res, n_events), normalized such that np.sum(abs(fd)**2) is photon spectral density, i.e: np.sum(abs(fd)**2)*spec_res = N_photons
    s - colunm of longitudinal positions along the pulse in yime domain in um
    td - matrix of radiation in time domain with shape ((spec_range[2]-spec_range[1])/spec_res, n_events), normalized such that abs(td)**2 = radiation_power
    '''
    
    
    if spec_range == (None,None):
        spec_range = (spec_center - spec_width*spec_extend, spec_center + spec_width*spec_extend)
    elif spec_center == None:
        spec_center = (spec_range[1] + spec_range[0]) / 2
    
    pulse_length_sigm = pulse_length / (2*np.sqrt(2*np.log(2)))
    spec_width_sigm = spec_width / (2*np.sqrt(2*np.log(2)))
    
    fd_scale = np.arange(spec_range[0], spec_range[1], spec_res)
    n_points = len(fd_scale)
    _logger.debug(ind_str + 'N_points * N_events = %i * %i' %(n_points, n_events))

    fd_env = np.exp(-(fd_scale - spec_center)**2 / 2 / spec_width_sigm**2)
    td_scale = np.linspace(0, 2*np.pi / (fd_scale[1] - fd_scale[0]) * hr_eV_s * speed_of_light, n_points)
    
    if flattop:
        td_env = np.zeros_like(td_scale)
        il = find_nearest_idx(td_scale, np.mean(td_scale)-pulse_length * 1e-6 / 2)
        ir = find_nearest_idx(td_scale, np.mean(td_scale)+pulse_length * 1e-6 / 2)
        td_env[il:ir]=1
    else:
        s0 = np.mean(td_scale)
        td_env = np.exp(-(td_scale-s0)**2 / 2 / (pulse_length_sigm * 1e-6)**2)
        
    result = imitate_1d_sase_like(td_scale, td_env, fd_scale, fd_env, phen0 = spec_center, en_pulse = en_pulse, fit_scale = 'fd', n_events = n_events)
    
    return result
