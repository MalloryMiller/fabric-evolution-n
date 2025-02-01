# N. M. Rathmann <rathmann@nbi.ku.dk>, 2020-2023

import sys, os, code # code.interact(local=locals())

import numpy as np
import scipy.special as sp
from netCDF4 import Dataset

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import constants as sfconst

#----------------------
# Experiment selection
#----------------------

    
exp     = None
exptype = None
ijstr   = None

if len(sys.argv) == 2: 
    exp     = sys.argv[1]
    exptype = exp[0:2]
    ijstr   = exp[-2:]

# DDRX values
Gamma = 4 # between 1 and 5
lambd = 0.15 # https://nam12.safelinks.protection.outlook.com/?url=https%3A%2F%2Fdoi.org%2F10.1016%2Fj.epsl.2020.116718&data=05%7C02%7Cmillemal%40iu.edu%7Ccca8159aade6400c8c4308dd3c9a22fe%7C1113be34aed14d00ab4bcdd02510be91%7C0%7C0%7C638733354878925519%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=vTp%2FaTh44HHpY1OiGwbYXNLLgXfznvvike7pDPEwgy8%3D&reserved=0


def generate_solution(exp, exptype, ijstr, Gamma=Gamma, lambd=lambd, timesteps = 50, truncation = 8, deformation = None):

    ii_to_axis  = lambda ij: 0 if (ij=='xx') else (1 if (ij=='yy') else 2)
    ij_to_plane = lambda ij: 0 if (ij=='yz') else (1 if (ij=='xz') else 2)

    if exptype == 'ue': # Uniaxial extension
        mod = dict(type='ps', T=-1, r=0,  axis=ii_to_axis(ijstr))
        strain_target = 3

    if exptype == 'uc': # Uniaxial compression
        mod = dict(type='ps', T=1, r=0,  axis=ii_to_axis(ijstr))
        strain_target = -0.98

    if exptype == 'cc': # Confined compression
        mod = dict(type='ps', T=1, r=+1, axis=ii_to_axis(ijstr))
        strain_target = -0.98

    if exptype == 'ss': # Simple shear
        mod = dict(type='ss', T=1, plane=ij_to_plane(ijstr))
        strain_target = 70 # in degrees
        
    if exptype == 'rr': # Ridgid rotation
        mod = dict(type='rr', T=1, plane=ij_to_plane(ijstr))
        strain_target = 90 # in degrees


    if deformation:
        strain_target = deformation
    
    strain_to_use = strain_target
    if exptype == 'ss' or exptype == 'rr': 
        strain_to_use = np.deg2rad(strain_target)






    #----------------------
    # Model integration
    #----------------------

    Nt = timesteps # Number of time steps
    L  = truncation # Spectral truncation, no higher than 20

    lm, nlm_len = sf.init(L)
    nlm, F, time, ugrad = sfint.lagrangianparcel(sf, mod, strain_to_use, Nt=Nt, iota=1, nu=1, Gamma0=Gamma, Lambda=lambd)

    #----------------------
    # Determine eigenvalues, principal directions, and enhancement factors
    #----------------------

    # Empty structure to fill
    vecdim = (Nt+1,3)
    eigvals  = np.zeros(vecdim)
    m1,m2,m3 = np.zeros(vecdim),np.zeros(vecdim),np.zeros(vecdim)
    p1,p2,p3 = np.zeros(vecdim),np.zeros(vecdim),np.zeros(vecdim)
    vecdim = (Nt+1,6)
    Eij_lin, Eij_nlin   = np.zeros(vecdim), np.zeros(vecdim)
    Epij_lin, Epij_nlin = np.zeros(vecdim), np.zeros(vecdim)

    # Grain parameters
    (Eij_grain_lin,  alpha_lin,  n_grain_lin)  = sfconst.ice['viscoplastic']['linear']
    (Eij_grain_nlin, alpha_nlin, n_grain_nlin) = sfconst.ice['viscoplastic']['nonlinear']

    viscparams = (Eij_grain_lin, alpha_lin, n_grain_lin)

    for tt in np.arange(0,Nt+1):

        c = nlm[tt,:]
        
        m1[tt,:],m2[tt,:],m3[tt,:], eigvals[tt,:] = sf.frame(c, 'e')
        p1[tt,:],p2[tt,:],p3[tt,:], _             = sf.frame(c, 'p')

        # Linear (n'=1) mixed Taylor--Sachs enhancements            
        Eij_lin[tt,:]  = sf.Eij_tranisotropic(c, m1[tt,:],m2[tt,:],m3[tt,:], Eij_grain_lin, alpha_lin, n_grain_lin)
        Epij_lin[tt,:] = sf.Eij_tranisotropic(c, p1[tt,:],p2[tt,:],p3[tt,:], Eij_grain_lin, alpha_lin, n_grain_lin)
        
        # Nonlinear (n'=3) Sachs enhancements
        Eij_nlin[tt,:]  = sf.Eij_tranisotropic(c, m1[tt,:],m2[tt,:],m3[tt,:], Eij_grain_nlin, alpha_nlin, n_grain_nlin)
        Epij_nlin[tt,:] = sf.Eij_tranisotropic(c, p1[tt,:],p2[tt,:],p3[tt,:], Eij_grain_nlin, alpha_nlin, n_grain_nlin)
        
    #----------------------
    # Save
    #----------------------

    fname = generate_file_name(exp, Gamma, lambd, strain_target)

    os.makedirs("/".join (fname.split("/")[:-1]), exist_ok=True)

    ncfile = Dataset(fname, mode='w',format='NETCDF3_CLASSIC') 

    # Config
    dt = time[1]-time[0]
    ncfile.tsteps, ncfile.dt, ncfile.L = Nt, dt, L
    ncfile.Ecc_lin,  ncfile.Eca_lin,  ncfile.alpha_lin  = Eij_grain_lin[0],  Eij_grain_lin[1],  alpha_lin
    ncfile.Ecc_nlin, ncfile.Eca_nlin, ncfile.alpha_nlin = Eij_grain_nlin[0], Eij_grain_nlin[1], alpha_nlin

    ncfile.ugrad = ugrad.flatten()

    # Dimensions
    c_did       = ncfile.createDimension('DOF',     nlm_len)
    time_did    = ncfile.createDimension('tstep',   Nt+1)
    eig_did     = ncfile.createDimension('eigval',  3)
    dim_did     = ncfile.createDimension('dim',     3)
    dim6_did    = ncfile.createDimension('dim6',    6)
    pair_did    = ncfile.createDimension('pair',    2)

    # Variables
    myint = np.int32
    myflt = np.float64
    dimarr_vec = ('lm','lm','lmdyn')

    f_lm   = ncfile.createVariable('lm', myint, ('DOF','pair')) 
    f_c_re = ncfile.createVariable('c_re', myflt, ('tstep','DOF')) 
    f_c_im = ncfile.createVariable('c_im', myflt, ('tstep','DOF')) 
    f_eigvals = ncfile.createVariable('eigvals', myflt, ('tstep','eigval'))

    f_Eij_lin   = ncfile.createVariable('Eij_lin',   myflt, ('tstep','dim6'))
    f_Epij_lin  = ncfile.createVariable('Epij_lin',  myflt, ('tstep','dim6'))
    f_Eij_nlin  = ncfile.createVariable('Eij_nlin',  myflt, ('tstep','dim6'))
    f_Epij_nlin = ncfile.createVariable('Epij_nlin', myflt, ('tstep','dim6'))

    mkvec = lambda field: ncfile.createVariable(field, myflt, ('tstep','dim'))
    f_m1, f_m2, f_m3 = mkvec('m1'),  mkvec('m2'),  mkvec('m3')
    f_p1, f_p2, f_p3 = mkvec('p1'),  mkvec('p2'),  mkvec('p3')

    f_lm[:,:], f_c_re[:,:], f_c_im[:,:] = lm.T, np.real(nlm), np.imag(nlm)
    f_Eij_lin[:,:],  f_Eij_nlin[:,:]  = Eij_lin,  Eij_nlin
    f_Epij_lin[:,:], f_Epij_nlin[:,:] = Epij_lin, Epij_nlin
    f_m1[:,:], f_m2[:,:], f_m3[:,:] = m1, m2, m3
    f_p1[:,:], f_p2[:,:], f_p3[:,:] = p1, p2, p3
    f_eigvals[:,:] = eigvals
    ncfile.close(); 
    print('Solution dumped in solutions/%s'%fname)
    print('Plot result: python3 plot-fabric-evolution-latrot.py %s'%(exp))
    
    return fname


def generate_file_name(exp, Gamma, lambd, strain_target):
    folders = f'{exp}'
    if Gamma != None:
        folders += f'/gam{reduce(Gamma)}'
    if lambd != None:
        folders += f'/lam{reduce(lambd)}'

    os.makedirs("solutions/" + folders, exist_ok=True)

    fname = folders +  f'/target{reduce(strain_target)}' + '.nc'
    return "solutions/" + fname

def reduce(n):
    if (n%1 == 0):
        return int(n)
    else:
        return n
    

if len(sys.argv) == 2:
    generate_solution(exp, exptype, ijstr, Gamma, lambd)

