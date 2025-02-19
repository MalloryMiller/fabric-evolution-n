import numpy as np
import os

from netCDF4 import Dataset

from specfabpy import specfab as sf
from specfabpy import integrator as sfint
from specfabpy import constants as sfconst


from specfabpy import plotting as sfplt
from scipy.spatial.transform import Rotation
FS = sfplt.setfont_tex()




GAMMA_FORMULA = [ 0.188, 5.94] # [0] + [1]temp


def apply_formula(formula, x):
    m = formula[0]
    c = formula[1]
    return (x * m) + c


class Experiment():

    

    def __init__(self, exptype, ijstr, temp = None, lambd = None, timesteps = 50, truncation = 8, deformation = None):
        self.fname = None
        self.exptype = exptype
        
        self.temp = temp
        self.lamd = lambd
        self.timesteps = timesteps
        self.truncation = truncation


        self.Gamma = None
        self.lamd = None
        if temp != None:
            self.Gamma = round(apply_formula(GAMMA_FORMULA, self.temp), 3)
            if lambd != None and lambd != False:
                self.lamd = 0.15


        ii_to_axis  = lambda ij: 0 if (ij=='xx') else (1 if (ij=='yy') else 2)
        ij_to_plane = lambda ij: 0 if (ij=='yz') else (1 if (ij=='xz') else 2)

        self.mod = -1
        self.strain_target = -1

        if self.exptype == 'ue': # Uniaxial extension
            self.mod = dict(type='ps', T=-1, r=0,  axis=ii_to_axis(ijstr))
            self.strain_target = 5

        if self.exptype == 'uc': # Uniaxial compression
            self.mod = dict(type='ps', T=1, r=0,  axis=ii_to_axis(ijstr))
            self.strain_target = -0.98

        if self.exptype == 'cc': # Confined compression
            self.mod = dict(type='ps', T=1, r=+1, axis=ii_to_axis(ijstr))
            self.strain_target = -0.98

        if self.exptype == 'ss': # Simple shear
            self.mod = dict(type='ss', T=1, plane=ij_to_plane(ijstr))
            self.strain_target = 70 # in degrees
            
            
        if self.exptype == 'rr': # Ridgid rotation
            self.mod = dict(type='rr', T=1, plane=ij_to_plane(ijstr))
            self.strain_target = 90 # in degrees


        if deformation:
            self.strain_target = deformation
        
        self.strain_to_use = self.strain_target
        if exptype == 'ss' or exptype == 'rr': 
            self.strain_to_use = np.deg2rad(self.strain_target)



        pass



    def init_file_name(self):
        folders = f'{self.exptype}'
        if self.Gamma and self.temp:
            folders += f'/temp{reduce(self.temp)}'
        if self.lamd != None:
            folders += f'/lam'

        os.makedirs("solutions/" + folders, exist_ok=True)

        self.fname = "solutions/" + folders +  f'/L{self.truncation}' + '.nc'
        return  self.fname
    


    def generate_file(self, remake = False, rectify = False):

        self.init_file_name()


        if not remake and os.path.isfile(self.fname):
            print(self.fname, "already exists. Run with remake=True to overwrite.")
            return False

        
        #----------------------
        # Model integration
        #----------------------

        Nt = self.timesteps # Number of time steps
        L  = self.truncation # Spectral truncation, no higher than 20

        lm, nlm_len = sf.init(L)
        nlm, F, time, ugrad = sfint.lagrangianparcel(sf, self.mod, self.strain_to_use, Nt=self.timesteps, 
                                                     iota=1, nu=1, Gamma0=self.Gamma, Lambda=self.lamd)
        
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


        for tt in np.arange(0,Nt+1):

            c = nlm[tt,:]
            #c[:] = sf.rotate_nlm(sf.rotate_nlm(c[:], -np.pi/3, 0), 0 ,0)
            
            m1[tt,:],m2[tt,:],m3[tt,:], eigvals[tt,:] = sf.frame(c, 'e')
            p1[tt,:],p2[tt,:],p3[tt,:], _             = sf.frame(c, 'p')  # TODO shear
            

            if (rectify):
                m1[tt],m2[tt],m3[tt] = angle_snap(m1[tt,:], m2[tt,:], m3[tt,:])
            


            # Linear (n'=1) mixed Taylor--Sachs enhancements            
            Eij_lin[tt,:]  = sf.Eij_tranisotropic(c, m1[tt,:],m2[tt,:],m3[tt,:], Eij_grain_lin, alpha_lin, n_grain_lin)
            Epij_lin[tt,:] = sf.Eij_tranisotropic(c, p1[tt,:],p2[tt,:],p3[tt,:], Eij_grain_lin, alpha_lin, n_grain_lin)
            
            # Nonlinear (n'=3) Sachs enhancements
            Eij_nlin[tt,:]  = sf.Eij_tranisotropic(c, m1[tt,:],m2[tt,:],m3[tt,:], Eij_grain_nlin, alpha_nlin, n_grain_nlin)
            Epij_nlin[tt,:] = sf.Eij_tranisotropic(c, p1[tt,:],p2[tt,:],p3[tt,:], Eij_grain_nlin, alpha_nlin, n_grain_nlin)
        

        #----------------------
        # Save
        #----------------------

        ncfile = Dataset(self.fname, mode='w',format='NETCDF3_CLASSIC') 

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

        return True
    


    def load_solution(self, specific_steps = None, remake = False, isolate = False):
        self.init_file_name()


        if not os.path.isfile(self.fname):
            print(self.fname, "does not exist. Generate it using generate_file first.")
            return

        if os.path.isdir("solutions/frames/{}".format(self.fname[10:-3])) and not remake:
            print(self.fname, "already has an image folder. Run with remake=True to overwrite.")
            return


        fh = Dataset('%s'%(self.fname), mode='r')
        loadvar = lambda field: np.array(fh.variables[field][:])

        # Model config
        Nt, dt, L = fh.getncattr('tsteps'), fh.getncattr('dt'), fh.getncattr('L')

        # Grain parameters
        Eca_lin,  Ecc_lin  = fh.getncattr('Eca_lin'),  fh.getncattr('Ecc_lin')
        Eca_nlin, Ecc_nlin = fh.getncattr('Eca_nlin'), fh.getncattr('Ecc_nlin')
        alpha_lin, alpha_nlin = fh.getncattr('alpha_lin'), fh.getncattr('alpha_nlin')

        # CPO state
        lm, c = loadvar('lm'), loadvar('c_re') + 1j*loadvar('c_im') 
        lm = np.array(lm).T
        eigvals = loadvar('eigvals')
        m1,m2,m3 = loadvar('m1'), loadvar('m2'), loadvar('m3')
        p1,p2,p3 = loadvar('p1'), loadvar('p2'), loadvar('p3')

        # Enhancement factors
        Eij_lin, Eij_nlin = loadvar('Eij_lin'), loadvar('Eij_nlin') 
        Epij_lin, Epij_nlin = loadvar('Epij_lin'), loadvar('Epij_nlin') 

        #------------------
        # Plot
        #------------------

        tsteps = specific_steps # time steps to plot
        if specific_steps == None:
            tsteps = list(range(0, self.timesteps + 1)) # time steps to plot

        lw0, lw1, lw2 = 2.5,2.25,2.0

        import matplotlib.pyplot as plt
        from matplotlib import rcParams, rc
        import matplotlib.ticker as mticker
        import matplotlib.gridspec as gridspec
        import cartopy.crs as ccrs
        from matplotlib.offsetbox import AnchoredText
            
        rotation=55 
        inclination=45
        geo, prj = sfplt.getprojection(rotation=rotation, inclination=inclination)

        def plot_vec(ax, v, lbl, color, ls='-', lw=2):
            
            ax.plot([0, +v[0]],[0, +v[1]],[0,+v[2]], color=color, ls=ls, lw=lw, label=lbl)
            ax.plot([0, -v[0]],[0, -v[1]],[0,-v[2]], color=color, ls=ls, lw=lw)

        for tt in tsteps:

            #----------------------
            # Figure setup
            #----------------------

            dpi, scale = 200, 3.3
            fig = plt.figure(figsize=(3/2*2.1*scale,2.3*scale))
            gs = gridspec.GridSpec(2,3, height_ratios=[1,1.2], width_ratios=[1,1,1])
            gs.update(left=-0.03, right=1-0.06/3, top=0.97, bottom=0.20, wspace=0.015*18, hspace=0.35)

            ax_ODF     = plt.subplot(gs[0, 0], projection=prj)
            ax_mi      = plt.subplot(gs[0, 1], projection='3d')
            ax_eigvals = plt.subplot(gs[0, 2])
            ax_Elin    = plt.subplot(gs[1, 1])
            ax_Enlin   = plt.subplot(gs[1, 2])
            
            ax_ODF.set_global() # be sure to show entire S^2
                                    
            # CORRECT ROTATION (SHEAR)

            vertical_component = get_vertical_component(m1[tt],m2[tt],m3[tt])


            #Eij_lin[tt,:] = sf.Eij_tranisotropic(c, m1[tt,:], m2[tt,:], m3[tt,:], Eij_grain_lin, alpha_lin, n_grain_lin)
            
            # Nonlinear (n'=3) Sachs enhancements
            #Eij_nlin[tt,:]  = sf.Eij_tranisotropic(c, m1[tt,:], m2[tt,:], m3[tt,:], Eij_grain_nlin, alpha_nlin, n_grain_nlin)
           
            #----------------------
            # ODF (orb)
            #----------------------

            ax = ax_ODF
            sfplt.plotODF(c[tt,:], lm, ax, cmap='Greys')
            sfplt.plotcoordaxes(ax, geo, axislabels='vuxi')

            #----------------------
            # Eigenvalues (top right)
            #----------------------

            ax_eigvals.plot([self.get_pressure(tt),self.get_pressure(tt)],[0,1],':k', lw=2)
            
            steps = np.arange(len(eigvals[:,0]))
            steps = self.get_pressure(steps)

            if (vertical_component == 0 or not isolate):
                ax_eigvals.plot(steps,eigvals[:,0], '-', c=sfplt.c_red,   label='$a_{1}$', lw=lw0)
            if (vertical_component == 1 or  not isolate):
                ax_eigvals.plot(steps,eigvals[:,1], '-', c=sfplt.c_green, label='$a_{2}$', lw=lw1)
            if (vertical_component == 2 or  not isolate):
                ax_eigvals.plot(steps,eigvals[:,2], '-', c=sfplt.c_blue,  label='$a_{3}$', lw=lw2)
            
            ax_eigvals.set_ylim([0,1])
            ax_eigvals.set_xlim([0, self.get_pressure(Nt+1)])
            
            if self.exptype == "ss":
                ax_eigvals.set_xlabel('Target Angle')
            else:
                ax_eigvals.set_xlabel('Target Change')

            ax_eigvals.set_ylabel('$a_{i}$')
            ax_eigvals.grid()              
            ax_eigvals.legend(handlelength=1, ncol=1, labelspacing=0.3, fancybox=False)#, loc=2)
                    
            #----------------------
            # Principal frame (top middle)
            #----------------------

            ax_mi.view_init(elev=90-inclination, azim=rotation) # same as ODF plot
            
            ax_mi.set_xlabel('$x$'); ax_mi.set_xlim([-1,1])
            ax_mi.set_ylabel('$y$'); ax_mi.set_ylim([-1,1])
            ax_mi.set_zlabel('$z$'); ax_mi.set_zlim([0,1])
            
            plot_vec(ax_mi,m1[tt,:], r'$\vb{m}_1$', sfplt.c_red)
            plot_vec(ax_mi,m2[tt,:], r'$\vb{m}_2$', sfplt.c_green)
            plot_vec(ax_mi,m3[tt,:], r'$\vb{m}_3$', sfplt.c_blue)
            lwpq=1
            plot_vec(ax_mi,p1[tt,:], r'$\vb{p}_{1}$', sfplt.c_red,   ls='--', lw=lwpq)
            plot_vec(ax_mi,p2[tt,:], r'$\vb{p}_{2}$', sfplt.c_green, ls='--', lw=lwpq)
            plot_vec(ax_mi,p3[tt,:], r'$\vb{p}_{3}$', sfplt.c_blue,  ls='--', lw=lwpq)
            
            ax_mi.legend(handlelength=1, bbox_to_anchor=(1.17,1), fancybox=False)#, loc=1)

            #----------------------
            # Enhancement factors (bottom charts)
            #----------------------
            
            lblm = lambda ii,jj: '$E_{m_%i m_%i}$'%(ii+1,jj+1)
            lblp = lambda ii,jj: '$E_{p_%i p_%i}$'%(ii+1,jj+1)

            def plot_enhancements(ax, Eij, Epij):

                
                if not isolate or vertical_component == 0:
                    ax.semilogy(steps, Eij[:,0], '-', c=sfplt.c_red,   label=lblm(0,0), lw=lw0)
                if not isolate or vertical_component == 1:
                    ax.semilogy(steps, Eij[:,1], '-', c=sfplt.c_green, label=lblm(1,1), lw=lw1)
                if not isolate or vertical_component == 2:
                    ax.semilogy(steps, Eij[:,2], '-', c=sfplt.c_blue,  label=lblm(2,2), lw=lw2)    
                    
                if not isolate:
                    ax.semilogy(steps, Eij[:,3], ':', c=sfplt.c_lgreen, label=lblm(1,2), lw=lw2)
                    ax.semilogy(steps, Eij[:,4], ':', c=sfplt.c_lblue,  label=lblm(0,2), lw=lw1)
                    ax.semilogy(steps, Eij[:,5], ':', c=sfplt.c_lred,   label=lblm(0,1), lw=lw0)

                    Eratio = np.divide(Eij[:,5], Epij[:,5])
                    ax.semilogy(steps, Epij[:,5], '-', c=sfplt.c_gray, label=lblp(0,1), lw=lw2) 
                    ax.semilogy(steps, Eratio, '--', c=sfplt.c_gray, label=lblm(0,1)+'/'+lblp(0,1), lw=lw2)

                xlims=[0, self.get_pressure(Nt+1)]
                ax.semilogy(xlims, [2.5,2.5], '--k', lw=1) 
                ax.semilogy(xlims, [4.375,4.375], '--k', lw=1)
                ax.set_xlim(xlims)
                ax.set_ylim(1e-1, 1e1)
                #ax.set_ylim([np.amin([1e-1, np.amin(Eij[:]), np.amin(Epij[:])]), np.amax([1e+1, np.amax(Eij[:]), np.amax(Epij[:]), np.amax(Eratio)])])
                
                if self.exptype == "ss":
                    ax.set_xlabel('Target Angle')
                else:
                    ax.set_xlabel('Target Change')

                ax.set_ylabel('$E_{vw}$')
                ax.grid()
            
            plot_enhancements(ax_Enlin, Eij_nlin, Epij_nlin)
            ax_Enlin.set_title(r'Nonlinear Sachs')
            
            plot_enhancements(ax_Elin,  Eij_lin,  Epij_lin)
            ax_Elin.set_title(r'Linear mixed Taylor--Sachs')
            
            ax_Enlin.add_artist(AnchoredText('\n'.join( (r"$n' = %i$"%(3), r'$E_{cc} = %.1f$'%(Ecc_nlin), r'$E_{ca} = %.0e$'%(Eca_nlin), r'$\alpha = %.4f$'%(alpha_nlin)) ), loc='lower left', frameon=True,))
            ax_Elin.add_artist( AnchoredText('\n'.join( (r"$n' = %i$"%(1), r'$E_{cc} = %.1f$'%(Ecc_lin),  r'$E_{ca} = %.0e$'%(Eca_lin),  r'$\alpha = %.4f$'%(alpha_lin)) ),  loc='lower left', frameon=True,))
            
            ax_Enlin.legend(fontsize=FS+1, ncol=3, handlelength=2, columnspacing=1.2, labelspacing=0.3, fancybox=False, bbox_to_anchor=(1.05,-0.55), loc=4)   
            
            #----------------------
            # Model config
            #----------------------
            
            props = dict(boxstyle='square', facecolor='wheat', alpha=0.5)
            textstr1 = '\n'.join(( r'"%s"'%(self.fname.replace('_', '\_')), r'$L = %i$'%(L) ))
            ax_eigvals.text(-1.6, -0.1, textstr1, transform=ax_eigvals.transAxes, fontsize=FS, bbox=props)

            #----------------------
            # Save figure
            #----------------------
            
            os.makedirs("solutions/frames/{}".format(self.fname[10:-3]), exist_ok=True)
            fout = 'solutions/frames/%s/%f.png'%(self.fname[10:-3], reduce(self.get_pressure(tt)))
            print('Saving %s'%(fout))
            plt.savefig(fout, dpi=dpi)
            plt.close('all')



    def get_pressure(self, n):
        amount = float(self.strain_target) / float(self.timesteps)
        return n * amount




def reduce(n):
    '''
    If the given number is whole, returns as an integer with no floating point. Otherwise returns floored value
    '''
    if (n == 0):
        return 0 # it was returning -0 sometimes???
    else:
        return round(n, 1)
    


def angle_snap(e1, e2, e3):
    print(e1, e2, e3)
    a = [1,0,0] #top
    b = [-0,1,0]
    c = [0,0,1]
    og = np.array([e1, e2, e3])
    options = [a, b, c]
    order = []

    for x in range(len(og)):
        order.append(np.argmax(np.abs(og[:,x])))
        
    print([options[order[0]], options[order[1]], options[order[2]]])
    return [options[2], options[1], options[0]]



    
def get_vertical_component(m1,m2,m3):

    if (m1[2] == 1):
        return 0
    if (m2[2] == 1):
        return 1
    if (m3[2] == 1):
        return 2
    return None
