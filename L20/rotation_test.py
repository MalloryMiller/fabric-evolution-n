import numpy as np
from specfabpy import specfab as sf
import matplotlib.pyplot as plt
from specfabpy import plotting as sfplt
import matplotlib.gridspec as gridspec
from scipy.spatial.transform import Rotation
FS = sfplt.setfont_tex()
lm, nlm_len = sf.init(8) 

### Construct an arbitrary fabric to rotate
a2 = np.diag([0, 0, 1]) # arbitrary second-order structure tensor, a^(2)
nlm = np.zeros((nlm_len), dtype=np.complex64) # array of expansion coefficients
nlm[:sf.L2len] = sf.a2_to_nlm(a2) # l<=2 expansion coefficients of corresponding ODF

### Rotate ODF
# Note: assumes L=<12 (rotation for larger L is not implemented)
theta_d = '35'
phi_d = '10'
theta = np.deg2rad(int(theta_d))
phi   = np.deg2rad(int(phi_d))
nlm_rot1 = sf.rotate_nlm(nlm, theta, 0)    # first rotate around y axis in xz plane
nlm_rot2 = sf.rotate_nlm(nlm_rot1, 0, phi) # next  rotate around z axis in xy plane 

# See "plotting" pages on how to plot the resulting ODFs



def plot_vec(ax, v, lbl, color, ls='-', lw=2):
    
    ax.plot([0, +v[0]],[0, +v[1]],[0,+v[2]], color=color, ls=ls, lw=lw, label=lbl)
    ax.plot([0, -v[0]],[0, -v[1]],[0,-v[2]], color=color, ls=ls, lw=lw)


def print_ms(m1, m2, m3):
    y_offset, z_offset = angle_offsets(m1, m2, m3)
    print("ANGLES,",np.rad2deg(y_offset), np.rad2deg(z_offset))
    





def angle_offsets(m1, m2, m3):
    '''
    first return is the y rotation offset
    second return is the z rotation offset
    '''
    #y = m1 
    #z = m2
    y = m2 
    z = m3

    if (np.abs(z[0]) + np.abs(z[2]) == 0):
        y_r = 0
    else:
        y_r = (np.arctan(z[2]/np.sqrt((z[0]**2) + ( z[1]**2))) + np.deg2rad(90))


    if (np.abs(y[0]) + np.abs(y[2]) == 0):
        z_r = 0
    else:
        z_r = - ((np.arctan(y[1]/np.sqrt((y[0]**2) + ( y[2]**2)))) + np.deg2rad(90))

    return (y_r, z_r)


def arrange_m(m1,m2,m3):
    options = [m1,m2,m3]
    x_v = [m1[0], m2[0], m3[0]]
    y_v = [m1[1], m2[1], m3[1]]
    z_v = [m1[2], m2[2], m3[2]]

    x = -1
    y = -1
    z = -1

    x = x_v.index(np.max(x_v))
    y_v[x] = -10
    z_v[x] = -10

    y = y_v.index(np.max(y_v))
    z_v[y] = -10

    z = z_v.index(np.max(z_v))

    print(options[x], options[y], options[z])

    return (options[x], options[y], options[z]) 
    


def correct_angle(nlm):
    m1,m2,m3, eigvals = sf.frame(nlm, 'e')
    print(m1,m2,m3)

    #m1,m2,m3 = arrange_m(m1,m2,m3)
    y_offset, z_offset = angle_offsets(m1, m2, m3)
    return sf.rotate_nlm(nlm, y_offset, z_offset) # rotate around z axis in xy plane 






def plot_ms(ax, m1, m2, m3, rotation=55, inclination=45):
    #print(m1, m2, m3)

    ax.view_init(elev=90-inclination, azim=rotation) # same as ODF plot
    
    ax.set_xlabel('$x$'); ax.set_xlim([-1,1])
    ax.set_ylabel('$y$'); ax.set_ylim([-1,1])
    ax.set_zlabel('$z$'); ax.set_zlim([0,1])

    plot_vec(ax,[0,0,1], r'$\vb{c}_1$', sfplt.c_lred)
    plot_vec(ax,[0,1,0], r'$\vb{c}_2$', sfplt.c_lgreen)
    plot_vec(ax,[1,0,0], r'$\vb{c}_3$', sfplt.c_lblue)
    
    plot_vec(ax,m1[:], r'$\vb{m}_1$', sfplt.c_red)
    plot_vec(ax,m2[:], r'$\vb{m}_2$', sfplt.c_green)
    plot_vec(ax,m3[:], r'$\vb{m}_3$', sfplt.c_blue)
    lwpq=1
    
    ax.legend(handlelength=1, bbox_to_anchor=(1.17,1), fancybox=False)#, loc=1)





def rot_test():

    dpi, scale = 200, 3.3

    fig = plt.figure(figsize=(4*scale,scale))
    gs = gridspec.GridSpec(1,4, height_ratios=[1.5], width_ratios=[1,1,1,1])
    #gs.update(left=, right=1-0.06/3, top=0.97, bottom=0.20, wspace=0.015*18, hspace=0.35)
    ax_mi      = plt.subplot(gs[0, 0], projection='3d')
    ax_y = plt.subplot(gs[0, 1], projection='3d')
    ax_z = plt.subplot(gs[0, 2], projection='3d')
    ax_m_true      = plt.subplot(gs[0, 3], projection='3d')


    m1,m2,m3, eigvals = sf.frame(nlm_rot2, 'e')

    print_ms(m1, m2, m3)

    plot_ms(ax_mi, m1, m2, m3)
    ax_mi.set_title(f'Rotated Fabric Orientation')


    plot_ms(ax_z, m1, m2, m3, -90, 0)
    ax_y.set_title(f'{theta_d}° y-axis Rotation')
    plot_ms(ax_y, m1, m2, m3, -90, 90)
    ax_z.set_title(f'{phi_d}° z-axis Rotation')


    nlm_rot3 = correct_angle(nlm_rot2) # next  rotate around z axis in xy plane 
    m1,m2,m3, eigvals = sf.frame(nlm_rot3, 'e')
    plot_ms(ax_m_true, m1, m2, m3)
    ax_m_true.set_title(f'Corrected Rotation')

    fig.savefig('rotation_test.png', dpi=dpi)

rot_test()