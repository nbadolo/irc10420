#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 11:36:54 2022

@author: nbadolo
"""

"""
Code simplifié 
"""

import numpy as np
from astropy.io import fits
from scipy import optimize
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import math
import matplotlib.colors as colors
from matplotlib.pyplot import Figure, subplot
#%%
fdir='/home/nbadolo/Bureau/Aymard/Donnees_sph/'
fdir_star=fdir+'star/IRC10420/visible/data/'
fdir_psf=fdir+'psf/'
fname1='zpl_p23_make_polar_maps-ZPL_SCIENCE_P23_REDUCED'
fname2='-zpl_science_p23_REDUCED'
file_I_star= fdir_star +fname1+'_I'+fname2+'_I.fits'
file_PI_star= fdir_star+fname1+'_PI'+fname2+'_PI.fits'
file_DOLP_star= fdir_star+fname1+'_DOLP'+fname2+'_DOLP.fits'
file_AOLP_star= fdir_star+ fname1+'_AOLP'+fname2+'_AOLP.fits'

# file_I_psf= fdir_psf+ fname1+'_I'+fname2+'_I.fits'
# file_PI_psf= fdir_psf+fname1+'_PI'+fname2+'_PI.fits'
# file_DOLP_psf= fdir_psf+ fname1+'_DOLP'+fname2+'_DOLP.fits'
# file_AOLP_psf= fdir_psf+fname1+'_AOLP'+fname2+'_AOLP.fits'
#%%
file_lst = [file_I_star,file_PI_star,file_DOLP_star,file_AOLP_star]
file_lst_ = [file_I_star,file_PI_star,file_DOLP_star,file_AOLP_star]
nFrames = len(file_lst)
nFrames_ = len(file_lst_)
#%%
nDim=1024
nSubDim = 200 # plage de pixels que l'on veut afficher
nDimfigj=[9,10,11]
nDimfigk=[0,1,2]

vmin0=3.5
vmax0=15

x_min=-3.5*nSubDim//2
x_max=3.5*(nSubDim//2-1)
y_min=-3.5*nSubDim//2
y_max=3.5*(nSubDim//2-1)

x_min_=-3.5*nDim//2
x_max_=3.5*(nDim//2-1)
y_min_=-3.5*nDim//2
y_max_=3.5*(nDim//2-1)
#%%
#mean_sub_v_arr = np.empty((nFrames,nSubDim//2-1))
sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))
sub_v_arr_ = np.empty((nFrames,nDim,nDim))
im_name_lst = ['IRC+10420 I','IRC+10420 PI','IRC+10420 DOLP','IRC+10420 Vectors and PI',
               'IRC+10420 PI_', 'IRC+10420 Vectors and PI_']
Vmin = np.empty((nFrames))
Vmax = np.empty((nFrames))
Vmin_ = np.empty((nFrames))
Vmax_ = np.empty((nFrames_))
#%%
pix2mas=3.5  #en mas/pix
position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)
size_ = (nDim, nDim)

x, y = np.meshgrid(np.arange(nSubDim),np.arange(nSubDim)) #cree un tableau 
X, Y= np.meshgrid(np.linspace(-100,99,200), np.linspace(-100,99,200))
X_, Y_= np.meshgrid(np.linspace(-nDim/2,nDim/2-1,nDim), np.linspace(-nDim/2,nDim/2-1,nDim))
R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1)
r_mas=pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde
X *= pix2mas
Y *= pix2mas
X_ *= pix2mas
Y_ *= pix2mas
#%%
"""
### Filtre utilisé: I_PRIM 
"""
#%%
for i in range(nFrames):
    hdu = fits.open(file_lst[i])   
    data = hdu[0].data   
    i_v = data[0,:,:]
    print(np.shape(data))
    cutout = Cutout2D(i_v, position = position, size = size)
    zoom_hdu = hdu.copy()
    sub_v = cutout.data 
    sub_v_arr[i] = sub_v
    if i == 3:
        Vmin[i] = np.min(sub_v_arr[i])
        Vmax[i] = np.max(sub_v_arr[i]) 
    
    else:
        Vmin[i] = np.min(np.log10(sub_v_arr[i]))
        Vmax[i] = np.max(np.log10(sub_v_arr[i]))  
U = sub_v_arr[2]*np.cos(np.pi*sub_v_arr[3]/180)
V = sub_v_arr[2]*np.sin(np.pi*sub_v_arr[3]/180)
#%%
for i in range(nFrames_):
    hdu = fits.open(file_lst_[i])   
    data= hdu[0].data   
    i_v_= data[0,:,:]
    print(np.shape(data))
    sub_v_arr_[i] = i_v_
    
    if i == 1:
        Vmin_[i]=np.min(np.log10(sub_v_arr_[i]))
        Vmax_[i]=np.max(np.log10(sub_v_arr_[i])) 
    
    else:
        Vmin_[i]=np.min(sub_v_arr_[i])
        Vmax_[i]=np.max(sub_v_arr_[i])  
U_= sub_v_arr_[2]*np.cos(np.pi*sub_v_arr_[3]/180)
V_= sub_v_arr_[2]*np.sin(np.pi*sub_v_arr_[3]/180)

#%%
# [::X_step,::X_step]
X_step = 10
X_step_ = 50
plt.figure('I_PIM(visible)')
plt.clf()
for i in range (nFrames+2):   
    plt.subplot(3,2,i+1)
    if i==3:
       plt.imshow(np.log10(sub_v_arr[1]), cmap ='inferno', origin='lower',vmin=Vmin[1], 
                   vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])   
       plt.colorbar(label='ADU in log$_{10}$ scale')       
       q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U[::X_step,::X_step], V[::X_step,::X_step])
       plt.quiverkey(q, X = -0.05, Y = 1.03, U = 0.03, label='pol. degree vector norm scale 0.03 ', labelpos='E')
       plt.text(-17*size[0]//25., 3*size[1]//2, im_name_lst[3], color='cyan',
                fontsize='x-small', ha='center')
    else:
        if i == 4:
            plt.imshow(np.log10(sub_v_arr_[1]), cmap='inferno', origin='lower',
                       vmin=Vmin[1], vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])
            plt.colorbar(label='')
            plt.text(-17*size[0]//17., 3*size[1]//2, im_name_lst[4], color='cyan',
                     fontsize='x-small', ha='center')
        else:
            if i == 5:
                plt.imshow(np.log10(sub_v_arr_[1]), cmap='inferno', origin='lower',
                           vmin=Vmin[1], vmax=Vmax[1], extent = [x_min_ , x_max_, y_min_ , y_max_])
                plt.colorbar(label='ADU in log$_{10}$ scale')
                q_ = plt.quiver(X_[::X_step_,::X_step_], Y_[::X_step_,::X_step_], U_[::X_step_,::X_step_], V_[::X_step_,::X_step_])
                plt.quiverkey(q_, X = -0.05, Y = 1.03, U = 0.1, label='pol. degree vector norm scale 0.1 ', labelpos='E')
                plt.text(-17*size_[0]//25., 3*size_[1]//2, im_name_lst[5], color='cyan',
                         fontsize='x-small', ha='center')
            else:
                plt.imshow(np.log10(sub_v_arr[i]), cmap='inferno', origin='lower',
                           vmin=Vmin[i], vmax=Vmax[i], extent = [x_min , x_max, y_min , y_max])
                if i == 1:
                    plt.colorbar(label='ADU in log$_{10}$ scale')
                else:
                    plt.colorbar(label='') 
                plt.text(-17*size[0]//17., 3*size[1]//2, im_name_lst[i], color='cyan',
             fontsize='x-small', ha='center')
    #plt.colorbar(label='ADU in log$_{10}$ scale')
    #plt.clim(Vmin[i],Vmax[i])
    # if i == 0 or i == 1:
    #     plt.ylabel('Relative Dec.(mas)', size=10)
    # else:    
    #     plt.xlabel('Relative R.A.(mas)', size=10)   
    #     plt.ylabel('Relative Dec.(mas)', size=10)
    
    if i != 1 and i != 3:
        if i ==4:
            plt.xlabel('Relative R.A.(mas)', size=10) 
            plt.ylabel('Relative Dec.(mas)', size=10)
        else:
            if i == 0 or i == 2:
                plt.ylabel('Relative Dec.(mas)', size=10)
            else:
                plt.xlabel('Relative R.A.(mas)', size=10) 
#%%
"""
### Filtre utilisé: R_PRIM 
"""
#%%
for i in range(nFrames):
    hdu = fits.open(file_lst[i])   
    data = hdu[0].data   
    i_v = data[1,:,:]
    print(np.shape(data))
    cutout = Cutout2D(i_v, position = position, size = size)
    zoom_hdu = hdu.copy()
    sub_v = cutout.data 
    sub_v_arr[i] = sub_v
    if i == 3:
        Vmin[i] = np.min(sub_v_arr[i])
        Vmax[i] = np.max(sub_v_arr[i]) 
    
    else:
        Vmin[i] = np.min(np.log10(sub_v_arr[i]))
        Vmax[i] = np.max(np.log10(sub_v_arr[i]))  
U = sub_v_arr[2]*np.cos(np.pi*sub_v_arr[3]/180)
V = sub_v_arr[2]*np.sin(np.pi*sub_v_arr[3]/180)
#%%
for i in range(nFrames_):
    hdu = fits.open(file_lst_[i])   
    data= hdu[0].data   
    i_v_= data[1,:,:]
    print(np.shape(data))
    sub_v_arr_[i] = i_v_
    
    if i == 1:
        Vmin_[i]=np.min(np.log10(sub_v_arr_[i]))
        Vmax_[i]=np.max(np.log10(sub_v_arr_[i])) 
    
    else:
        Vmin_[i]=np.min(sub_v_arr_[i])
        Vmax_[i]=np.max(sub_v_arr_[i])  
U_= sub_v_arr_[2]*np.cos(np.pi*sub_v_arr_[3]/180)
V_= sub_v_arr_[2]*np.sin(np.pi*sub_v_arr_[3]/180)

#%%
# [::X_step,::X_step]
X_step = 10
X_step_ = 50
plt.figure('R_PIM(visible)')
plt.clf()
for i in range (nFrames+2):   
    plt.subplot(3,2,i+1)
    if i==3:
       plt.imshow(np.log10(sub_v_arr[1]), cmap ='inferno', origin='lower',vmin=Vmin[1], 
                   vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])   
       plt.colorbar(label='ADU in log$_{10}$ scale')       
       q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U[::X_step,::X_step], V[::X_step,::X_step])
       plt.quiverkey(q, X = -0.05, Y = 1.03, U = 0.03, label='pol. degree vector norm scale 0.03 ', labelpos='E')
       plt.text(-17*size[0]//25., 3*size[1]//2, im_name_lst[3], color='cyan',
                fontsize='x-small', ha='center')
    else:
        if i == 4:
            plt.imshow(np.log10(sub_v_arr_[1]), cmap='inferno', origin='lower',
                       vmin=Vmin[1], vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])
            plt.colorbar(label='')
            plt.text(-17*size[0]//17., 3*size[1]//2, im_name_lst[4], color='cyan',
                     fontsize='x-small', ha='center')
        else:
            if i == 5:
                plt.imshow(np.log10(sub_v_arr_[1]), cmap='inferno', origin='lower',
                           vmin=Vmin[1], vmax=Vmax[1], extent = [x_min_ , x_max_, y_min_ , y_max_])
                plt.colorbar(label='ADU in log$_{10}$ scale')
                q_ = plt.quiver(X_[::X_step_,::X_step_], Y_[::X_step_,::X_step_], U_[::X_step_,::X_step_], V_[::X_step_,::X_step_])
                plt.quiverkey(q_, X = -0.05, Y = 1.03, U = 0.1, label='pol. degree vector norm scale 0.1 ', labelpos='E')
                plt.text(-17*size_[0]//25., 3*size_[1]//2, im_name_lst[5], color='cyan',
                         fontsize='x-small', ha='center')
            else:
                plt.imshow(np.log10(sub_v_arr[i]), cmap='inferno', origin='lower',
                           vmin=Vmin[i], vmax=Vmax[i], extent = [x_min , x_max, y_min , y_max])
                if i == 1:
                    plt.colorbar(label='ADU in log$_{10}$ scale')
                else:
                    plt.colorbar(label='') 
                plt.text(-17*size[0]//17., 3*size[1]//2, im_name_lst[i], color='cyan',
             fontsize='x-small', ha='center')
    #plt.colorbar(label='ADU in log$_{10}$ scale')
    #plt.clim(Vmin[i],Vmax[i])
    # if i == 0 or i == 1:
    #     plt.ylabel('Relative Dec.(mas)', size=10)
    # else:    
    #     plt.xlabel('Relative R.A.(mas)', size=10)   
    #     plt.ylabel('Relative Dec.(mas)', size=10)
    
    if i != 1 and i != 3:
        if i ==4:
            plt.xlabel('Relative R.A.(mas)', size=10) 
            plt.ylabel('Relative Dec.(mas)', size=10)
        else:
            if i == 0 or i == 2:
                plt.ylabel('Relative Dec.(mas)', size=10)
            else:
                plt.xlabel('Relative R.A.(mas)', size=10) 



    
# for j in range(len(nDimfigj)):      
#     plt.subplot(3,4,nDimfigj[j])
#     plt.plot(r_mas, np.log10(mean_sub_v_arr[j]), color='darkorange',linewidth=2, label='Mira') 
#     plt.plot(r_mas, np.log10(mean_sub_v_arr[j+4]),color='purple',linewidth=2, label='HD204971') 
#     plt.legend(loc=0) 
#     plt.xlabel('r (mas)', size=10) 
#     if j==0:
#         plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)
#%%