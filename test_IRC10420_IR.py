#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 16:27:40 2022

@author: nbadolo
"""

"""
Code simplifié 
"""
import numpy as np
from astropy.io import fits
from scipy import optimize
import cv2  # pour ameliorer la résolution de l'image
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import math
import matplotlib.colors as colors
from matplotlib.pyplot import Figure, subplot
#%%
"""
### Calling files
"""

fdir= '/home/nbadolo/Bureau/Aymard/Donnees_sph/'
fdir_star = fdir+'star/IRC10420/IR/data/'
fdir_psf = fdir+'psf/'
fname1 ='IRC10420_2016-07-21_'
fname2 ='_star_pol_subtr'
file_I_star = fdir_star + fname1 +'I_tot.fits'
file_PI_star = fdir_star + fname1 +'I_pol' + fname2 + '.fits'
file_U_phi_star = fdir_star + fname1 +'U_phi' + fname2 + '.fits'
file_Q_star = fdir_star + fname1+'Q' + fname2 + '.fits'
file_U_star = fdir_star + fname1+'U' + fname2 + '.fits'
file_AOLP_star = fdir_star + fname1+'AoLP' + fname2 + '.fits'

# file_I_psf= fdir_psf+ fname1+'_I'+fname2+'_I.fits'
# file_PI_psf= fdir_psf+fname1+'_PI'+fname2+'_PI.fits'
# file_DOLP_psf= fdir_psf+ fname1+'_DOLP'+fname2+'_DOLP.fits'
# file_AOLP_psf= fdir_psf+fname1+'_AOLP'+fname2+'_AOLP.fits'
#print(cv2.__version__)
#%%
"""
### Parameters
"""
nDim = 1024
nSubDim = 400 # plage de pixels que l'on veut afficher
nDimfigj = [9,10,11]
nDimfigk = [0,1,2]

pix2mas = 12.255  #en mas/pix
x_min = -pix2mas*nSubDim//2
x_max = pix2mas*(nSubDim//2-1)
y_min = -pix2mas*nSubDim//2
y_max = pix2mas*(nSubDim//2-1)
position = (nDim//2,nDim//2)
size = (nSubDim, nSubDim)
Size = (pix2mas*nSubDim, pix2mas*nSubDim)
X_step = 10
my_dpi = 96
#%%
"""
### Creation of lists
"""
file_lst=[file_U_phi_star,file_PI_star,file_AOLP_star, file_AOLP_star]
file_lst2=[file_Q_star, file_U_star, file_I_star]
im_name_lst = ['IRC+10420 U_phi','IRC+10420 PI','IRC+10420 DoLP','IRC+10420 Vectors and PI']
nFrames = len(file_lst)
nFrames2 = len(file_lst2)
#%%
"""
### OpenCV bloc for super resolution
"""
# sr = cv2.dnn_superres.DnnSuperResImpl_create()
# path = "LapSRN_x8.pb"
# sr.readModel(path)
#%%
"""
### Creation of empty arrays
"""
mean_sub_v_arr = np.empty((nFrames,nSubDim//2-1))
sub_v_arr = np.empty((nFrames,nSubDim,nSubDim))
sub_v_arr2 = np.empty((nFrames2,nSubDim,nSubDim))
Vmin = np.empty((nFrames))
Vmax = np.empty((nFrames))

#%%
"""
### creation of value grids for plots 
"""
x, y = np.meshgrid(np.arange(nSubDim),np.arange(nSubDim)) #cree un tableau 
X, Y = np.meshgrid(np.linspace(-nSubDim/2,nSubDim/2-1,nSubDim), np.linspace(-nSubDim/2,nSubDim/2-1,nSubDim))
X*= pix2mas
Y*= pix2mas
R = np.sqrt((x-nSubDim/2)**2+(y-nSubDim/2)**2)
r = np.linspace(1,nSubDim//2-1,nSubDim//2-1)
r_mas = pix2mas*r #  où r est en pixels et r_mas en millièmes d'arcseconde
#%%
"""
### Opning and trimming of files for PI calculation
"""
for j in range(nFrames2):
    hdu2 = fits.open(file_lst2[j])   
    data2 = hdu2[0].data
    print(np.shape(data2))
    i_v2 = data2
       
    cutout = Cutout2D(i_v2, position=position, size=size)
    zoom_hdu = hdu2.copy()
    sub_v2 = cutout.data
    sub_v_arr2[j] = sub_v2

Q = sub_v_arr2[0]
W = sub_v_arr2[1]
I = sub_v_arr2[2]
    
#%%
"""
### Opning and trimming of files for plots
"""
for i in range(nFrames):
    hdu = fits.open(file_lst[i])   
    data = hdu[0].data 
    print(np.shape(data))
    i_v = data
   
    cutout = Cutout2D(i_v, position=position, size=size)
    zoom_hdu = hdu.copy()
    sub_v = cutout.data    
    #sub_v = sr.upsample(sub_v_lr) # for super resolution
    
    #sub_v = cv2.resize(sub_v_lr,dsize=None,fx=8,fy=8) # On pouvait utiliser cette 
    # ligne en lieu et place de la precedente pour pour améliorer la resoltion

    if i == 2 :
        sub_v_arr[i] = np.sqrt(Q**2 + W**2)/I
    else:
        sub_v_arr[i] = sub_v
    if i== 0 or i == 3:
        Vmin[i] = np.min(sub_v_arr[i])
        Vmax[i] = np.max(sub_v_arr[i]) 
    
    else:
        Vmin[i] = np.min(np.log10(sub_v_arr[i]))
        Vmax[i] = np.max(np.log10(sub_v_arr[i]))  
      

U = sub_v_arr[2]*np.cos(sub_v[3])
V = sub_v_arr[2]*np.sin(sub_v[3])   


          
#%%
"""
### All  plots
"""
plt.figure('irc10420_IR')
#plt.figure('irc10420_IR', figsize=(nSubDim/my_dpi, nSubDim/my_dpi), dpi=my_dpi)

plt.clf()
for i in range (nFrames):   
    plt.subplot(2,2,i+1)
    if i==3:        
       plt.imshow(np.log10(sub_v_arr[1]), cmap='inferno', origin='lower',vmin=Vmin[1], 
                   vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])
       plt.colorbar(label='ADU in log$_{10}$ scale')
       q = plt.quiver(X[::X_step,::X_step],Y[::X_step,::X_step],U[::X_step,::X_step], V[::X_step,::X_step])
       plt.quiverkey(q, X = 0.12, Y = 1.03, U = 0.03, label='pol. degree vector norm scale 0.03 ', labelpos='E')
       plt.text(-1.1*pix2mas*size[0]//6, 2*pix2mas*size[1]//6, im_name_lst[3], color='w',
                 fontsize='large', ha='center')
    else:
        if i == 0 :
            plt.imshow(sub_v_arr[i], cmap='inferno', origin='lower',
                       vmin=Vmin[i], vmax=Vmax[i], extent = [x_min , x_max, y_min , y_max])
            plt.colorbar(label='')
            plt.text(-1.5*pix2mas*size[0]//6., 2*pix2mas*size[1]//6, im_name_lst[i], color='w',
                fontsize='large', ha='center')
        else:           
            plt.imshow(np.log10(sub_v_arr[i]), cmap='inferno', origin='lower',
                       vmin=Vmin[i], vmax=Vmax[i], extent = [x_min , x_max, y_min , y_max])
            if i == 1:
                plt.colorbar(label='ADU in log$_{10}$ scale')
            else:
                plt.colorbar(label='') 
            plt.text(-1.5*pix2mas*size[0]//6., 2*pix2mas*size[1]//6, im_name_lst[i], color='w',
                     fontsize='large', ha='center')
    
    #plt.colorbar(label='ADU in log$_{10}$ scale')
    #plt.clim(Vmin[i],Vmax[i])
    if i != 1:
        if i ==2:
            plt.xlabel('Relative R.A.(mas)', size=10) 
            plt.ylabel('Relative Dec.(mas)', size=10)
        else:
            if i == 0:
                plt.ylabel('Relative Dec.(mas)', size=10)
            else:
                plt.xlabel('Relative R.A.(mas)', size=10) 

#for j in range(len(nDimfigj)):      
#     plt.subplot(3,4,nDimfigj[j])
#     plt.plot(r_mas, np.log10(mean_sub_v_arr[j]), color='darkorange',linewidth=2, label='Mira') 
#     plt.plot(r_mas, np.log10(mean_sub_v_arr[j+4]),color='purple',linewidth=2, label='HD204971') 
#     plt.legend(loc=0) 
#     plt.xlabel('r (mas)', size=10) 
#     if j==0:
#         plt.ylabel(r'Intensity in log$_{10}$ scale', size=10)

plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/star/IRC10420/IR/plots/irc10420.pdf', 
            dpi=100, bbox_inches='tight')
plt.tight_layout()
#%%

plt.figure('irc10420_IR_pol_i')
plt.clf()
plt.imshow(np.log10(sub_v_arr[1]), cmap='inferno', origin='lower',
           vmin=Vmin[1], vmax=Vmax[1], extent = [x_min , x_max, y_min , y_max])
plt.colorbar(label='ADU in log$_{10}$ scale')
plt.text(-1.5*pix2mas*size[0]//6., 2*pix2mas*size[1]//5., 'IRC+10420 Polarized Intensity', color='w',
         fontsize='xx-large', ha='center')
plt.savefig('/home/nbadolo/Bureau/Aymard/Donnees_sph/star/IRC10420/IR/plots/irc10420_IR_pol_i.pdf', 
            dpi=300, bbox_inches='tight')
plt.tight_layout()
#%%