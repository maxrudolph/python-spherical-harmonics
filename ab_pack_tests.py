#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 12:45:27 2018
This file contains routines designed to test my ab_pack.py
There are currently three tests.
1) For each spherical harmonic degree and angular order m, a set of spherical
harmonic coefficients is generated. An expansion is performed using Thorsten Becker's sh_syn
code as well as my own routine, which uses the spherical harmonics functions implemented 
in scipy. The results are compared visually and for numerical accuracy through an assertion.
2) Each of the spatial fields in (1) is represented in spherical harmonic form using
least-squares fitting. The spherical harmonic coefficients are compared with the input
coefficients. The result is then re-expanded on to the same spatial grid and 
checked for consistency with the input.
3) A reference geoid is read from a .ab file and plotted using the plot_geoid() function in ab_pack
Then, the reference geoid is evaluated on a regular lat/lon grid. Finally, we recompute spherical harmonic
coefficients and compare with the input through an assertion.

@author: Max Rudolph (maxrudolph@ucdavis.edu)
"""

from ab_pack import *
import matplotlib.pyplot as plt
from rem3d import plots
import cartopy.util as cu

# test the consistency of my expansions with Backer's sh_syn C-program
for lperturb in range(0,3):
    for mperturb in range(-lperturb,lperturb+1):
        # create a synthetic field
        field = ab_field()
        layer = field.layers[0]
        lmax=7
        layer['lmax'] = lmax
        ll=0
        mm=0
        layer['depth'] = 0
        while ll<=lmax:
            mm=0
            while mm<=ll:
                layer['ll'].append(ll)
                layer['mm'].append(mm)
                if ll==lperturb and mm==abs(mperturb) :
                    if mperturb<0 :
                        layer['slm'].append(1.0)                
                        layer['clm'].append(0.0)
                    else:
                        layer['slm'].append(0.0)                
                        layer['clm'].append(1.0)
                else:
                    layer['clm'].append(0.0)
                    layer['slm'].append(0.0)
                mm+=1
            ll+=1
    
        # expand the synthetic field using sh_syn
        axs=[]
        xx,yy,zz = field.to_xyz(layer_number=0,method='sh_syn')
        latlonval = np.vstack((xx.reshape(xx.size),yy.reshape(yy.size),zz.reshape(zz.size))).transpose()
        dt = {'names':['lon', 'lat', 'val'], 'formats':[np.float, np.float, np.float]}
        plotmodel = np.zeros(len(latlonval), dtype=dt);
        plotmodel['lon'] = latlonval[:,0]; 
        plotmodel['lat'] = latlonval[:,1]; 
        plotmodel['val'] = latlonval[:,2]*100.;
        vmin=-1.5; vmax=1.5
        plt.figure()
        ax = plt.subplot(1,3,1,projection=ccrs.Robinson(central_longitude=0.))
        axs.append(ax)    
        thecmap = plots.standardcolorpalette()
        plt.contourf(xx,yy,zz, 60,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap='RdYlBu') 
        plt.colorbar()
        ax.coastlines()
        
        # expand the synthetic field using my routine
        xx1,yy1,zz1 = field.to_xyz(layer_number=0,method='scipy')

        latlonval = np.vstack((xx1.reshape(xx1.size),yy1.reshape(yy1.size),zz1.reshape(zz1.size))).transpose()
        dt = {'names':['lon', 'lat', 'val'], 'formats':[np.float, np.float, np.float]}
        plotmodel = np.zeros(len(latlonval), dtype=dt);
        plotmodel['lon'] = latlonval[:,0]; 
        plotmodel['lat'] = latlonval[:,1]; 
        plotmodel['val'] = latlonval[:,2]*100.;
        vmin=-1.5; vmax=1.5
        ax = plt.subplot(1,3,2,projection=ccrs.Robinson(central_longitude=0.))
        axs.append(ax)    
        thecmap = plots.standardcolorpalette()
        plt.contourf(xx1,yy1,zz1, 60,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap='RdYlBu') 
        plt.colorbar()
        ax.coastlines()

        
        assert(np.isclose(zz,zz1).all())
        
        # test the spherical harmonic expansion routines.
        new_field = ab_field(xyz=[xx,yy,zz],lmax=7)
        # expand the synthetic field using my routine
        xx1,yy1,zz1 = new_field.to_xyz(layer_number=0,method='scipy')

        latlonval = np.vstack((xx1.reshape(xx1.size),yy1.reshape(yy1.size),zz1.reshape(zz1.size))).transpose()
        dt = {'names':['lon', 'lat', 'val'], 'formats':[np.float, np.float, np.float]}
        plotmodel = np.zeros(len(latlonval), dtype=dt);
        plotmodel['lon'] = latlonval[:,0]; 
        plotmodel['lat'] = latlonval[:,1]; 
        plotmodel['val'] = latlonval[:,2]*100.;
        vmin=-1.5; vmax=1.5
        ax = plt.subplot(1,3,3,projection=ccrs.Robinson(central_longitude=0.))
        axs.append(ax)    
        thecmap = plots.standardcolorpalette()
        plt.contourf(xx1,yy1,zz1, 60,transform=ccrs.PlateCarree(),vmin=vmin,vmax=vmax,cmap='RdYlBu') 
        plt.colorbar()
        ax.coastlines()
        plt.show()
        # assert that new expansion is 'close' to old values
        assert(np.isclose(new_field.layers[0]['clm'],field.layers[0]['clm']).all())
        assert(np.isclose(new_field.layers[0]['slm'],field.layers[0]['slm']).all())

geoid = ab_field(file='data/ggm05_chambat_geoid.ab')
# expand the geoid and plot it.
plot_geoid(geoid)
[xg,yg,zg] = geoid.to_xyz(layer_number=0,method='scipy')
geoid2 = ab_field(xyz=[xg,yg,zg],lmax=15)
np.isclose(geoid.layers[0]['clm'],geoid2.layers[0]['clm'])
np.isclose(geoid.layers[0]['slm'],geoid2.layers[0]['slm'])