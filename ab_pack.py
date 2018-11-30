#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 14:37:07 2017
Routines related to spherical harmonic representations of data.
The routines here are designed to work with Thorsten Becker's .ab file format
which is used with the HC mantle circulation code. The spherical harmonic representations
use the Dahlen and Tromp normalization convention.

A test suite is included in ab_pack_tests.py.


@author: max
"""

import numpy as np
import scipy
import scipy.optimize
import os
import subprocess
from io import StringIO
import copy
import cartopy.crs as ccrs
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
sh_syn_cmd  = './sh_syn'

def empty_layer():
    return {'lmax':None, 'll':[],'mm':[],'clm':[],'slm':[],'depth':[]}
    
class ab_field:
    def __init__(self,file=None,xyz=None,lmax=None):
        self.layers = []
        self.layers.append( empty_layer() )
        self.nlayer = 0
        self.layer_depths = []
        if file != None:
            self.read(file)
        elif xyz != None:
            if lmax == None:
                raise ValueError('lmax must be defined for the spherical harmonic expansion')
            self.xyz2plm(xyz,lmax)
                
    def format_ab(self,ilayer):
        layer = self.layers[ilayer]
        """convert to Thorsten Becker's ab-file format as a string."""
        ab_string = str(layer['lmax']) + " 0 0 1 1 0\n"
        ilm=0
        l=0
        while l <= layer['lmax']:
            m=0
            while m <= l:
                ab_string += ('{:e}'.format(layer['clm'][ilm]) + ' ' + '{:e}'.format(layer['slm'][ilm]) + "\n")
                ilm+=1
                m+=1
            l+=1
        return ab_string
    def write(self,output_ab_file):
        file=open(output_ab_file,'w')
        file.write(self.format_ab())
        file.close()        
    def read(self,input_ab_file):
        file = open(input_ab_file,'r')
        # read the header line
        ilayer=0
        ilm=-1
        l=0
        m=0
   
        for line in file:
            # 1st line is header
            if ilm == -1:
                header = line
                fields = header.split()
                layer=self.layers[ilayer]
                layer['lmax']   = int(fields[0])
                layer['ilayer'] = int(fields[1])
                layer['depth']  = float(fields[2])
                self.nlayer = int(fields[3])
                self.layer_depths.append(layer['depth'])
                ilm += 1 
            else:
                tmp = line.split()
                layer['clm'].append(float(tmp[0]))
                layer['slm'].append(float(tmp[1]))
                layer['ll'].append(l)
                layer['mm'].append(m)
                if( m<l ):
                    m+=1
                else:
                    l+=1
                    m=0
                ilm+=1
                if l>layer['lmax']:
                    # move on to the next layer
                    if ilayer == self.nlayer-1:
                        # This allows poorly-made single-layer ab-files to be truncated at the specified lmax.
                        file.close()
                        return
                    self.layers.append( empty_layer() )
                    ilm=-1
                    l=0
                    m=0
                    ilayer += 1
                    
        file.close()
    def filt(self,lkeep=range(0,360)):
        for layer in self.layers:
            l=0
            m=0
            ilm=0
            while l<= layer['lmax']:
                m=0
                while m <= l:
                    if l not in lkeep:
                        layer['clm'][ilm]=0.0
                        layer['slm'][ilm]=0.0
                        ilm+=1
                    m+=1
                l+=1
        return self
    def print(self):
        for ilayer in range(len(self.layers)):
            layer = self.layers[ilayer]
            l=0
            m=0
            ilm=0
            while l<=layer['lmax']:
                m=0
                while m<=l:
                    print(layer['ll'][ilm],layer['mm'][ilm],layer['clm'][ilm],layer['slm'][ilm])
                    ilm+=1
                    m+=1
                l+=1
    def to_xyz(self,method='scipy',spacing=1.0,layer_number=0):
        """Expand the one slice of the field out to a regular lat/lon grid using Thorsten Becker's sh_syn command"""
        if method == 'sh_syn':
            if spacing==1.0:
                cmd = 'echo \"' + self.format_ab(layer_number) + '\" | ' + sh_syn_cmd + ' 0 0 0.0 360.0 -90.0 90.0 1.0 '
                nlon=361
                nlat=181
            elif spacing==2.0:
                cmd = 'echo \"' + self.format_ab(layer_number) + '\" | ' + sh_syn_cmd + ' 0 0 0.5 359.5 -89.5 89.5 2.0 '
                nlon=180 # must be consistent with line above.
                nlat=90
            else:
                raise ValueError
            # call sh_syn
            result = subprocess.run( cmd ,shell=True,stdout=subprocess.PIPE )
            f = StringIO(result.stdout.decode())
            # load the xyz output
            x,y,z = np.loadtxt(f,usecols=(0,1,2),unpack=True)
            shp = (nlat,nlon)
            x = x.reshape(shp)
            y = y.reshape(shp)
            z = z.reshape(shp)
            return x,y,z
        elif method=='scipy':
            # use scipy to compute the spherical harmonic expansion
            if spacing==1.0:
                lon = np.linspace(0.0,360.0,361)
                lat = np.linspace(90.0,-90.0,181)
            elif spacing ==2.0: # cell-centered output
                lon = np.linspace(0.5,359.5,360)
                lat = np.linspace(-89.5,89.5,180)
            [lon,lat] = np.meshgrid(lon,lat)
            theta = lon*np.pi/180.   # scipy convention - azimuthal angle
            phi = (90.-lat)*np.pi/180. # scipy convention - polar angle
            z = np.zeros_like(lon)
            # evaluate each spherical harmonic
            layer = self.layers[0]
            for ilm in range(len(layer['ll'])):
                m = layer['mm'][ilm]
                assert(m>=0)
                l = layer['ll'][ilm]
                clm = layer['clm'][ilm]
                slm = layer['slm'][ilm]
                # generate real spherical harmonics from complex harmonics inmplemented in scipy.
                if( m==0):
                    if not np.isclose(clm,0.0):
                        # note - should be real anyways for m=0
                        z += clm*np.real( scipy.special.sph_harm(m,l,theta,phi) )
                else:
                    # don't bother evaluating anything unless clm or slm is nonzero
                    if not (np.isclose(clm,0.0) and np.isclose(slm,0.0) ):
                        ylm = scipy.special.sph_harm(m,l,theta,phi)
                        z +=  clm *np.sqrt(2.0)*np.real( ylm )
                        z +=  slm *np.sqrt(2.0)*np.imag( ylm )
            return lon,lat,z
    def xyz2plm(self,xyz,lmax,layer_number=0):
        xx=xyz[0].flatten()
        yy=xyz[1].flatten()
        zz=xyz[2].flatten()
        # calculate spherical harmonic coefficients for input longitude, latitude, z values
        theta = xx*np.pi/180.0
        phi = (90.-yy)*np.pi/180.0
        nlm=0
        for i in range(lmax+1):
            nlm += 2*i+1

        nx = xx.size
        # calculate spherical harmonic weights at each input point
        weights = np.zeros((nx,nlm))

        ilm=0
        ll = []
        mm = []
        for l in range(0,lmax+1):
            for m in range(0,l+1):
                if m==0:
                    ylm = np.real( scipy.special.sph_harm(m,l,theta,phi) )
                    weights[:,ilm] = ylm
                    ll.append(l)
                    mm.append(m)
                    ilm+=1
                else:                    
                    ylm = scipy.special.sph_harm(m,l,theta,phi)
                    weights[:,ilm] = np.sqrt(2.0)*np.real( ylm )
                    ll.append(l)
                    mm.append(-m)
                    ilm += 1
                    weights[:,ilm] = np.sqrt(2.0)*np.imag( ylm )
                    ll.append(l)
                    mm.append(m)
                    ilm+=1
        layer = self.layers[layer_number]
        # compute the coefficients using least squares
        coefs = scipy.optimize.lsq_linear(weights,zz)
        print(coefs.message) # This can't hurt!
        # retrieve the coefficients from the inversion
        ilm=0
        for l in range(0,lmax+1):
            for m in range(0,l+1):
                if m==0:
                      layer['ll'].append(l)
                      layer['mm'].append(m)
                      layer['clm'].append(coefs.x[ilm])
                      layer['slm'].append(0.0)
                      ilm += 1
                else:
                      layer['ll'].append(l)
                      layer['mm'].append(m)
                      layer['clm'].append(coefs.x[ilm])
                      ilm += 1
                      layer['slm'].append(coefs.x[ilm])
                      ilm += 1
        
    
    def __sub__(self,y):
        lmax_min = min(self.layers[0]['lmax'],y.layers[0]['lmax'])
        assert( len(self.layers) == len(y.layers) )
        result = copy.deepcopy(self)
        for ilayer in range(len(self.layers)):
            l=0
            m=0
            ilm=0
            while l<=lmax_min:
                m=0
                while m <= l:
                    result.clm[ilm] -= y.clm[ilm]
                    result.slm[ilm] -= y.slm[ilm]
                    ilm+=1
                    m+=1
                l+=1
        return result
    
    import matplotlib.pyplot as plt


def plot_geoid(field,vmin=None,vmax=None,layer_number=0):
    plot_ab(field,vmin,vmax,label='Geoid height (m)',layer_number=layer_number)
    
def plot_ab(field,vmin=None,vmax=None,label='',flipcb=False,file='',layer_number=0):
    x,y,z = field.to_xyz(layer_number=layer_number)
    fig = plt.figure()
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=0.0),position=(0.05,0.05,0.8,0.8))
    cax1 = fig.add_axes([0, 0, 0.1, 0.1])

    #lons, lats, data = sample_data()
    if vmin == None:
        vmin = -np.max(np.abs(z))
    if vmax == None:
        vmax = np.max(np.abs(z))
    if flipcb:
        cmap = cm.coolwarm_r
    else:
        cmap = cm.coolwarm
    
    geoid = ax.contourf(x, y, z,128,
                transform=ccrs.PlateCarree(),
                cmap=cmap,
                vmin=vmin,
                vmax=vmax)
    ax.coastlines()
    ax.set_global()

    cb = plt.colorbar(geoid,cax=cax1,cmap=geoid.cmap,label=label)
    plt.draw()
    posn = ax.get_position()
    cax1.set_position([posn.x0 + posn.width + 0.01, posn.y0,
                          0.04, posn.height])
    
    #plt.colorbar()
    if file != '':
        plt.savefig(file,boundingbox_inches='tight')
    plt.show()
    