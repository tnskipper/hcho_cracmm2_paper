#!/usr/bin/env python3

"""
Created 20211213
Updated 20230824 - updated several method and variable names
                   updated to use either netCDF4 or xarray Dataset

@author: Nash Skipper
"""

import numpy as np
import cartopy.crs as ccrs


class cmaqfile(object):
    
    """
    Interface to spatial data from CMAQ file metadata.
    Mostly for helping with making plots using cartopy.
    Note: Only Lambert Conformal Conic projection and North
          Polar Sterographic projection tested.

    Parameters
    ----------
    dataset : xarray or netCDF4 Dataset

    Returns
    -------
    ds : cmaqfile dataset object
        Original Dataset with new methods specific to CMAQ files.
    """
    
    
    def __init__(self, dataset):
        self.ds = dataset
    
    
    def getXYcenters(self):

        """
        Get X and Y grid cell centers (in m) from CMAQ file metadata.

        Returns
        -------
        X : ndarray of float
            X direction grid cell centers in m
        Y : ndarray of float
            Y direction grid cell centers in m

        """
        
        if self.ds.GDTYP==1:
            print('ERROR: Cannot use getXYcenters with lat-lon projection.')
            return
        
        xcell = self.ds.XCELL
        xstart = self.ds.XORIG
        xend = xstart + xcell*self.ds.NCOLS
        ycell = self.ds.YCELL
        ystart = self.ds.YORIG
        yend = ystart + ycell*self.ds.NROWS
        
        # grid cell centers
        X = np.arange(xstart+xcell/2., xend, xcell)
        Y = np.arange(ystart+ycell/2., yend, ycell)
        
        return X, Y
    
    
    def getXYcorners(self):

        """
        Get X and Y grid cell corners (in m) from CMAQ file metadata.

        Returns
        -------
        Xcorners : ndarray of float
            X direction grid cell corners in m
        Ycorners : ndarray of float
            Y direction grid cell corners in m
    
        """
        
        if self.ds.GDTYP==1:
            print('ERROR: Cannot use getXYcorners with lat-lon projection.')
            return
        
        xcell = self.ds.XCELL
        xstart = self.ds.XORIG
        xend = xstart + xcell*self.ds.NCOLS
        ycell = self.ds.YCELL
        ystart = self.ds.YORIG
        yend = ystart + ycell*self.ds.NROWS
        
        # grid cell corners
        X = np.arange(xstart, xend+xcell, xcell)
        Y = np.arange(ystart, yend+ycell, ycell)
        Xcorners, Ycorners = np.meshgrid(X, Y)
        
        return Xcorners, Ycorners
    
    
    def getCMAQproj(self, radius=6370000.):

        """
        Generate cartopy map projection from CMAQ file metadata.
        Note: Only Lambert Conformal Conic projection and North
              Polar Sterographic projection implemented.
        
        Parameters
        ----------
        radius : float, optional
            Assumed radius (in m) of the earth. The default is 6370000.
        
        Returns
        -------
        projection : cartopy crs
            Cartopy coordinate reference system generated from CMAQ file metadata.
    
        """
        
        # Lambert conformal conic
        if self.ds.GDTYP==2:
            centlon = self.ds.XCENT
            centlat = self.ds.YCENT
            stdpar = (self.ds.P_ALP, self.ds.P_BET)
            cmaqglobe = ccrs.Globe(
                ellipse=None,
                semimajor_axis=radius,
                semiminor_axis=radius
            )
            projection = ccrs.LambertConformal(
                central_longitude=centlon,
                central_latitude=centlat,
                standard_parallels=stdpar,
                globe=cmaqglobe
            )
        # north polar stereographic
        elif self.ds.GDTYP==6:
            centlon = self.ds.XCENT
            lat_true_scale = self.ds.P_BET
            cmaqglobe = ccrs.Globe(
                ellipse=None,
                semimajor_axis=radius,
                semiminor_axis=radius
            )
            projection = ccrs.NorthPolarStereo(
                central_longitude=centlon,
                true_scale_latitude=lat_true_scale,
                globe=cmaqglobe
            )
        # other projections not implemented
        else:
            print(
                f'''
                Projection GDTYP {self.ds.GDTYP} not supported.
                Only Lambert Conformal Conic (GDTYP=2) and 
                North Polar Stereogrpahic (GDTYP=6) projections
                have been implemented.
                '''
            )
            return
        
        return projection
    
    
    def ll2xy(self, lons, lats):

        """
        Get the X, Y coordinates of longitude, latitude points.

        Parameters
        ----------
        lons : float OR iterable of float
            Longitudes to be transformed.
        lats : float OR iterable of float
            Latitudes to be transformed.

        Returns
        -------
        xpts : list of float
            Coordinates of west-east dimension (X).
        ypts : list of float
            Coordinates of south-north dimension (Y).

        """
        
        if self.ds.GDTYP==1:
            print('ERROR: GDTYPE=1. Cannot use ll2xy with lat-lon projection.')
            return
        
        lons = np.asarray(lons)
        lats = np.asarray(lats)
        if lons.size != lats.size:
            print('ERROR: Number of longitude and latitude points must be equal.')
            return
        
        cmaqproj = self.getCMAQproj()
        xpts, ypts, _ = cmaqproj.transform_points(ccrs.PlateCarree(), lons, lats).T
        
        return xpts, ypts
    
    
    def xy2ll(self, X, Y):
        
        """
        Get the longitude, latitude coordinates of X, Y points.

        Parameters
        ----------
        X : float OR iterable of float
            Longitudes to be transformed.
        Y : float OR iterable of float
            Latitudes to be transformed.

        Returns
        -------
        lonpts : list of float
            Longitude coordinates.
        latpts : list of float
            Latitude coordinates.

        """
        
        if self.ds.GDTYP==1:
            print('ERROR: Cannot use ll2xy with lat-lon projection.')
            return
        
        X = np.asarray(X)
        Y = np.asarray(Y)
        if X.size != Y.size:
            print('ERROR: Number of X and Y points must be equal.')
            return
        
        cmaqproj = self.getCMAQproj()
        lonpts, latpts, _ = ccrs.PlateCarree().transform_points(cmaqproj, X, Y).T
        
        return lonpts, latpts
    
    
    def ll2ij(self, lons, lats):

        """
        Get the indices of longitude, latitude points.

        Parameters
        ----------
        lons : float OR iterable of float
            Longitudes to be transformed.
        lats : float OR iterable of float
            Latitudes to be transformed.

        Returns
        -------
        i : list of int
            Indices of west-east dimension. NaN if outside the domain.
        j : list of int
            Indices of south-north dimension. NaN if outside the domain.

        """
        
        if self.ds.GDTYP==1:
            print('ERROR: Cannot use ll2ij with lat-lon projection.')
            return
        
        lons = np.asarray(lons)
        lats = np.asarray(lats)
        if lons.size != lats.size:
            print('ERROR: Number of longitude and latitude points must be equal.')
            return
        
        Xcenters, Ycenters = self.getXYcenters()
        Xcorners, Ycorners = self.getXYcorners()
        xpts, ypts = self.ll2xy(lons, lats)
        
        i=[]
        j=[]
        for x, y in zip(xpts, ypts):
            if x < Xcorners.min() or x > Xcorners.max(): # outside the domain
                i.append(np.nan)
                j.append(np.nan)
            elif y < Ycorners.min() or y > Ycorners.max(): # outside the domain
                i.append(np.nan)
                j.append(np.nan)
            else:
                i.append(np.argmin(np.abs(Xcenters-x)))
                j.append(np.argmin(np.abs(Ycenters-y)))
        return i, j
