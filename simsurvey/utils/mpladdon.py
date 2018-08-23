#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" This module gather the add-on methods for matplotlib axes and figure"""

import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.transforms  import Bbox
import warnings
from .tools import kwargs_update
from .decorators import make_method

__all__ = ["specplot","skyplot","figout"]

# ========================== #
# =  Axes Add-on           = #
# ========================== #
# --------------------- #
# - Matrix           - #
# --------------------- #
@make_method(mpl.Axes)
def corrmatrix(ax, matrix, npoints=None,
               cmap=None, significance_range=[0,6], ncolors=6,
               matrixlabels=None,
               **kwargs):
    """
    In Dev
    """
    # --------------
    # Init Settings
    
    shape = np.shape(matrix)
    def get_significance(rho, npoints):
        """ """
        rho = np.abs(rho)
        return rho * np.sqrt( (npoints-2.) /(1.-rho**2) )

    prop = dict(interpolation="nearest", origin="upper",
                alpha=0.9)
    
    if npoints is None:
        color_matrix = matrix
        cmap = mpl.cm.Spectral if cmap is None else \
          mpl.cm.get_cmap(cmap) if type(cmap) == str else cmap

    else:
        color_matrix = get_significance(matrix, npoints)
        cmap = mpl.cm.get_cmap("YlGnBu_r",ncolors) if cmap is None else cmap
        prop["vmin"],prop["vmax"] = significance_range
        cticks = np.linspace(prop["vmin"],prop["vmax"],ncolors+1)


    # - imshow
    propimshow = kwargs_update(prop,**kwargs)
    
    pl = mpl.imshow(color_matrix, cmap=cmap,**propimshow)

    # - ColorBar
    cb = mpl.colorbar(pl, ticks=cticks if npoints is not None else None)
    if npoints is not None:
        cb.ax.set_yticklabels([r'$%d$'%v if i<(len(cticks)-1) else r'$>%d$'%v for i,v in enumerate(cticks)],
                                fontsize="large")
        cb.set_label(r"$\mathrm{Corr.\ Significance\ [in\,\sigma]}$", fontsize="x-large")

    # --------
    # - Ticks
    ax.xaxis.set_ticks(np.linspace(0,shape[1]-1,shape[1]))
    ax.yaxis.set_ticks(np.linspace(0,shape[0]-1,shape[0]))


    
    # - Text
    
    _ = [[ax.text(i,j,"%.2f"%matrix[i][j],
            va="center", ha="center",
            fontproperties="black",color="k")
            for i in range(shape[0])] for j in range(shape[1])]
        
    if matrixlabels is not None:
        _ = ax.set_xticklabels(matrixlabels,
                    fontsize="xx-large")
        _ = ax.set_yticklabels(matrixlabels,
                                fontsize="xx-large")

    # - No ticks 
    [ax_.tick_params(length=0) for ax_ in [ax,cb.ax]]
    return pl

    
    



# --------------------- #
# - Scatter           - #
# --------------------- #
@make_method(mpl.Axes)
def errorscatter(ax,x,y,dx=None,dy=None,**kwargs):
    """
    In Dev
    """
    if dx is None and dy is None:
        return
    prop = kwargs_update({"ls":"None","marker":None,"zorder":2,
                          "ecolor":"0.7"},**kwargs)
    return ax.errorbar(x,y,xerr=dx,yerr=dy,**prop)

# --------------------- #
# - Spectrum          - #
# --------------------- #
@make_method(mpl.Axes)
def specplot(ax,x,y,var=None,
             color=None,bandprop={},
             err_onzero=False,**kwargs):
    """This function in a build-in axes method that enable to quickly and
    easily plot a spectrum.
    """
    # -----------------------
    # - Properties of plot
    default_kwargs = dict(
        color=mpl.cm.Blues(0.8),
        ls="-",lw=1,marker=None,zorder=6,
        )
    if color is not None:
        default_kwargs["color"] = color
    propplot = kwargs_update(default_kwargs,**kwargs)
    # -- Plot 
    pl = ax.plot(x,y,**propplot)
    
    # -----------------------
    # - Properties of band
    if var is not None:
        default_band   = dict(
            color=propplot["color"],alpha=0.3,
            zorder=3,label="_no_legend_"
            )
        bandprop = kwargs_update(default_band,**bandprop)
        # -- Band
        if not err_onzero:
            fill = ax.fill_between(x,y+np.sqrt(var),y-np.sqrt(var),
                            **bandprop)
        else:
            fill = ax.fill_between(x,np.sqrt(var),-np.sqrt(var),
                            **bandprop)
    else:
        fill = None
        
    return pl,fill

# --------------------- #
# - Skyplot           - #
# --------------------- #
@make_method(mpl.Axes)
def skyplot(ax, ra, dec, color=None, **kwargs):
    """
    Build-in axes method for easy plotting of points on the sky.

    ax: [matplotlib.axes]          where the data will be displayed
                                   should be using Mollweide or Hammer
                                   projection (not strictly enforced)

    ra, dec: [N array, N array]    arrays of sky coordinates (RA, Dec)
                                   in degrees

    - options -

    color: [string]                color to be used for the plot

    - kwargs goes to ax.plot -

    Return
    ------
    pl (output of ax.plot)
    """
    from .skyplot import convert_radec_azel
    # -----------------------
    # - Properties of plot
    default_kwargs = dict(marker='o', markersize=5, linestyle='none')
    if color is not None:
        default_kwargs["color"] = color
    propplot = kwargs_update(default_kwargs,**kwargs)

    az, el = convert_radec_azel(ra, dec)
    # -- Plot 
    pl = ax.plot(az, el, **propplot)
    
    return pl


@make_method(mpl.Axes)
def skyscatter(ax, ra, dec, **kwargs):
    """
    Build-in axes method for scatter plots of points on the sky.

    ax: [matplotlib.axes]          where the data will be displayed
                                   should be using Mollweide or Hammer
                                   projection (not strictly enforced)

    ra, dec: [N array, N array]    arrays of sky coordinates (RA, Dec)
                                   in degrees

    - kwargs goes to ax.scatter -

    Return
    ------
    sc (output of ax.scatter)
    """
    from .skyplot import convert_radec_azel
    # -----------------------
    # - Properties of plot
    default_kwargs = dict(marker='o', s=30,
                          cmap=mpl.cm.RdYlBu_r)

    propplot = kwargs_update(default_kwargs,**kwargs)

    az, el = convert_radec_azel(ra, dec)
    # -- Plot 
    sc = ax.scatter(az, el, **propplot)
    
    return sc

@make_method(mpl.Axes)
def skyhist(ax, ra=None, dec=None, values=None, bins=None, steps=None, max_stepsize=5, edge=1e-6,
            vmin=None, vmax=None, cmap=mpl.cm.Blues, cblabel=None, **kwargs):
    """This function in a build-in axes method that makes a sky histogram of the 
    coordinates.
    
    Build-in axes method for 2d-histograms of points on the sky.

    ax: [matplotlib.axes]          where the data will be displayed
                                   should be using Mollweide or Hammer
                                   projection (not strictly enforced)

    - options -

    ra, dec: [N array, N array]    arrays of sky coordinates (RA, Dec)
                                   in degrees

    values: [N array]              array of values for each bin (can be used 
                                   instead of ra and dec; must number of bins as
                                   length)
    
    bins [Bins object]             object that bins the coordinates and
                                   provides the boundaries for drawing
                                   polygons 
                                   (see the documentation of 
                                   utils/plot/skybins.py)

    steps [int]                    number of steps between bin verices
                                   (if None, it will be determined based on
                                   max_stepsize)
    
    max_stepsize [float]           maximal stepsize used to determined
                                   steps if None
    
    edges [float]                  edge to be left near RA = 180 deg

    vmin,vmax: [float,float]       minimal and maximal values for the colormapping.

    cmap: [mpl.cm]                 a colormap
    
    cblabel: [string]              colorbar label

    - kwargs goes to matplotlib.collections.PolyCollection -

    Return
    ------
    collection (output of matplotlib.collections.PolyCollection)
    cbar       (output of ax.colorbar)
    """
    # -----------------------
    # - Properties of plot
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection, PolyCollection

    from .skybins import SkyBins
    from .skyplot import convert_radec_azel
    
    if bins is None:
        bins = SkyBins()

    if cblabel is None:
        cblabel = ''

    if values is None:
        values = bins.hist(ra, dec)
    else:
        values = np.array(values)

    patches = []
    p_idx = []
    for k in range(len(values)):
        radec_bd = bins.boundary(k, steps=steps, max_stepsize=max_stepsize,
                                 edge=edge)
        for r, d in radec_bd:
            coord_bd = np.asarray(convert_radec_azel(r, d, edge=edge)).T
            patches.append(coord_bd)
            p_idx.append(k)

    c = np.asarray(values[np.asarray(p_idx)])
    vmin = c.min() if vmin is None else vmin
    vmax = c.max() if vmax is None else vmax
    color = cmap((c-vmin)/float(vmax-vmin))

    collec = PolyCollection(patches, facecolors=color, **kwargs)
    collec.set_edgecolor('face')
    
    # -- Plot 
    ax.add_collection(collec)

    axcar = ax.insert_ax('bottom',shrunk=0.85,space=0,axspace=0.08)

    cbar = axcar.colorbar(cmap, vmin=vmin, vmax=vmax, label=cblabel)

    return collec, cbar


# --------------------- #
# - WCS Plots         - #
# --------------------- #
@make_method(mpl.Axes)
def wcsplot(ax, wcs, exp_order=10,
            fc=mpl.cm.Blues(0.6,0.1),ec=mpl.cm.binary(0.8,1),
            draw_corner=False,
            **kwargs):
    """
    Comments TO BE DONE
    """
    from matplotlib.patches import Polygon
    # -----------------
    # - verticles
    if "has_contours" not in dir(wcs) or not wcs.has_contours():
        npoints = 2+exp_order
        width = np.linspace(0,wcs._naxis1,npoints)
        heigh = np.linspace(0,wcs._naxis2,npoints)
        v1 = np.asarray([np.ones(npoints-1)*0, width[:-1]]).T
        v2 = np.asarray([heigh[:-1], np.ones(npoints-1)*wcs._naxis1]).T
        v3 = np.asarray([np.ones(npoints-1)*wcs._naxis2, width[::-1][:-1]]).T
        v4 = np.asarray([heigh[::-1][:-1], np.ones(npoints-1)*0]).T
        v = np.asarray([wcs.pix2world(i,j)
                        for i,j in np.concatenate([v1,v2,v3,v4],axis=0)])
    else:
        from .shape import polygon_to_vertices
        v = polygon_to_vertices(wcs.contours)
        
    poly = Polygon(v,fc=fc,ec=ec,lw=1,**kwargs)
    # ------------------
    # - Draw
    # The point used
    pl = ax.plot(v.T[0],v.T[1],ls="None",marker="o",mfc=fc,mec=ec,
            visible=draw_corner)
    # The actual Patch
    ax.add_patch(poly)
    # ------------------
    # - Returns
    return pl, poly
    
# --------------------------- #
# -    Color Bar            - #
# --------------------------- #
@make_method(mpl.Axes)
def colorbar(ax,cmap,vmin=0,vmax=1,label="",
             fontsize="x-large",**kwargs):
    """ Set a colorbar in the given axis

    Parameters
    -----------
    ax: [mpl's Axes]
        Axis in which the colorbar will be drawn

    cmap: [mpl's colormap]
        A matplotlib colormap

    vmin, vmax: [float,float] -optional-
        Extend of the colormap, values of the upper and lower colors

    label, fontsize: [string, string/float] -optional-
        Label of the colorbar and its associated size
     
    **kwargs goes to matplotlib.colobar.ColorbarBase

    Return
    ------
    colorbar
    """
    import matplotlib
    if "orientation" not in kwargs.keys():
        bbox = ax.get_position()
        orientiation = "vertical" if bbox.xmax - bbox.xmin < bbox.ymax - bbox.ymin \
          else "horizontal"
        kwargs["orientation"] = orientiation

    norm    = matplotlib.colors.Normalize(vmin=vmin,vmax=vmax)
    c_bar   = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
                              norm=norm,**kwargs)
    
    c_bar.set_label(label,fontsize=fontsize)
    if "ticks" in kwargs.keys() and "ticklabels" not in kwargs.keys():
        c_bar.ax.set_xticklabels([r"%s"%v for v in kwargs["ticks"]])
        
    ax.set_xticks([]) if kwargs["orientation"] == "vertical" \
      else ax.set_yticks([])
    return c_bar


"""
@make_method(mpl.Axes)
def insert_ax(ax,space,pad,fraction=None,
              location="right"):
    ""
    Split the given axis to insert a new one at this location.
    This is typically useful to add colorbar for instance.

    Parameters
    ----------
    
    ax: [mpl.axes]                  The axis that will be split.

    space: [float]                  Extra space between the *ax* and the *newax* that
                                    will be created from the freed space.

    pad: [float]                    Extra space at the right/top of the *newax*

    fraction: [float /None]         The location where the axes will be
                                    split, given in faction of *ax* size.
                                    e.g., 0.5 means in the middle, while
                                    0.8 is at the far right / top (see *location*)
                                    If None: the new axis + given axis will fill exactly
                                    the same space as the given ax.
    - options -
    
    location: [string]              where you whish to add the *newax*. Only
                                    right and top are defined.

    Return
    ------
    axes (the *newax*, the input *ax* is reshape)
    ""
    known_location = ["right","top","left","bottom"]
    if location not in known_location:
        raise ValueError("'%s' is not a correct location"%location," These are: "+", ".join(known_location))

    if fraction is None:
        fraction = 1- (space + pad)
        
    if location == "right":
        bboxax,_,bboxnew,_ = ax.get_position().splitx(fraction,fraction+space,
                                                      fraction+space+pad)
    elif location == "left":
        bboxnew,_,bboxax,_ = ax.get_position().splitx(pad,pad+space,
                                                      fraction+space+pad)
    elif location == "top":
        bboxax,_,bboxnew,_ = ax.get_position().splity(fraction,fraction+space,
                                                      fraction+space+pad)
    else: # bottom
        bboxnew,_,bboxax,_ = ax.get_position().splity(pad,pad+space,
                                                      fraction+space+pad)

    ax.set_position(bboxax)
    return ax.figure.add_axes(bboxnew)
"""
    
# ========================== #
# =  axes manipulation     = #
# ========================== #
@make_method(mpl.Axes)
def insert_ax(ax,location,shrunk=0.7,space=.05,
              axspace=0.02,shareax=False,**kwargs):
    """ insert an axis at the requested location

              
    The new axis will share the main axis x-axis (location=top or bottom) or
    the y-axis (location=left or right).

    Parameters:
    -----------
    location: [string]
       top/bottom/left/right, i.e. where new axis will be set

    shrunk: [float]
        the main axis will be reduced by so much (0.7 = 70%).
        the new axis will take the room

    space: [float]
        extra space new axis does not use between it and the edge of
        the figure. (in figure unit, i.e., [0,1])

    axspace: [float]
        extra space new axis does not use between it and the input
        axis. (in figure unit, i.e., [0,1])

    shareax: [bool]
        The new axis will share the main axis x-axis (location=top or bottom) or
        the y-axis (location=left or right). If so, the axis ticks will be cleaned.
                           
    **kwargs goes to figure.add_axes() for the new axis

    Returns:
    --------
    axes (the new axis)
    """
    # --------------------
    # hist x
    # -------------------- #
    # -- keep trace of the original axes
    bboxorig = ax.get_position().frozen()

    if location in ["top","bottom"]:
        axhist = ax.figure.add_axes([0.1,0.2,0.3,0.4],sharex=ax if shareax else None,
                                    **kwargs) # This will be changed
        _bboxax = ax.get_position().shrunk(1,shrunk)
        _bboxhist = Bbox([[_bboxax.xmin, _bboxax.ymax+axspace ],
                          [_bboxax.xmax, bboxorig.ymax-space]])
        
        if location == "bottom":
            tanslate = _bboxhist.height + space+axspace
            _bboxhist = _bboxhist.translated(0, bboxorig.ymin-_bboxhist.ymin+space)
            _bboxax = _bboxax.translated(0,tanslate)
            
    # --------------------
    # hist y
    # -------------------- #            
    elif location in ["right","left"]:
        axhist = ax.figure.add_axes([0.5,0.1,0.2,0.42],sharey=ax if shareax else None,
                                    **kwargs) # This will be changed
        _bboxax = ax.get_position().shrunk(shrunk,1)
        _bboxhist = Bbox([[_bboxax.xmax+axspace, _bboxax.ymin ],
                          [bboxorig.xmax-space, _bboxax.ymax]])
        if location == "left":
            tanslate = _bboxhist.width + space + axspace
            _bboxhist = _bboxhist.translated(bboxorig.xmin-_bboxhist.xmin+space, 0)
            _bboxax = _bboxax.translated(tanslate,0)
        
    else:
        raise ValueError("location must be 'top'/'bottom'/'left' or 'right'")


    axhist.set_position(_bboxhist)
    ax.set_position(_bboxax)

    # ---------------------
    # remove their ticks
    if shareax:
        if location in ["top","right"]:
            [[label.set_visible(False) for label in lticks]
            for lticks in [axhist.get_xticklabels(),axhist.get_yticklabels()]]
        elif location == "bottom":
            [[label.set_visible(False) for label in lticks]
            for lticks in [ax.get_xticklabels(),axhist.get_yticklabels()]]
        elif location == "left":
            [[label.set_visible(False) for label in lticks]
            for lticks in [ax.get_yticklabels(),axhist.get_xticklabels()]]
    
    return axhist
        
# ========================== #
# =  Improved methods      = #
# ========================== #
@make_method(mpl.Axes)
def vline(ax,value,ymin=None,ymax=None,
                **kwargs):
     """ use this to help mpld3 """
     ymin,ymax = _read_bound_(ax.get_ylim(),ymin,ymax)
     ax.plot([value,value],[ymin,ymax],
             scalex=False,scaley=False,
             **kwargs)
     
@make_method(mpl.Axes)      
def hline(ax,value,xmin=None,xmax=None,
                **kwargs):
    """ use this to help mpld3 """
    xmin,xmax = _read_bound_(ax.get_xlim(),xmin,xmax)
    ax.plot([xmin,xmax],[value,value],
            scalex=False,scaley=False,
            **kwargs)

@make_method(mpl.Axes)
def hspan(ax,minvalue,maxvalue,
            xmin=None,xmax=None,
            **kwargs):
    """ use this to help mpld3 """
    lims = ax.get_xlim()
    xmin,xmax = _read_bound_(lims,xmin,xmax)
    ax.fill_betweenx([minvalue,maxvalue],xmax,x2=xmin,
                     **kwargs)
    ax.set_xlim(lims)

@make_method(mpl.Axes)
def vspan(ax,minvalue,maxvalue,ymin=None,ymax=None,
          **kwargs):
    """use this to help mpld3"""
    lims = ax.get_ylim()
    ymin,ymax = _read_bound_(lims,ymin,ymax)
    ax.fill_betweenx([ymin,ymax],[minvalue,minvalue],[maxvalue,maxvalue],
                     **kwargs)
    ax.set_ylim(lims)


# ========================== #
# =  Figure Add-on         = #
# ========================== #
@make_method(mpl.Figure)
def figout(fig,savefile=None,show=True,add_thumbnails=False,
           dpi=200):
    """This methods parse the show/savefile to know if the figure
    shall the shown or saved."""
    
    if savefile in ["dont_show","_dont_show_","_do_not_show_"]:
        show = False
        savefile = None

    if savefile is not None:
        fig.savefig(savefile+'.png',dpi=dpi)
        fig.savefig(savefile+'.pdf')
        if add_thumbnails:
            fig.savefig(savefile+"_thumb"+'.png',dpi=dpi/10.)
            
    elif show:
        fig.canvas.draw()
        fig.show()
        

@make_method(mpl.Figure)
def add_threeaxes(figure,xhist=True,yhist=True,
                  rect=[0.1,0.1,0.8,0.8],
                  shrunk=0.7, space=0, axspace=0.02,
                  **kwargs):
    """ Create an axis using the usual add_axes of matplotlib, but in addition
    add top and left axes if xhist and yhist are True, respectively
    **kwargs goes to mpl.Figure.add_axes()

    Parameters:
    -----------
    shrunk,space,axspace: [floats/2d-float]  the inserting parameters for the
                          histograms (see insert_ax). They could be 1d/float values
                          or 2d-arrays. In the first case, both axis will share the
                          same value, otherwise the first entry will be for x, the
                          second for y
                          
    """
    
    # ================
    # Input Parsing
    # ================ 
    if "__iter__" not in dir(shrunk):
        shrunk = [shrunk,shrunk]
    elif len(shrunk) == 1:
        shrunk = [shrunk[0],shrunk[0]]
    elif len(shrunk)>2:
        raise ValueError("shrunk cannot have more than 2 entries (x,y)")

    if "__iter__" not in dir(space):
        space = [space,space]
    elif len(space) == 1:
        space = [space[0],space[0]]
    elif len(space)>2:
        raise ValueError("space cannot have more than 2 entries (x,y)")
    
    if "__iter__" not in dir(axspace):
        axspace = [axspace,axspace]
    elif len(axspace) == 1:
        axspace = [axspace[0],axspace[0]]
    elif len(axspace)>2:
        raise ValueError("axspace cannot have more than 2 entries (x,y)")

    # ================
    # Axis Creation
    # ================ 
    ax = figure.add_axes(rect,**kwargs)
    # -- x axis
    if xhist:
        axhistx = ax.insert_ax("top", shrunk=shrunk[0],space=space[0],
                            axspace=axspace[0],shareax=True)
    else:
        axhistx = None
    # -- y axis
    if yhist:
        axhisty = ax.insert_ax("right", shrunk=shrunk[1],space=space[1],
                            axspace=axspace[1],shareax=True)
    else:
        axhisty = None
        
    if xhist and yhist:
        axhistx.set_position(axhistx.get_position().shrunk(shrunk[1],1))
        
    return ax,axhistx,axhisty

        

# ========================== #
# =  Tools                 = #
# ========================== #
def _read_bound_(lims,xmin,xmax):
    """
    """
    if xmin is None or xmax is None:
        size = lims[1] - lims[0]
    if xmin is None:
        xmin = lims[0] - size*3
    if xmax is None:
        xmax = lims[1] + size*3
    return xmin,xmax



def voronoi_grid(ax,xy):
    """
    TO BE REMOVED
    """
    def _adjust_bounds(ax, points):
        ptp_bound = points.ptp(axis=0)
        ax.set_xlim(points[:,0].min() - 0.1*ptp_bound[0],
                    points[:,0].max() + 0.1*ptp_bound[0])
        ax.set_ylim(points[:,1].min() - 0.1*ptp_bound[1],
                    points[:,1].max() + 0.1*ptp_bound[1])

    from scipy.spatial import Voronoi
    from matplotlib.patches import Polygon
    vor = Voronoi(xy)

    
    if vor.points.shape[1] != 2:
        raise ValueError("Voronoi diagram is not 2-D")

    ax.plot(vor.points[:,0], vor.points[:,1], '.')

    for simplex in vor.ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            ax.plot(vor.vertices[simplex,0], vor.vertices[simplex,1], 'k-')
                        
    ptp_bound = vor.points.ptp(axis=0)

    center = vor.points.mean(axis=0)
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.any(simplex < 0):
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex
            print(i,simplex, simplex[simplex >= 0][0])
            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[i] + direction * ptp_bound.max()

            ax.plot([vor.vertices[i,0], far_point[0]],
                    [vor.vertices[i,1], far_point[1]], 'k--')

    _adjust_bounds(ax, vor.points)
