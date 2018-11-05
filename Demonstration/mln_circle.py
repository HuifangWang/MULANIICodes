# -*- coding: utf-8 -*-
"""
Created on Tue May 26 09:30:35 2015
adaptive from mne codes by authors belows
@author: Huifang Wang
Functions to plot on circle as for connectivity
"""
from __future__ import print_function

# Authors: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Denis Engemann <denis.engemann@gmail.com>
#          Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: Simplified BSD


from itertools import cycle
from functools import partial
import matplotlib.pyplot as plt

import numpy as np
import sys, os
from mne.externals.six import string_types
from mne.fixes import tril_indices, normalize_colors
sys.path.append( os.path.expanduser('/Users/huifangwang/MULANIII/Codes/BasicTools/') )
import mln_tools as mt

import matplotlib.colors as mcolors
def make_colormap(seq):
    """Return a LinearSegmentedColormap
        seq: a sequence of floats and RGB-tuples. The floats should be increasing
        and in the interval (0,1).
        """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

#def mln_cmap_mulan():
c = mcolors.ColorConverter().to_rgb
rvb_mulan = make_colormap(
                          [c('Hotpink'),0.5,c('green')] ) # sometime we need to
plt.register_cmap(name='mulan', cmap=rvb_mulan)

#rvb_vmulan = make_colormap(
#                   [c('Pink'),0.2,c('Hotpink'),0.3,c('DeepPink'),0.4,c('Green'),0.7,c('darkgreen')] )
rvb_vmulan = make_colormap(
                                              [c('hotpink'),0.2,c('deeppink'),0.35,c('magenta'),0.5,c('mediumseagreen'),0.7,c('darkgreen')] )
                           
plt.register_cmap(name='vmulan', cmap=rvb_vmulan)

def mln_degree(con,flagcon,threshold=0.1):
    if flagcon == 'Tplot':
        aij=con>threshold
    else:
        aij=con
    indgree=np.sum(aij,axis=1)
    outdgree=np.sum(aij,axis=0)
   
    return indgree, outdgree

def circular_layout(node_names, node_order, start_pos=90, start_between=True,
                    group_boundaries=None, group_sep=10):
    """Create layout arranging nodes on a circle.

    Parameters
    ----------
    node_names : list of str
        Node names.
    node_order : list of str
        List with node names defining the order in which the nodes are
        arranged. Must have the elements as node_names but the order can be
        different. The nodes are arranged clockwise starting at "start_pos"
        degrees.
    start_pos : float
        Angle in degrees that defines where the first node is plotted.
    start_between : bool
        If True, the layout starts with the position between the nodes. This is
        the same as adding "180. / len(node_names)" to start_pos.
    group_boundaries : None | array-like
        List of of boundaries between groups at which point a "group_sep" will
        be inserted. E.g. "[0, len(node_names) / 2]" will create two groups.
    group_sep : float
        Group separation angle in degrees. See "group_boundaries".

    Returns
    -------
    node_angles : array, shape=(len(node_names,))
        Node angles in degrees.
    """
    n_nodes = len(node_names)

    if len(node_order) != n_nodes:
        raise ValueError('node_order has to be the same length as node_names')

    if group_boundaries is not None:
        boundaries = np.array(group_boundaries, dtype=np.int)
        if np.any(boundaries >= n_nodes) or np.any(boundaries < 0):
            raise ValueError('"group_boundaries" has to be between 0 and '
                             'n_nodes - 1.')
        if len(boundaries) > 1 and np.any(np.diff(boundaries) <= 0):
            raise ValueError('"group_boundaries" must have non-decreasing '
                             'values.')
        n_group_sep = len(group_boundaries)
    else:
        n_group_sep = 0
        boundaries = None

    # convert it to a list with indices
    node_order = [node_order.index(name) for name in node_names]
    node_order = np.array(node_order)
    if len(np.unique(node_order)) != n_nodes:
        raise ValueError('node_order has repeated entries')

    node_sep = (360. - n_group_sep * group_sep) / n_nodes

    if start_between:
        start_pos += node_sep / 2

        if boundaries is not None and boundaries[0] == 0:
            # special case when a group separator is at the start
            start_pos += group_sep / 2
            boundaries = boundaries[1:] if n_group_sep > 1 else None

    node_angles = np.ones(n_nodes, dtype=np.float) * node_sep
    node_angles[0] = start_pos
    if boundaries is not None:
        node_angles[boundaries] += group_sep

    node_angles = np.cumsum(node_angles)[node_order]

    return node_angles


def _plot_connectivity_circle_onpick(event, fig=None, axes=None, indices=None,
                                     n_nodes=0, node_angles=None, ylim=[9, 10]):
    """Isolates connections around a single node when user left clicks a node.

    On right click, resets all connections."""
    if event.inaxes != axes:
        return

    if event.button == 1:  # left click
        # click must be near node radius
        if not ylim[0] <= event.ydata <= ylim[1]:
            return

        # all angles in range [0, 2*pi]
        node_angles = node_angles % (np.pi * 2)
        node = np.argmin(np.abs(event.xdata - node_angles))

        patches = event.inaxes.patches
        for ii, (x, y) in enumerate(zip(indices[0], indices[1])):
            patches[ii].set_visible(node in [x, y])
        fig.canvas.draw()
    elif event.button == 3:  # right click
        patches = event.inaxes.patches
        for ii in xrange(np.size(indices, axis=1)):
            patches[ii].set_visible(True)
        fig.canvas.draw()



def plot_connectivity_circle_dir(con, node_names, indices=None, n_lines=None,
                             node_angles=None, node_width=None,
                             node_colors=None, facecolor='black',
                             textcolor='white', node_edgecolor='black',barhight=1.,
                             linewidth=1.5, levels=['con','din','dout'],
                             flagcon='Tplot',threshold=0.1,
                             bottom_lev1=9, colormap='hot', vmin=None,
                             levelmap=['None','deeppink','orangered'],
                             vmax=None, colorbar=True, title=None,
                             colorbar_size=0.3, colorbar_pos=(-0.3, 0.1),
                             fontsize_title=12, fontsize_names=8,
                             fontsize_colorbar=8, padding=6.,
                             fig=None, subplot=111, interactive=True):
    """Visualize connectivity as a directed circular graph.
        
        Note: This code is based on the circle graph example by Nicolas P. Rougier
        
        Parameters
        ----------
        con : array
        Connectivity scores. Can be a square matrix, or a 1D array. If a 1D
        array is provided, "indices" has to be used to define the connection
        indices.
        node_names : list of str
        Node names. The order corresponds to the order in con.
        indices : tuple of arrays | None
        Two arrays with indices of connections for which the connections
        strenghts are defined in con. Only needed if con is a 1D array.
        n_lines : int | None
        If not None, only the n_lines strongest connections (strength=abs(con))
        are drawn.
        node_angles : array, shape=(len(node_names,)) | None
        Array with node positions in degrees. If None, the nodes are equally
        spaced on the circle. See mne.viz.circular_layout.
        node_width : float | None
        Width of each node in degrees. If None, the minimum angle between any
        two nodes is used as the width.
        node_colors : list of tuples | list of str
        List with the color to use for each node. If fewer colors than nodes
        are provided, the colors will be repeated. Any color supported by
        matplotlib can be used, e.g., RGBA tuples, named colors.
        facecolor : str
        Color to use for background. See matplotlib.colors.
        textcolor : str
        Color to use for text. See matplotlib.colors.
        node_edgecolor : str
        Color to use for lines around nodes. See matplotlib.colors.
        linewidth : float
        Line width to use for connections.
        colormap : str
        Colormap to use for coloring the connections.
        vmin : float | None
        Minimum value for colormap. If None, it is determined automatically.
        vmax : float | None
        Maximum value for colormap. If None, it is determined automatically.
        colorbar : bool
        Display a colorbar or not.
        title : str
        The figure title.
        colorbar_size : float
        Size of the colorbar.
        colorbar_pos : 2-tuple
        Position of the colorbar.
        fontsize_title : int
        Font size to use for title.
        fontsize_names : int
        Font size to use for node names.
        fontsize_colorbar : int
        Font size to use for colorbar.
        padding : float
        Space to add around figure to accommodate long labels.
        fig : None | instance of matplotlib.pyplot.Figure
        The figure to use. If None, a new figure with the specified background
        color will be created.
        subplot : int | 3-tuple
        Location of the subplot when creating figures with multiple plots. E.g.
        121 or (1, 2, 1) for 1 row, 2 columns, plot 1. See
        matplotlib.pyplot.subplot.
        interactive : bool
        When enabled, left-click on a node to show only connections to that
        node. Right-click shows all connections.
        
        Returns
        -------
        fig : instance of matplotlib.pyplot.Figure
        The figure handle.
        axes : instance of matplotlib.axes.PolarAxesSubplot
        The subplot handle.
        """
    import matplotlib.pyplot as plt
    import matplotlib.path as m_path
    import matplotlib.patches as m_patches
    import matplotlib as mpl
    
    n_nodes = len(node_names)
    
    if node_angles is not None:
        if len(node_angles) != n_nodes:
            raise ValueError('node_angles has to be the same length '
                             'as node_names')
        # convert it to radians
        node_angles = node_angles * np.pi / 180
    else:
        # uniform layout on unit circle
        node_angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    
    if node_width is None:
        # widths correspond to the minimum angle between two nodes
        dist_mat = node_angles[None, :] - node_angles[:, None]
        dist_mat[np.diag_indices(n_nodes)] = 1e9
        node_width = np.min(np.abs(dist_mat))
    else:
        node_width = node_width * np.pi / 180
    
    if node_colors is not None:
        if len(node_colors) < n_nodes:
            node_colors = cycle(node_colors)
    else:
        # assign colors using colormap
        node_colors = [plt.cm.spectral(i / float(n_nodes))
                       for i in range(n_nodes)]
    
    # handle 1D and 2D connectivity information
    if con.ndim == 1:
        if indices is None:
            raise ValueError('indices has to be provided if con.ndim == 1')
    elif con.ndim == 2:
        if con.shape[0] != n_nodes or con.shape[1] != n_nodes:
            raise ValueError('con has to be 1D or a square matrix')
        # we use the lower-triangular part
        con = con
    else:
        raise ValueError('con has to be 1D or a square matrix')
    
    # get the colormap
    if isinstance(colormap, string_types):
        colormap = plt.get_cmap(colormap)

    # Make figure background the same colors as axes
    if fig is None:
        fig = plt.figure(figsize=(8, 8), facecolor=facecolor)
    
    # Use a polar axes
    if not isinstance(subplot, tuple):
        subplot = (subplot,)
    axes = plt.subplot(*subplot, polar=True, axisbg=facecolor)

# No ticks, we'll put our own
    plt.xticks([])
    plt.yticks([])
    
    # Set y axes limit, add additonal space if requested
    plt.ylim(0, 10 + padding)
    
    # Remove the black axes border which may obscure the labels
    axes.spines['polar'].set_visible(False)
    
    # Draw lines between connected nodes, only draw the strongest connections
    if n_lines is not None and len(con)*(len(con)-1) > n_lines:
        con_thresh = np.sort(np.abs(con).ravel())[-n_lines]
    else:
        con_thresh = threshold
    #print(n_lines,len(con),con_thresh)
    # get the connections which we are drawing and sort by connection strength
    # this will allow us to draw the strongest connections first

    con_abs = np.abs(con)
    con_draw_idx = np.where(con_abs >= con_thresh)
#print con_draw_idx, con_thresh,

    con = con[con_draw_idx]
    con_abs = con_abs[con_draw_idx]
#indices = [ind[con_draw_idx] for ind in indices]
#print con
# now sort them
    sort_idx = np.argsort(con_abs)
    con_abs = con_abs[sort_idx]
    con = con[sort_idx]
    indices = [ind[sort_idx] for ind in con_draw_idx]
    #print(indices,con)
    # Get vmin vmax for color scaling
    if vmin is None:
        vmin = np.min(con[np.abs(con) >= con_thresh])
    if vmax is None:
        vmax = np.max(con)
    vrange = vmax - vmin
    
    # We want to add some "noise" to the start and end position of the
    # edges: We modulate the noise with the number of connections of the
    # node and the connection strength, such that the strongest connections
    # are closer to the node center
    nodes_n_con = np.zeros((n_nodes), dtype=np.int)
    for i, j in zip(indices[0], indices[1]):
        nodes_n_con[i] += 1
        nodes_n_con[j] += 1

    # initalize random number generator so plot is reproducible
    rng = np.random.mtrand.RandomState(seed=0)

    n_con = len(indices[0])
    noise_max = 0.25 * node_width
    start_noise = rng.uniform(-noise_max, noise_max, n_con)
    end_noise = rng.uniform(-noise_max, noise_max, n_con)
    
    nodes_n_con_seen = np.zeros_like(nodes_n_con)
    for i, (start, end) in enumerate(zip(indices[0], indices[1])):
        nodes_n_con_seen[start] += 1
        nodes_n_con_seen[end] += 1
        
        start_noise[i] *= ((nodes_n_con[start] - nodes_n_con_seen[start])
                           / float(nodes_n_con[start]))
        end_noise[i] *= ((nodes_n_con[end] - nodes_n_con_seen[end])
                                            / float(nodes_n_con[end]))

# scale connectivity for colormap (vmin<=>0, vmax<=>1)
    con_val_scaled = (con - vmin) / vrange
    
    indgree = np.zeros(n_nodes)
    outdgree = np.zeros(n_nodes)
    # Finally, we draw the connections
    for pos, (i, j) in enumerate(zip(indices[0], indices[1])):
        indgree[i] = indgree[i]+1
        outdgree[j] = outdgree[j]+1
        
        # Start point
        t0, r0 = node_angles[j], 10
        
        # End point
        t1, r1 = node_angles[i], 10
        
        # Some noise in start and end point
        t0 += start_noise[pos]
        t1 += end_noise[pos]
        
        verts = [(t0, r0), (t0, 5), (t1, 5), (t1, r1)]
        codes = [m_path.Path.MOVETO, m_path.Path.CURVE4, m_path.Path.CURVE4,
                 m_path.Path.LINETO]
        path = m_path.Path(verts, codes)
                 
                 #print(verts,codes)
        #print(i,j)         
        color = colormap(con_val_scaled[pos])
                 #plt.plot(verts[3][0],9,marker='^',color=color)
        plt.arrow(t1,r1-1,0,0.1,alpha = 1,head_width=0.05,
                           facecolor = color,edgecolor = color,head_length=0.3)
                 # Actual line
        patch = m_patches.PathPatch(path, fill=False, edgecolor=color,
                                             linewidth=linewidth, alpha=1.)
                 #patch = m_patches.ConnectionPatch(path,arrowstyle=u'-')
                 
        axes.add_patch(patch)
    #axes.arrow(0,0,
    #          0.5,0.5, head_width=0.05, head_length=1, ec='k'
    #         )

    height = np.ones(n_nodes) * barhight

    for indl, ilevel in enumerate(levels):
        bottom = bottom_lev1+indl*1.4*barhight# Draw ring with colored nodes
        
        bars = axes.bar(node_angles, height, width=node_width, bottom=bottom,
                        edgecolor=node_edgecolor, lw=2, facecolor='.9',
                        align='center')
                        # Draw ring with different colored nodes
                        
        if ilevel == 'con':
            node_colors_lev=node_colors
        elif ilevel == 'din':
            #indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],indgree)

        elif ilevel == 'dout':
#indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],outdgree)

        for bar, color in zip(bars, node_colors_lev):
            bar.set_facecolor(color)

# Draw node labels
    angles_deg = 180 * node_angles / np.pi
    for name, angle_rad, angle_deg in zip(node_names, node_angles, angles_deg):
        if angle_deg >= 270:
            ha = 'left'
        else:
            # Flip the label, so text is always upright
            angle_deg += 180
            ha = 'right'
        
        axes.text(angle_rad, bottom+1.5, name, size=fontsize_names,
                  rotation=angle_deg, rotation_mode='anchor',
                  horizontalalignment=ha, verticalalignment='center',
                  color=textcolor)
    
    if title is not None:
        plt.title(title, color=textcolor, fontsize=fontsize_title,
                  axes=axes)

    if colorbar:
    # for links
          norm = normalize_colors(vmin=vmin, vmax=vmax)
          sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
          sm.set_array(np.linspace(vmin, vmax))
          cb = plt.colorbar(sm, ax=axes, use_gridspec=False,
                              shrink=colorbar_size,
                              anchor=colorbar_pos)
          cb_yticks = plt.getp(cb.ax.axes, 'yticklabels')
          cb.ax.tick_params(labelsize=fontsize_colorbar)
          plt.setp(cb_yticks, color=textcolor)
          
          # for indegree
          for indl,idegree in enumerate(['in','out']):
              if idegree == 'in':
                  degrees = indgree
                  deLabel = 'indegree'
              else:
                  degrees = outdgree
                  deLabel = 'outdegree'
                                      
              ax3 = fig.add_axes([0.05+indl*0.4, 0.1, 0.3, 0.02])
              listmap=mt.mln_generateColormap(levelmap[indl+1])
              cmap = mpl.colors.ListedColormap(listmap)
              bounds = np.linspace(np.min(degrees),np.max(degrees),10).astype(int)
                                                      #print(bounds)
              norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
              cb3 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,
                                              norm=norm,
                                              # boundaries=[-10]+bounds+[10],
                                              extend='both',
                                              # Make the length of each extension
                                              # the same as the length of the
                                              # interior colors:
                                              ticks=bounds,
                                              spacing='uniform',
                                              orientation='horizontal')
              cb3.set_label(deLabel)
    #Add callback for interaction
    if interactive:
        callback = partial(_plot_connectivity_circle_onpick, fig=fig,
                           axes=axes, indices=indices, n_nodes=n_nodes,
                           node_angles=node_angles)
                           
        fig.canvas.mpl_connect('button_press_event', callback)

    return fig, axes

def plot_connectivity_circle_dir_2para(con, node_names, indices=None, n_lines=None,threshold=0.1,
                             node_angles=None, node_width=None,
                             node_colors=None, facecolor='black',
                             textcolor='white', node_edgecolor='black',barhight=1.,
                             linewidth=1.5, levels=['con','din','dout'],
                             flagcon='Tplot',
                             bottom_lev1=9, colormap='hot', vmin=None,
                             levelmap=['None','deeppink','orangered'],
                             vmax=None, colorbar=True, title=None,
                             colorbar_size=0.3, colorbar_pos=(0.3, 0.1),
                             fontsize_title=12, fontsize_names=8,
                             fontsize_colorbar=8, padding=6.,
                             fig=None, subplot=121, interactive=True,color_given=0):
    """Visualize connectivity as a directed circular graph.
        
        Note: This code is based on the circle graph example by Nicolas P. Rougier
        
        Parameters
        ----------
        con : array
        Connectivity scores. Can be a square matrix, or a 1D array. If a 1D
        array is provided, "indices" has to be used to define the connection
        indices.
        node_names : list of str
        Node names. The order corresponds to the order in con.
        indices : tuple of arrays | None
        Two arrays with indices of connections for which the connections
        strenghts are defined in con. Only needed if con is a 1D array.
        n_lines : int | None
        If not None, only the n_lines strongest connections (strength=abs(con))
        are drawn.
        node_angles : array, shape=(len(node_names,)) | None
        Array with node positions in degrees. If None, the nodes are equally
        spaced on the circle. See mne.viz.circular_layout.
        node_width : float | None
        Width of each node in degrees. If None, the minimum angle between any
        two nodes is used as the width.
        node_colors : list of tuples | list of str
        List with the color to use for each node. If fewer colors than nodes
        are provided, the colors will be repeated. Any color supported by
        matplotlib can be used, e.g., RGBA tuples, named colors.
        facecolor : str
        Color to use for background. See matplotlib.colors.
        textcolor : str
        Color to use for text. See matplotlib.colors.
        node_edgecolor : str
        Color to use for lines around nodes. See matplotlib.colors.
        linewidth : float
        Line width to use for connections.
        colormap : str
        Colormap to use for coloring the connections.
        vmin : float | None
        Minimum value for colormap. If None, it is determined automatically.
        vmax : float | None
        Maximum value for colormap. If None, it is determined automatically.
        colorbar : bool
        Display a colorbar or not.
        title : str
        The figure title.
        colorbar_size : float
        Size of the colorbar.
        colorbar_pos : 2-tuple
        Position of the colorbar.
        fontsize_title : int
        Font size to use for title.
        fontsize_names : int
        Font size to use for node names.
        fontsize_colorbar : int
        Font size to use for colorbar.
        padding : float
        Space to add around figure to accommodate long labels.
        fig : None | instance of matplotlib.pyplot.Figure
        The figure to use. If None, a new figure with the specified background
        color will be created.
        subplot : int | 3-tuple
        Location of the subplot when creating figures with multiple plots. E.g.
        121 or (1, 2, 1) for 1 row, 2 columns, plot 1. See
        matplotlib.pyplot.subplot.
        interactive : bool
        When enabled, left-click on a node to show only connections to that
        node. Right-click shows all connections.
        
        Returns
        -------
        fig : instance of matplotlib.pyplot.Figure
        The figure handle.
        axes : instance of matplotlib.axes.PolarAxesSubplot
        The subplot handle.
        """
    import matplotlib.pyplot as plt
    import matplotlib.path as m_path
    import matplotlib.patches as m_patches
    import matplotlib as mpl
    
    n_nodes = len(node_names)
    con_cp=con
    if node_angles is not None:
        if len(node_angles) != n_nodes:
            raise ValueError('node_angles has to be the same length '
                             'as node_names')
        # convert it to radians
        node_angles = node_angles * np.pi / 180
    else:
        # uniform layout on unit circle
        node_angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    
    if node_width is None:
        # widths correspond to the minimum angle between two nodes
        dist_mat = node_angles[None, :] - node_angles[:, None]
        dist_mat[np.diag_indices(n_nodes)] = 1e9
        node_width = np.min(np.abs(dist_mat))
    else:
        node_width = node_width * np.pi / 180
    
    if node_colors is not None:
        if len(node_colors) < n_nodes:
            node_colors = cycle(node_colors)
    else:
        # assign colors using colormap
        node_colors = [plt.cm.spectral(i / float(n_nodes))
                       for i in range(n_nodes)]
    
    # handle 1D and 2D connectivity information
    if con.ndim == 1:
        if indices is None:
            raise ValueError('indices has to be provided if con.ndim == 1')
    elif con.ndim == 2:
        if con.shape[0] != n_nodes or con.shape[1] != n_nodes:
            raise ValueError('con has to be 1D or a square matrix')
        # we use the lower-triangular part
        con = con
    else:
        raise ValueError('con has to be 1D or a square matrix')
    
    # get the colormap
    if isinstance(colormap, string_types):
        colormap = plt.get_cmap(colormap)
    
    # Make figure background the same colors as axes
    if fig is None:
        fig = plt.figure(figsize=(8, 8), facecolor=facecolor)
    
    # Use a polar axes
    if not isinstance(subplot, tuple):
        subplot = (subplot,)
    axes = plt.subplot(*subplot, polar=True, axisbg=facecolor)

# No ticks, we'll put our own
    plt.xticks([])
    plt.yticks([])
    
    # Set y axes limit, add additonal space if requested
    plt.ylim(0, 10 + padding)
    
    # Remove the black axes border which may obscure the labels
    axes.spines['polar'].set_visible(False)
    
    # Draw lines between connected nodes, only draw the strongest connections
    if n_lines is not None and len(con)*(len(con)-1) > n_lines:
        con_thresh = np.sort(np.abs(con).ravel())[-n_lines]
    else:
        con_thresh = threshold
    #print(n_lines,len(con),con_thresh)
    # get the connections which we are drawing and sort by connection strength
    # this will allow us to draw the strongest connections first

    con_abs = np.abs(con)
    con_draw_idx = np.where(con_abs >= con_thresh)
#print con_draw_idx, con_thresh,

    con = con[con_draw_idx]
    con_abs = con_abs[con_draw_idx]
#indices = [ind[con_draw_idx] for ind in indices]
#print con
# now sort them
    sort_idx = np.argsort(con_abs)
    con_abs = con_abs[sort_idx]
    con = con[sort_idx]
    indices = [ind[sort_idx] for ind in con_draw_idx]
    #print(indices,con)
    # Get vmin vmax for color scaling
    if vmin is None:
        vmin = np.min(con[np.abs(con) >= con_thresh])
    if vmax is None:
        vmax = np.max(con)
    vrange = vmax - vmin
    
    # We want to add some "noise" to the start and end position of the
    # edges: We modulate the noise with the number of connections of the
    # node and the connection strength, such that the strongest connections
    # are closer to the node center
    nodes_n_con = np.zeros((n_nodes), dtype=np.int)
    for i, j in zip(indices[0], indices[1]):
        nodes_n_con[i] += 1
        nodes_n_con[j] += 1

    # initalize random number generator so plot is reproducible
    rng = np.random.mtrand.RandomState(seed=0)

    n_con = len(indices[0])
    noise_max = 0.25 * node_width
    start_noise = rng.uniform(-noise_max, noise_max, n_con)
    end_noise = rng.uniform(-noise_max, noise_max, n_con)
    
    nodes_n_con_seen = np.zeros_like(nodes_n_con)
    for i, (start, end) in enumerate(zip(indices[0], indices[1])):
        nodes_n_con_seen[start] += 1
        nodes_n_con_seen[end] += 1
        
        start_noise[i] *= ((nodes_n_con[start] - nodes_n_con_seen[start])
                           / float(nodes_n_con[start]))
        end_noise[i] *= ((nodes_n_con[end] - nodes_n_con_seen[end])
                                            / float(nodes_n_con[end]))

# scale connectivity for colormap (vmin<=>0, vmax<=>1)
    con_val_scaled = (con - vmin) / vrange
    
    indgree = np.zeros(n_nodes)
    outdgree = np.zeros(n_nodes)
    # Finally, we draw the connections
    for pos, (i, j) in enumerate(zip(indices[0], indices[1])):
        indgree[i] = indgree[i]+1
        outdgree[j] = outdgree[j]+1
        
        # Start point
        t0, r0 = node_angles[j], 10
        
        # End point
        t1, r1 = node_angles[i], 10
        
        # Some noise in start and end point
        t0 += start_noise[pos]
        t1 += end_noise[pos]
        
        verts = [(t0, r0), (t0, 5), (t1, 5), (t1, r1)]
        codes = [m_path.Path.MOVETO, m_path.Path.CURVE4, m_path.Path.CURVE4,
                 m_path.Path.LINETO]
        path = m_path.Path(verts, codes)
                 
                 #print(verts,codes)
        #print(i,j)         

        if color_given==0:
            color = colormap(con_val_scaled[pos])
            
        else:
            color = np.array([48,128,20,255])/255.
                 #plt.plot(verts[3][0],9,marker='^',color=color)
        plt.arrow(t1,r1-1,0,0.1,alpha = 1,head_width=0.05,
                           facecolor = color,edgecolor = color,head_length=0.3)
                 # Actual line
        patch = m_patches.PathPatch(path, fill=False, edgecolor=color,
                                             linewidth=linewidth, alpha=1.)
                 #patch = m_patches.ConnectionPatch(path,arrowstyle=u'-')
                 
        axes.add_patch(patch)
    #axes.arrow(0,0,
    #          0.5,0.5, head_width=0.05, head_length=1, ec='k'
    #         )
    

    height = np.ones(n_nodes) * barhight

    for indl, ilevel in enumerate(levels):
        bottom = bottom_lev1+indl*1.4*barhight# Draw ring with colored nodes
        
        bars = axes.bar(node_angles, height, width=node_width, bottom=bottom,
                        edgecolor=node_edgecolor, lw=2, facecolor='.9',
                        align='center')
                        # Draw ring with different colored nodes
                        
        if ilevel == 'con':
            node_colors_lev=node_colors
        elif ilevel == 'din':
            #indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],indgree)

        elif ilevel == 'dout':
#indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],outdgree)

        for bar, color in zip(bars, node_colors_lev):
            bar.set_facecolor(color)

# Draw node labels
    angles_deg = 180 * node_angles / np.pi
    for name, angle_rad, angle_deg in zip(node_names, node_angles, angles_deg):
        if angle_deg >= 270:
            ha = 'left'
        else:
            # Flip the label, so text is always upright
            angle_deg += 180
            ha = 'right'
        
        axes.text(angle_rad, bottom+1.5, name, size=fontsize_names,
                  rotation=angle_deg, rotation_mode='anchor',
                  horizontalalignment=ha, verticalalignment='center',
                  color=textcolor)
    
    if title is not None:
        plt.title(title, color=textcolor, fontsize=fontsize_title,
                  axes=axes)

    if colorbar:
    # for links
          norm = normalize_colors(vmin=vmin, vmax=vmax)
          sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
          sm.set_array(np.linspace(vmin, vmax))
          cb = plt.colorbar(sm, ax=axes, use_gridspec=False,
                              shrink=colorbar_size,
                              anchor=colorbar_pos)
          cb_yticks = plt.getp(cb.ax.axes, 'yticklabels')
          cb.ax.tick_params(labelsize=fontsize_colorbar)
          plt.setp(cb_yticks, color=textcolor)
          
          # for indegree
          for indl,idegree in enumerate(['in','out']):
              if idegree == 'in':
                  degrees = indgree
                  deLabel = 'indegree'
              else:
                  degrees = outdgree
                  deLabel = 'outdegree'
              pos_axes=axes.get_position()
              
              ax3 = fig.add_axes([pos_axes.x0, 0.17-indl*0.05, 0.27, 0.02])
              
              #ax3 = fig.add_axes([0.05+indl*0.4, 0.1, 0.3, 0.02])
              listmap=mt.mln_generateColormap(levelmap[indl+1])
              cmap = mpl.colors.ListedColormap(listmap)
              bounds = np.linspace(np.min(degrees),np.max(degrees),10).astype(int)
                                                      #print(bounds)
              norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
              cb3 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,
                                              norm=norm,
                                              # boundaries=[-10]+bounds+[10],
                                              extend='both',
                                              # Make the length of each extension
                                              # the same as the length of the
                                              # interior colors:
                                              ticks=bounds,
                                              spacing='uniform',
                                              orientation='horizontal')
                                              #cb3.set_label(deLabel)
              cb3.ax.tick_params(labelsize=fontsize_names)
              ax3.text(pos_axes.x0+0.01,0.17-indl*0.05,deLabel,size=fontsize_names)
    #Add callback for interaction
    if interactive:
        callback = partial(_plot_connectivity_circle_onpick, fig=fig,
                           axes=axes, indices=indices, n_nodes=n_nodes,
                           node_angles=node_angles)
                           
        fig.canvas.mpl_connect('button_press_event', callback)
    
    ax = fig.add_axes([0.4, 0.3, 0.5, 0.5])
    con_cp[con_cp<con_thresh]=0
    cim=plt.imshow(con_cp,cmap='Blues',interpolation='nearest')
    xytickindex=np.arange(0,n_nodes)
    ax.set_xticks(xytickindex)
    ax.set_yticks(xytickindex)
    ax.set_xticklabels(node_names,size=fontsize_names-2,rotation=90)
    ax.set_yticklabels(node_names,size=fontsize_names-2)
    cb = plt.colorbar(cim, ax=ax, use_gridspec=False,
                              shrink=colorbar_size*2,
                              anchor=colorbar_pos)
    
    return fig, axes




"""
    only output a cicle"""


def plot_connectivity_circle_dir_signal(con, node_names, indices=None, n_lines=None,
                                       node_angles=None, node_width=None,
                                       node_colors=None, facecolor='black',
                                       textcolor='white', node_edgecolor='black',barhight=1.,
                                       linewidth=1.5, levels=['con','din','dout'],
                                       flagcon='Tplot',threshold=0.1,
                                       bottom_lev1=9, colormap='hot', vmin=None,
                                       levelmap=['None','deeppink','orangered'],
                                       vmax=None, colorbar=True, title=None,
                                       colorbar_size=0.3, colorbar_pos=(-0.3, 0.1),
                                       fontsize_title=18, fontsize_names=8,
                                       fontsize_colorbar=8, padding=6.,
                                       fig=None, subplot=111,
                                       interactive=True,color_given=0):
    """Visualize connectivity as a directed circular graph.
        
        Note: This code is based on the circle graph example by Nicolas P. Rougier
        
        Parameters
        ----------
        con : array
        Connectivity scores. Can be a square matrix, or a 1D array. If a 1D
        array is provided, "indices" has to be used to define the connection
        indices.
        node_names : list of str
        Node names. The order corresponds to the order in con.
        indices : tuple of arrays | None
        Two arrays with indices of connections for which the connections
        strenghts are defined in con. Only needed if con is a 1D array.
        n_lines : int | None
        If not None, only the n_lines strongest connections (strength=abs(con))
        are drawn.
        node_angles : array, shape=(len(node_names,)) | None
        Array with node positions in degrees. If None, the nodes are equally
        spaced on the circle. See mne.viz.circular_layout.
        node_width : float | None
        Width of each node in degrees. If None, the minimum angle between any
        two nodes is used as the width.
        node_colors : list of tuples | list of str
        List with the color to use for each node. If fewer colors than nodes
        are provided, the colors will be repeated. Any color supported by
        matplotlib can be used, e.g., RGBA tuples, named colors.
        facecolor : str
        Color to use for background. See matplotlib.colors.
        textcolor : str
        Color to use for text. See matplotlib.colors.
        node_edgecolor : str
        Color to use for lines around nodes. See matplotlib.colors.
        linewidth : float
        Line width to use for connections.
        colormap : str
        Colormap to use for coloring the connections.
        vmin : float | None
        Minimum value for colormap. If None, it is determined automatically.
        vmax : float | None
        Maximum value for colormap. If None, it is determined automatically.
        colorbar : bool
        Display a colorbar or not.
        title : str
        The figure title.
        colorbar_size : float
        Size of the colorbar.
        colorbar_pos : 2-tuple
        Position of the colorbar.
        fontsize_title : int
        Font size to use for title.
        fontsize_names : int
        Font size to use for node names.
        fontsize_colorbar : int
        Font size to use for colorbar.
        padding : float
        Space to add around figure to accommodate long labels.
        fig : None | instance of matplotlib.pyplot.Figure
        The figure to use. If None, a new figure with the specified background
        color will be created.
        subplot : int | 3-tuple
        Location of the subplot when creating figures with multiple plots. E.g.
        121 or (1, 2, 1) for 1 row, 2 columns, plot 1. See
        matplotlib.pyplot.subplot.
        interactive : bool
        When enabled, left-click on a node to show only connections to that
        node. Right-click shows all connections.
        
        Returns
        -------
        fig : instance of matplotlib.pyplot.Figure
        The figure handle.
        axes : instance of matplotlib.axes.PolarAxesSubplot
        The subplot handle.
        """
    import matplotlib.pyplot as plt
    import matplotlib.path as m_path
    import matplotlib.patches as m_patches
    import matplotlib as mpl
    
    n_nodes = len(node_names)
    con_cp=con
    if node_angles is not None:
        if len(node_angles) != n_nodes:
            raise ValueError('node_angles has to be the same length '
                             'as node_names')
        # convert it to radians
        node_angles = node_angles * np.pi / 180
    else:
        # uniform layout on unit circle
        node_angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    
    if node_width is None:
        # widths correspond to the minimum angle between two nodes
        dist_mat = node_angles[None, :] - node_angles[:, None]
        dist_mat[np.diag_indices(n_nodes)] = 1e9
        node_width = np.min(np.abs(dist_mat))
    else:
        node_width = node_width * np.pi / 180
    
    if node_colors is not None:
        if len(node_colors) < n_nodes:
            node_colors = cycle(node_colors)
    else:
        # assign colors using colormap
        node_colors = [plt.cm.spectral(i / float(n_nodes))
                       for i in range(n_nodes)]
    
    # handle 1D and 2D connectivity information
    if con.ndim == 1:
        if indices is None:
            raise ValueError('indices has to be provided if con.ndim == 1')
    elif con.ndim == 2:
        if con.shape[0] != n_nodes or con.shape[1] != n_nodes:
            raise ValueError('con has to be 1D or a square matrix')
        # we use the lower-triangular part
        con = con
    else:
        raise ValueError('con has to be 1D or a square matrix')
    
    # get the colormap
    if isinstance(colormap, string_types):
        colormap = plt.get_cmap(colormap)
    

    
    # Make figure background the same colors as axes
    if fig is None:
        fig = plt.figure(figsize=(8, 8), facecolor=facecolor)
    
    # Use a polar axes
    if not isinstance(subplot, tuple):
        subplot = (subplot,)
    axes = plt.subplot(*subplot, polar=True, axisbg=facecolor)

# No ticks, we'll put our own
    plt.xticks([])
    plt.yticks([])
    
    # Set y axes limit, add additonal space if requested
    plt.ylim(0, 10 + padding)
    
    # Remove the black axes border which may obscure the labels
    axes.spines['polar'].set_visible(False)
    
    # Draw lines between connected nodes, only draw the strongest connections
    if n_lines is not None and len(con)*(len(con)-1) > n_lines:
        con_thresh = np.sort(np.abs(con).ravel())[-n_lines]
    else:
        con_thresh = threshold
    #print(n_lines,len(con),con_thresh)
    # get the connections which we are drawing and sort by connection strength
    # this will allow us to draw the strongest connections first

    con_abs = np.abs(con)
    con_draw_idx = np.where(con_abs >= con_thresh)
#print con_draw_idx, con_thresh,

    con = con[con_draw_idx]
    con_abs = con_abs[con_draw_idx]
#indices = [ind[con_draw_idx] for ind in indices]
#print con
# now sort them
    sort_idx = np.argsort(con_abs)
    con_abs = con_abs[sort_idx]
    con = con[sort_idx]
    indices = [ind[sort_idx] for ind in con_draw_idx]
    #print(indices,con)
    # Get vmin vmax for color scaling
    if vmin is None:
        vmin = np.min(con[np.abs(con) >= con_thresh])
    if vmax is None:
        vmax = np.max(con)
    vrange = vmax - vmin
    
    # We want to add some "noise" to the start and end position of the
    # edges: We modulate the noise with the number of connections of the
    # node and the connection strength, such that the strongest connections
    # are closer to the node center
    nodes_n_con = np.zeros((n_nodes), dtype=np.int)
    for i, j in zip(indices[0], indices[1]):
        nodes_n_con[i] += 1
        nodes_n_con[j] += 1

    # initalize random number generator so plot is reproducible
    rng = np.random.mtrand.RandomState(seed=0)

    n_con = len(indices[0])
    noise_max = 0.25 * node_width
    start_noise = rng.uniform(-noise_max, noise_max, n_con)
    end_noise = rng.uniform(-noise_max, noise_max, n_con)
    
    nodes_n_con_seen = np.zeros_like(nodes_n_con)
    for i, (start, end) in enumerate(zip(indices[0], indices[1])):
        nodes_n_con_seen[start] += 1
        nodes_n_con_seen[end] += 1
        
        start_noise[i] *= ((nodes_n_con[start] - nodes_n_con_seen[start])
                           / float(nodes_n_con[start]))
        end_noise[i] *= ((nodes_n_con[end] - nodes_n_con_seen[end])
                                            / float(nodes_n_con[end]))

# scale connectivity for colormap (vmin<=>0, vmax<=>1)
    con_val_scaled = (con - vmin) / vrange
    
    indgree = np.zeros(n_nodes)
    outdgree = np.zeros(n_nodes)
    # Finally, we draw the connections
    for pos, (i, j) in enumerate(zip(indices[0], indices[1])):
        indgree[i] = indgree[i]+1
        outdgree[j] = outdgree[j]+1
        
        # Start point
        t0, r0 = node_angles[j], 10
        
        # End point
        t1, r1 = node_angles[i], 10
        
        # Some noise in start and end point
        t0 += start_noise[pos]
        t1 += end_noise[pos]
        
        verts = [(t0, r0), (t0, 5), (t1, 5), (t1, r1)]
        codes = [m_path.Path.MOVETO, m_path.Path.CURVE4, m_path.Path.CURVE4,
                 m_path.Path.LINETO]
        path = m_path.Path(verts, codes)
                 
                 #print(verts,codes)
        #print(i,j)         

        if color_given==0:
            color = colormap(con_val_scaled[pos])
        
        else:
#color = np.array([48,128,20,255])/255. % green
#color = np.array([173,255,47,255])/255.
            color = np.array([34,139,34,255])/255.

                 #plt.plot(verts[3][0],9,marker='^',color=color)
        plt.arrow(t1,r1-1,0,0.1,alpha = 1,head_width=0.05,
                           facecolor = color,edgecolor = color,head_length=0.3)
                 # Actual line
        patch = m_patches.PathPatch(path, fill=False, edgecolor=color,
                                             linewidth=linewidth, alpha=1.)
                 #patch = m_patches.ConnectionPatch(path,arrowstyle=u'-')
                 
        axes.add_patch(patch)
    #axes.arrow(0,0,
    #          0.5,0.5, head_width=0.05, head_length=1, ec='k'
    #         )
    

    height = np.ones(n_nodes) * barhight

    for indl, ilevel in enumerate(levels):
        bottom = bottom_lev1+indl*1.4*barhight# Draw ring with colored nodes
        
        bars = axes.bar(node_angles, height, width=node_width, bottom=bottom,
                        edgecolor=node_edgecolor, lw=2, facecolor='.9',
                        align='center')
                        # Draw ring with different colored nodes
                        
        if ilevel == 'con':
            node_colors_lev=node_colors
        elif ilevel == 'din':
            #indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],indgree)

        elif ilevel == 'dout':
#indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],outdgree)

        for bar, color in zip(bars, node_colors_lev):
            bar.set_facecolor(color)

# Draw node labels
    angles_deg = 180 * node_angles / np.pi
    for name, angle_rad, angle_deg in zip(node_names, node_angles, angles_deg):
        if angle_deg >= 270:
            ha = 'left'
        else:
            # Flip the label, so text is always upright
            angle_deg += 180
            ha = 'right'
        
        axes.text(angle_rad, bottom+1.5, name, size=fontsize_names,
                  rotation=angle_deg, rotation_mode='anchor',
                  horizontalalignment=ha, verticalalignment='center',
                  color=textcolor)
    
    if title is not None:
        plt.title(title, color=textcolor, fontsize=fontsize_title,
                  axes=axes)

    if colorbar:
    # for links
          norm = normalize_colors(vmin=vmin, vmax=vmax)
          sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
          sm.set_array(np.linspace(vmin, vmax))
          cb = plt.colorbar(sm, ax=axes, use_gridspec=False,
                              shrink=colorbar_size,
                              anchor=colorbar_pos)
          cb_yticks = plt.getp(cb.ax.axes, 'yticklabels')
          cb.ax.tick_params(labelsize=fontsize_colorbar)
          plt.setp(cb_yticks, color=textcolor)
          
          # for indegree
          for indl,idegree in enumerate(['in','out']):
              if idegree == 'in':
                  degrees = indgree
                  deLabel = 'indegree'
              else:
                  degrees = outdgree
                  deLabel = 'outdegree'
              pos_axes=axes.get_position()
                                      
              ax3 = fig.add_axes([pos_axes.x0, 0.18-indl*0.05, 0.27, 0.02])
              listmap=mt.mln_generateColormap(levelmap[indl+1])
              cmap = mpl.colors.ListedColormap(listmap)
              bounds = np.linspace(np.min(degrees),np.max(degrees),10).astype(int)
                                                      #print(bounds)
              norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
              cb3 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,
                                              norm=norm,
                                              # boundaries=[-10]+bounds+[10],
                                              extend='both',
                                              # Make the length of each extension
                                              # the same as the length of the
                                              # interior colors:
                                              ticks=bounds,
                                              spacing='uniform',
                                              orientation='horizontal')
              cb3.ax.tick_params(labelsize=fontsize_names)
              #cb3.set_label(deLabel,size=fontsize_names)
              ax3.text(pos_axes.x0+0.01,0.17-indl*0.05,deLabel,size=fontsize_names)
    #Add callback for interaction
    if interactive:
        callback = partial(_plot_connectivity_circle_onpick, fig=fig,
                           axes=axes, indices=indices, n_nodes=n_nodes,
                           node_angles=node_angles)
                           
        fig.canvas.mpl_connect('button_press_event', callback)
    
    return fig, axes


"""
    especially for the long labels which ocuppied huge spaces"""

def plot_connectivity_circle_dir_2para_ll(con, node_names, indices=None, n_lines=None,threshold=0.1,
                             node_angles=None, node_width=None,
                             node_colors=None, facecolor='black',
                             textcolor='white', node_edgecolor='black',barhight=1.,
                             linewidth=1.5, levels=['con','din','dout'],
                             flagcon='Tplot',
                             bottom_lev1=9, colormap='hot', vmin=None,
                             levelmap=['None','deeppink','orangered'],
                             vmax=None, colorbar=True, title=None,
                             colorbar_size=0.3, colorbar_pos=(0.3, 0.1),
                             fontsize_title=12, fontsize_names=8,
                             fontsize_colorbar=8, padding=6.,
                             fig=None, subplot=121, interactive=True,color_given=0):
    """Visualize connectivity as a directed circular graph.
        
        Note: This code is based on the circle graph example by Nicolas P. Rougier
        
        Parameters
        ----------
        con : array
        Connectivity scores. Can be a square matrix, or a 1D array. If a 1D
        array is provided, "indices" has to be used to define the connection
        indices.
        node_names : list of str
        Node names. The order corresponds to the order in con.
        indices : tuple of arrays | None
        Two arrays with indices of connections for which the connections
        strenghts are defined in con. Only needed if con is a 1D array.
        n_lines : int | None
        If not None, only the n_lines strongest connections (strength=abs(con))
        are drawn.
        node_angles : array, shape=(len(node_names,)) | None
        Array with node positions in degrees. If None, the nodes are equally
        spaced on the circle. See mne.viz.circular_layout.
        node_width : float | None
        Width of each node in degrees. If None, the minimum angle between any
        two nodes is used as the width.
        node_colors : list of tuples | list of str
        List with the color to use for each node. If fewer colors than nodes
        are provided, the colors will be repeated. Any color supported by
        matplotlib can be used, e.g., RGBA tuples, named colors.
        facecolor : str
        Color to use for background. See matplotlib.colors.
        textcolor : str
        Color to use for text. See matplotlib.colors.
        node_edgecolor : str
        Color to use for lines around nodes. See matplotlib.colors.
        linewidth : float
        Line width to use for connections.
        colormap : str
        Colormap to use for coloring the connections.
        vmin : float | None
        Minimum value for colormap. If None, it is determined automatically.
        vmax : float | None
        Maximum value for colormap. If None, it is determined automatically.
        colorbar : bool
        Display a colorbar or not.
        title : str
        The figure title.
        colorbar_size : float
        Size of the colorbar.
        colorbar_pos : 2-tuple
        Position of the colorbar.
        fontsize_title : int
        Font size to use for title.
        fontsize_names : int
        Font size to use for node names.
        fontsize_colorbar : int
        Font size to use for colorbar.
        padding : float
        Space to add around figure to accommodate long labels.
        fig : None | instance of matplotlib.pyplot.Figure
        The figure to use. If None, a new figure with the specified background
        color will be created.
        subplot : int | 3-tuple
        Location of the subplot when creating figures with multiple plots. E.g.
        121 or (1, 2, 1) for 1 row, 2 columns, plot 1. See
        matplotlib.pyplot.subplot.
        interactive : bool
        When enabled, left-click on a node to show only connections to that
        node. Right-click shows all connections.
        
        Returns
        -------
        fig : instance of matplotlib.pyplot.Figure
        The figure handle.
        axes : instance of matplotlib.axes.PolarAxesSubplot
        The subplot handle.
        """
    import matplotlib.pyplot as plt
    import matplotlib.path as m_path
    import matplotlib.patches as m_patches
    import matplotlib as mpl
    
    n_nodes = len(node_names)
    con_cp=con
    if node_angles is not None:
        if len(node_angles) != n_nodes:
            raise ValueError('node_angles has to be the same length '
                             'as node_names')
        # convert it to radians
        node_angles = node_angles * np.pi / 180
    else:
        # uniform layout on unit circle
        node_angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    
    if node_width is None:
        # widths correspond to the minimum angle between two nodes
        dist_mat = node_angles[None, :] - node_angles[:, None]
        dist_mat[np.diag_indices(n_nodes)] = 1e9
        node_width = np.min(np.abs(dist_mat))
    else:
        node_width = node_width * np.pi / 180
    
    if node_colors is not None:
        if len(node_colors) < n_nodes:
            node_colors = cycle(node_colors)
    else:
        # assign colors using colormap
        node_colors = [plt.cm.spectral(i / float(n_nodes))
                       for i in range(n_nodes)]
    
    # handle 1D and 2D connectivity information
    if con.ndim == 1:
        if indices is None:
            raise ValueError('indices has to be provided if con.ndim == 1')
    elif con.ndim == 2:
        if con.shape[0] != n_nodes or con.shape[1] != n_nodes:
            raise ValueError('con has to be 1D or a square matrix')
        # we use the lower-triangular part
        con = con
    else:
        raise ValueError('con has to be 1D or a square matrix')
    
    # get the colormap
    if isinstance(colormap, string_types):
        colormap = plt.get_cmap(colormap)
    
    # Make figure background the same colors as axes
    if fig is None:
        fig = plt.figure(figsize=(8, 8), facecolor=facecolor)
    
    # Use a polar axes
    if not isinstance(subplot, tuple):
        subplot = (subplot,)
    axes = plt.subplot(*subplot, polar=True, axisbg=facecolor)

# No ticks, we'll put our own
    plt.xticks([])
    plt.yticks([])
    
    # Set y axes limit, add additonal space if requested
    plt.ylim(0, 10 + padding)
    
    # Remove the black axes border which may obscure the labels
    axes.spines['polar'].set_visible(False)
    
    # Draw lines between connected nodes, only draw the strongest connections
    if n_lines is not None and len(con)*(len(con)-1) > n_lines:
        con_thresh = np.sort(np.abs(con).ravel())[-n_lines]
    else:
        con_thresh = threshold
    #print(n_lines,len(con),con_thresh)
    # get the connections which we are drawing and sort by connection strength
    # this will allow us to draw the strongest connections first

    con_abs = np.abs(con)
    con_draw_idx = np.where(con_abs >= con_thresh)
#print con_draw_idx, con_thresh,

    con = con[con_draw_idx]
    con_abs = con_abs[con_draw_idx]
#indices = [ind[con_draw_idx] for ind in indices]
#print con
# now sort them
    sort_idx = np.argsort(con_abs)
    con_abs = con_abs[sort_idx]
    con = con[sort_idx]
    indices = [ind[sort_idx] for ind in con_draw_idx]
    #print(indices,con)
    # Get vmin vmax for color scaling
    if vmin is None:
        vmin = np.min(con[np.abs(con) >= con_thresh])
    if vmax is None:
        vmax = np.max(con)
    vrange = vmax - vmin
    
    # We want to add some "noise" to the start and end position of the
    # edges: We modulate the noise with the number of connections of the
    # node and the connection strength, such that the strongest connections
    # are closer to the node center
    nodes_n_con = np.zeros((n_nodes), dtype=np.int)
    for i, j in zip(indices[0], indices[1]):
        nodes_n_con[i] += 1
        nodes_n_con[j] += 1

    # initalize random number generator so plot is reproducible
    rng = np.random.mtrand.RandomState(seed=0)

    n_con = len(indices[0])
    noise_max = 0.25 * node_width
    start_noise = rng.uniform(-noise_max, noise_max, n_con)
    end_noise = rng.uniform(-noise_max, noise_max, n_con)
    
    nodes_n_con_seen = np.zeros_like(nodes_n_con)
    for i, (start, end) in enumerate(zip(indices[0], indices[1])):
        nodes_n_con_seen[start] += 1
        nodes_n_con_seen[end] += 1
        
        start_noise[i] *= ((nodes_n_con[start] - nodes_n_con_seen[start])
                           / float(nodes_n_con[start]))
        end_noise[i] *= ((nodes_n_con[end] - nodes_n_con_seen[end])
                                            / float(nodes_n_con[end]))

# scale connectivity for colormap (vmin<=>0, vmax<=>1)
    con_val_scaled = (con - vmin) / vrange
    
    indgree = np.zeros(n_nodes)
    outdgree = np.zeros(n_nodes)
    # Finally, we draw the connections
    for pos, (i, j) in enumerate(zip(indices[0], indices[1])):
        indgree[i] = indgree[i]+1
        outdgree[j] = outdgree[j]+1
        
        # Start point
        t0, r0 = node_angles[j], 10
        
        # End point
        t1, r1 = node_angles[i], 10
        
        # Some noise in start and end point
        t0 += start_noise[pos]
        t1 += end_noise[pos]
        
        verts = [(t0, r0), (t0, 5), (t1, 5), (t1, r1)]
        codes = [m_path.Path.MOVETO, m_path.Path.CURVE4, m_path.Path.CURVE4,
                 m_path.Path.LINETO]
        path = m_path.Path(verts, codes)
                 
                 #print(verts,codes)
        #print(i,j)         

        if color_given==0:
            color = colormap(con_val_scaled[pos])
            
        else:
            color = np.array([48,128,20,255])/255.
                 #plt.plot(verts[3][0],9,marker='^',color=color)
        plt.arrow(t1,r1-1,0,0.1,alpha = 1,head_width=0.05,
                           facecolor = color,edgecolor = color,head_length=0.3)
                 # Actual line
        patch = m_patches.PathPatch(path, fill=False, edgecolor=color,
                                             linewidth=linewidth, alpha=1.)
                 #patch = m_patches.ConnectionPatch(path,arrowstyle=u'-')
                 
        axes.add_patch(patch)
    #axes.arrow(0,0,
    #          0.5,0.5, head_width=0.05, head_length=1, ec='k'
    #         )
    

    height = np.ones(n_nodes) * barhight

    for indl, ilevel in enumerate(levels):
        bottom = bottom_lev1+indl*1.4*barhight# Draw ring with colored nodes
        
        bars = axes.bar(node_angles, height, width=node_width, bottom=bottom,
                        edgecolor=node_edgecolor, lw=2, facecolor='.9',
                        align='center')
                        # Draw ring with different colored nodes
                        
        if ilevel == 'con':
            node_colors_lev=node_colors
        elif ilevel == 'din':
            #indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],indgree)

        elif ilevel == 'dout':
#indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],outdgree)

        for bar, color in zip(bars, node_colors_lev):
            bar.set_facecolor(color)

# Draw node labels
    angles_deg = 180 * node_angles / np.pi
    for name, angle_rad, angle_deg in zip(node_names, node_angles, angles_deg):
        if angle_deg >= 270:
            ha = 'left'
        else:
            # Flip the label, so text is always upright
            angle_deg += 180
            ha = 'right'
        
        axes.text(angle_rad, bottom+1.5, name, size=fontsize_names,
                  rotation=angle_deg, rotation_mode='anchor',
                  horizontalalignment=ha, verticalalignment='center',
                  color=textcolor)
    

    if colorbar:
    # for links
          norm = normalize_colors(vmin=vmin, vmax=vmax)
          sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
          sm.set_array(np.linspace(vmin, vmax))
          cb = plt.colorbar(sm, ax=axes, use_gridspec=False,
                              shrink=colorbar_size,
                              anchor=colorbar_pos)
          cb_yticks = plt.getp(cb.ax.axes, 'yticklabels')
          cb.ax.tick_params(labelsize=fontsize_colorbar)
          plt.setp(cb_yticks, color=textcolor)
          
          # for indegree
          for indl,idegree in enumerate(['in','out']):
              if idegree == 'in':
                  degrees = indgree
                  deLabel = 'indegree'
              else:
                  degrees = outdgree
                  deLabel = 'outdegree'
              pos_axes=axes.get_position()
              
              ax3 = fig.add_axes([pos_axes.x0+indl*0.4, 0.1, 0.27, 0.02])
              
              #ax3 = fig.add_axes([0.05+indl*0.4, 0.1, 0.3, 0.02])
              listmap=mt.mln_generateColormap(levelmap[indl+1])
              cmap = mpl.colors.ListedColormap(listmap)
              bounds = np.linspace(np.min(degrees),np.max(degrees),10).astype(int)
                                                      #print(bounds)
              norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
              cb3 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,
                                              norm=norm,
                                              # boundaries=[-10]+bounds+[10],
                                              extend='both',
                                              # Make the length of each extension
                                              # the same as the length of the
                                              # interior colors:
                                              ticks=bounds,
                                              spacing='uniform',
                                              orientation='horizontal')
                                              #cb3.set_label(deLabel)
              cb3.ax.tick_params(labelsize=fontsize_names)
              ax3.text(pos_axes.x0+0.01,0.17-indl*0.05,deLabel,size=fontsize_names)
    #Add callback for interaction
    if interactive:
        callback = partial(_plot_connectivity_circle_onpick, fig=fig,
                           axes=axes, indices=indices, n_nodes=n_nodes,
                           node_angles=node_angles)
                           
        fig.canvas.mpl_connect('button_press_event', callback)
    
    ax = fig.add_axes([0.4, 0.3, 0.5, 0.5])
    con_cp[con_cp<con_thresh]=0
    cim=plt.imshow(con_cp,cmap='Blues',interpolation='nearest')
    xytickindex=np.arange(0,n_nodes)
    ax.set_xticks(xytickindex)
    ax.set_yticks(xytickindex)
    ax.set_xticklabels(node_names,size=fontsize_names-2,rotation=90)
    ax.set_yticklabels(node_names,size=fontsize_names-2)
    cb = plt.colorbar(cim, ax=ax, use_gridspec=False,
                              shrink=colorbar_size*2,
                              anchor=colorbar_pos)
    if title is not None:
        plt.title(title, color=textcolor, fontsize=fontsize_title,
                  axes=ax)
    return fig, axes

def plot_connectivity_circle_dir_2para_ll_ts(con,eeg, times,node_names, indices=None, n_lines=None,threshold=0.1,
                             node_angles=None, node_width=None,
                             node_colors=None, facecolor='black',
                             textcolor='white', node_edgecolor='black',barhight=1.,
                             linewidth=1.5, levels=['con','din','dout'],
                             flagcon='Tplot',xtickstep=3
                                             ,
                             bottom_lev1=9, colormap='hot', vmin=None,
                             levelmap=['None','deeppink','orangered'],
                             vmax=None, colorbar=True, title=None,
                             colorbar_size=0.3, colorbar_pos=(0.3, 0.1),
                             fontsize_title=12, fontsize_names=8,GivenfixScale =200,
                             fontsize_colorbar=8, padding=6.,
                             fig=None, subplot=221, interactive=True,color_given=0,isxtick=0):
    """Visualize connectivity as a directed circular graph.
        
        HF: also plot the time series below the graph 
        Note: This code is based on the circle graph example by Nicolas P. Rougier
        
        Parameters
        ----------
        con : array
        Connectivity scores. Can be a square matrix, or a 1D array. If a 1D
        array is provided, "indices" has to be used to define the connection
        indices.
        node_names : list of str
        Node names. The order corresponds to the order in con.
        indices : tuple of arrays | None
        Two arrays with indices of connections for which the connections
        strenghts are defined in con. Only needed if con is a 1D array.
        n_lines : int | None
        If not None, only the n_lines strongest connections (strength=abs(con))
        are drawn.
        node_angles : array, shape=(len(node_names,)) | None
        Array with node positions in degrees. If None, the nodes are equally
        spaced on the circle. See mne.viz.circular_layout.
        node_width : float | None
        Width of each node in degrees. If None, the minimum angle between any
        two nodes is used as the width.
        node_colors : list of tuples | list of str
        List with the color to use for each node. If fewer colors than nodes
        are provided, the colors will be repeated. Any color supported by
        matplotlib can be used, e.g., RGBA tuples, named colors.
        facecolor : str
        Color to use for background. See matplotlib.colors.
        textcolor : str
        Color to use for text. See matplotlib.colors.
        node_edgecolor : str
        Color to use for lines around nodes. See matplotlib.colors.
        linewidth : float
        Line width to use for connections.
        colormap : str
        Colormap to use for coloring the connections.
        vmin : float | None
        Minimum value for colormap. If None, it is determined automatically.
        vmax : float | None
        Maximum value for colormap. If None, it is determined automatically.
        colorbar : bool
        Display a colorbar or not.
        title : str
        The figure title.
        colorbar_size : float
        Size of the colorbar.
        colorbar_pos : 2-tuple
        Position of the colorbar.
        fontsize_title : int
        Font size to use for title.
        fontsize_names : int
        Font size to use for node names.
        fontsize_colorbar : int
        Font size to use for colorbar.
        padding : float
        Space to add around figure to accommodate long labels.
        fig : None | instance of matplotlib.pyplot.Figure
        The figure to use. If None, a new figure with the specified background
        color will be created.
        subplot : int | 3-tuple
        Location of the subplot when creating figures with multiple plots. E.g.
        121 or (1, 2, 1) for 1 row, 2 columns, plot 1. See
        matplotlib.pyplot.subplot.
        interactive : bool
        When enabled, left-click on a node to show only connections to that
        node. Right-click shows all connections.
        
        Returns
        -------
        fig : instance of matplotlib.pyplot.Figure
        The figure handle.
        axes : instance of matplotlib.axes.PolarAxesSubplot
        The subplot handle.
        """
    import matplotlib.pyplot as plt
    import matplotlib.path as m_path
    import matplotlib.patches as m_patches
    import matplotlib as mpl
    
    n_nodes = len(node_names)
    con_cp=con
    if node_angles is not None:
        if len(node_angles) != n_nodes:
            raise ValueError('node_angles has to be the same length '
                             'as node_names')
        # convert it to radians
        node_angles = node_angles * np.pi / 180
    else:
        # uniform layout on unit circle
        node_angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    
    if node_width is None:
        # widths correspond to the minimum angle between two nodes
        dist_mat = node_angles[None, :] - node_angles[:, None]
        dist_mat[np.diag_indices(n_nodes)] = 1e9
        node_width = np.min(np.abs(dist_mat))
    else:
        node_width = node_width * np.pi / 180
    
    if node_colors is not None:
        if len(node_colors) < n_nodes:
            node_colors = cycle(node_colors)
    else:
        # assign colors using colormap
        node_colors = [plt.cm.spectral(i / float(n_nodes))
                       for i in range(n_nodes)]
    
    # handle 1D and 2D connectivity information
    if con.ndim == 1:
        if indices is None:
            raise ValueError('indices has to be provided if con.ndim == 1')
    elif con.ndim == 2:
        if con.shape[0] != n_nodes or con.shape[1] != n_nodes:
            raise ValueError('con has to be 1D or a square matrix')
        # we use the lower-triangular part
        con = con
    else:
        raise ValueError('con has to be 1D or a square matrix')
    
    # get the colormap
    if isinstance(colormap, string_types):
        colormap = plt.get_cmap(colormap)
    
    # Make figure background the same colors as axes
    if fig is None:
        fig = plt.figure(figsize=(20, 16), facecolor=facecolor)
    
    # Use a polar axes
    if not isinstance(subplot, tuple):
        subplot = (subplot,)
    axes = plt.subplot(*subplot, polar=True, axisbg=facecolor)
    axes.set_position([0.05,0.54,0.45,0.4])
    
# No ticks, we'll put our own
    plt.xticks([])
    plt.yticks([])
    
    # Set y axes limit, add additonal space if requested
    plt.ylim(0, 10 + padding)
    
    # Remove the black axes border which may obscure the labels
    axes.spines['polar'].set_visible(False)
    
    # Draw lines between connected nodes, only draw the strongest connections
    if n_lines is not None and len(con)*(len(con)-1) > n_lines:
        con_thresh = np.sort(np.abs(con).ravel())[-n_lines]
    else:
        con_thresh = threshold
    #print(n_lines,len(con),con_thresh)
    # get the connections which we are drawing and sort by connection strength
    # this will allow us to draw the strongest connections first

    con_abs = np.abs(con)
    con_draw_idx = np.where(con_abs >= con_thresh)
#print con_draw_idx, con_thresh,

    con = con[con_draw_idx]
    con_abs = con_abs[con_draw_idx]
#indices = [ind[con_draw_idx] for ind in indices]
#print con
# now sort them
    sort_idx = np.argsort(con_abs)
    con_abs = con_abs[sort_idx]
    con = con[sort_idx]
    indices = [ind[sort_idx] for ind in con_draw_idx]
    #print(indices,con)
    # Get vmin vmax for color scaling
    if vmin is None:
        vmin = np.min(con[np.abs(con) >= con_thresh])
    if vmax is None:
        vmax = np.max(con)
    vrange = vmax - vmin
    
    # We want to add some "noise" to the start and end position of the
    # edges: We modulate the noise with the number of connections of the
    # node and the connection strength, such that the strongest connections
    # are closer to the node center
    nodes_n_con = np.zeros((n_nodes), dtype=np.int)
    for i, j in zip(indices[0], indices[1]):
        nodes_n_con[i] += 1
        nodes_n_con[j] += 1

    # initalize random number generator so plot is reproducible
    rng = np.random.mtrand.RandomState(seed=0)

    n_con = len(indices[0])
    noise_max = 0.25 * node_width
    start_noise = rng.uniform(-noise_max, noise_max, n_con)
    end_noise = rng.uniform(-noise_max, noise_max, n_con)
    
    nodes_n_con_seen = np.zeros_like(nodes_n_con)
    for i, (start, end) in enumerate(zip(indices[0], indices[1])):
        nodes_n_con_seen[start] += 1
        nodes_n_con_seen[end] += 1
        
        start_noise[i] *= ((nodes_n_con[start] - nodes_n_con_seen[start])
                           / float(nodes_n_con[start]))
        end_noise[i] *= ((nodes_n_con[end] - nodes_n_con_seen[end])
                                            / float(nodes_n_con[end]))

# scale connectivity for colormap (vmin<=>0, vmax<=>1)
    con_val_scaled = (con - vmin) / vrange
    
    indgree = np.zeros(n_nodes)
    outdgree = np.zeros(n_nodes)
    # Finally, we draw the connections
    for pos, (i, j) in enumerate(zip(indices[0], indices[1])):
        indgree[i] = indgree[i]+1
        outdgree[j] = outdgree[j]+1
        
        # Start point
        t0, r0 = node_angles[j], 10
        
        # End point
        t1, r1 = node_angles[i], 10
        
        # Some noise in start and end point
        t0 += start_noise[pos]
        t1 += end_noise[pos]
        
        verts = [(t0, r0), (t0, 5), (t1, 5), (t1, r1)]
        codes = [m_path.Path.MOVETO, m_path.Path.CURVE4, m_path.Path.CURVE4,
                 m_path.Path.LINETO]
        path = m_path.Path(verts, codes)
                 
                 #print(verts,codes)
        #print(i,j)         

        if color_given==0:
            color = colormap(con_val_scaled[pos])
            
        else:
            color = np.array([48,128,20,255])/255.
                 #plt.plot(verts[3][0],9,marker='^',color=color)
        plt.arrow(t1,r1-1,0,0.1,alpha = 1,head_width=0.05,
                           facecolor = color,edgecolor = color,head_length=0.3)
                 # Actual line
        patch = m_patches.PathPatch(path, fill=False, edgecolor=color,
                                             linewidth=linewidth, alpha=1.)
                 #patch = m_patches.ConnectionPatch(path,arrowstyle=u'-')
                 
        axes.add_patch(patch)
    #axes.arrow(0,0,
    #          0.5,0.5, head_width=0.05, head_length=1, ec='k'
    #         )
    

    height = np.ones(n_nodes) * barhight

    for indl, ilevel in enumerate(levels):
        bottom = bottom_lev1+indl*1.4*barhight# Draw ring with colored nodes
        
        bars = axes.bar(node_angles, height, width=node_width, bottom=bottom,
                        edgecolor=node_edgecolor, lw=2, facecolor='.9',
                        align='center')
                        # Draw ring with different colored nodes
                        
        if ilevel == 'con':
            node_colors_lev=node_colors
        elif ilevel == 'din':
            #indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],indgree)

        elif ilevel == 'dout':
#indgree, outdgree = mln_degree(con2D,flagcon,threshold)
            node_colors_lev=mt.mln_nodeColors(levelmap[indl],outdgree)

        for bar, color in zip(bars, node_colors_lev):
            bar.set_facecolor(color)

# Draw node labels
    angles_deg = 180 * node_angles / np.pi
    for name, angle_rad, angle_deg in zip(node_names, node_angles, angles_deg):
        if angle_deg >= 270:
            ha = 'left'
        else:
            # Flip the label, so text is always upright
            angle_deg += 180
            ha = 'right'
        
        axes.text(angle_rad, bottom+1.5, name, size=fontsize_names,
                  rotation=angle_deg, rotation_mode='anchor',
                  horizontalalignment=ha, verticalalignment='center',
                  color=textcolor)
    

    if colorbar:
    # for links
          norm = normalize_colors(vmin=vmin, vmax=vmax)
          sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
          sm.set_array(np.linspace(vmin, vmax))
          cb = plt.colorbar(sm, ax=axes, use_gridspec=False,
                              shrink=colorbar_size,
                              anchor=colorbar_pos)
          cb_yticks = plt.getp(cb.ax.axes, 'yticklabels')
          cb.ax.tick_params(labelsize=fontsize_colorbar)
          plt.setp(cb_yticks, color=textcolor)
          
          # for indegree
          for indl,idegree in enumerate(['in','out']):
              if idegree == 'in':
                  degrees = indgree
                  deLabel = 'indegree'
              else:
                  degrees = outdgree
                  deLabel = 'outdegree'
              pos_axes=axes.get_position()
              
              ax3 = fig.add_axes([pos_axes.x0+indl*0.4, 0.47, 0.27, 0.01])
              
              #ax3 = fig.add_axes([0.05+indl*0.4, 0.1, 0.3, 0.02])
              listmap=mt.mln_generateColormap(levelmap[indl+1])
              cmap = mpl.colors.ListedColormap(listmap)
              bounds = np.linspace(np.min(degrees),np.max(degrees),10).astype(int)
                                                      #print(bounds)
              norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
              cb3 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,
                                              norm=norm,
                                              # boundaries=[-10]+bounds+[10],
                                              extend='both',
                                              # Make the length of each extension
                                              # the same as the length of the
                                              # interior colors:
                                              ticks=bounds,
                                              spacing='uniform',
                                              orientation='horizontal')
                                              #cb3.set_label(deLabel)
              cb3.ax.tick_params(labelsize=fontsize_names)
              ax3.text(pos_axes.x0+0.01,1.1,deLabel,size=fontsize_names)
    #Add callback for interaction
    if interactive:
        callback = partial(_plot_connectivity_circle_onpick, fig=fig,
                           axes=axes, indices=indices, n_nodes=n_nodes,
                           node_angles=node_angles)
                           
        fig.canvas.mpl_connect('button_press_event', callback)
    ax = fig.add_axes([0.52, 0.54, 0.44, 0.4])
    con_cp[con_cp<con_thresh]=0
    cim=plt.imshow(con_cp,cmap='Blues',interpolation='nearest')
    xytickindex=np.arange(0,n_nodes)
    ax.set_xticks(xytickindex)
    ax.set_yticks(xytickindex)
    ax.set_xticklabels(node_names,size=fontsize_names-2,rotation=90)
    ax.set_yticklabels(node_names,size=fontsize_names-2)
    cb = plt.colorbar(cim, ax=ax, use_gridspec=False,
                              shrink=colorbar_size*2,
                              anchor=colorbar_pos)
    if title is not None:
        plt.title(title, color=textcolor, fontsize=fontsize_title,
                  axes=ax)

    ax1 = fig.add_axes([0.1, 0.03, 0.84
                        , 0.4])
    for i, y in enumerate(eeg[:, :]):
    
        ax1.plot(times,(y - y.mean())/GivenfixScale + i)
        ScaleMark=GivenfixScale
    
#xticks_given=np.arange(int(times[0]),int(max(times))+1,xtickstep)
    ax1.set_yticklabels(node_names,fontsize=20)
    ax1.set_yticks(range(0,len(node_names)))
#   ax1.set_xlim([times[0],max(times)])

    ax1.set_ylim([-1,len(node_names)])
    ax1.invert_yaxis()
    if isxtick==0:
        ax1.set_xticks([])
        removespine=['left', 'bottom', 'right', 'top']
        for ispines in removespine:
            ax1.spines[ispines].set_color('none')
    else:
        removespine=['left', 'right', 'top']
        for ispines in removespine:
            ax1.spines[ispines].set_color('none')
        timeX=np.arange(int(np.round(times[0])),int(np.round(times[-1]))+xtickstep,xtickstep)
        print(timeX)
        ax1.set_xticks(timeX)
        ax1.set_xlim([times[0],max(times)])
        ax1.set_xticklabels(timeX,fontsize=16)
    return fig, axes







    