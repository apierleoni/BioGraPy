'''
Created on 28/ott/2009

@author: andrea pierleoni

quote:
That's what I like in matplotlib: no matter how hard you try,
there's always a simpler solution you're not yet aware of...


'''

import matplotlib, warnings, operator
matplotlib.use('Agg')
import biograpy.tracks
from matplotlib.font_manager import FontProperties

warnings.simplefilter("ignore")

class Panel(object):
    '''
 
    
    The :class:`Panel` acts as a collector of :class:`~biograpy.tracks` objects,
    and accounts for their correct arrangement and  figure saving
    A :class:`Panel` must be initialized before the tracks and features are 
    created to link the drawn artists to the :class:`Panel`


    ``panel = Panel(fig_width=800, fig_height=600, fig_dpi=80)``

    gives a Panel object with a 800x600 pixel canvas at 80 dpi resolution

    Valid keyword arguments for :class:`Panel` are:

        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        fig_height            figure height in pixel, default is ``1500``
        fig_width             figure width in pixel, default is ``None`` and 
                              will be autodetermined
        fig_dpi               default is ``80``
        grid                  ``'major`` | ``'minor'`` | ``'both'`` | ``None`` 
                              default is ``'major'``
        join_tracks           ``True`` | ``False`` default is ``False``.
                              Will join all tracks showing just the last bottom
                              axis
        track_padding         track padding in pixel default is ``0``
        start_position        start int number for axis tick labels default is 
                              ``0``
        padding               padding of all the axis with respect to the figure
                              border  default is ``fig_width*.02``
        figure_bottom_space   define the distance of the last bottom axis from 
                              the figure bottom in pixel. Default is ``0`` . Use
                              to correct too long figure width coming from 
                              autodetermined `fig_width` 
        xmin                  int minimum value of the X axes, use to plot just 
                              a part of the drawing, default is ``None``
        xmax                  int maximum value of the X axes, use to plot just 
                              a part of the drawing, default is ``None``  
        ===================== ==================================================
        
    '''


    def __init__(self,  fig_width=1500, fig_height=None, fig_dpi=80, **kwargs):
        """ """
        self.fig_width = fig_width / float(fig_dpi)
        # fig_height are handled by matplotlib in inches not pixel
        if fig_height:
            self.fig_height = fig_height / float(fig_dpi)
        else:
            self.fig_height = fig_height
        self.dpi = float(fig_dpi)
        self.features = {}
        self.tracks = []
        self.grid = kwargs.get('grid', 'major') # can be "major", "minor", "both" or None
        self.join_tracks = kwargs.get('join_tracks', False)
        if self.join_tracks:
            self.track_padding = 0.
        else:
            self.track_padding = kwargs.get('track_padding', 0)
        self.start_position = kwargs.get('start_position', 0) #changes the start position by adding this value to all x coords
        self.padding = kwargs.get('padding', fig_width*.02)#pixel to be used as padding in all 4 directions
        self.hpadding = self.padding / float(fig_width)
        self.figure_bottom_space = kwargs.get('figure_bottom_space', 0)
        if self.fig_height:
            self.vpadding = self.padding / float(fig_height)
            self.vtrack_padding = self.track_padding / float(fig_height)
        else:
            self.vpadding = self.hpadding
        self.xmin = kwargs.get('xmin', None)
        self.xmax = kwargs.get('xmax', None)
        
            
        '''create figure object'''
        self.fig = matplotlib.pyplot.figure(1, figsize = (self.fig_width, self.fig_width), dpi = self.dpi, frameon = False)
        self.ax = self.fig.add_subplot(111)#needed to make it invisible
        ''' '''

    def add_track(self, track):
        """
        Accepts any :class:`~biograpy.tracks` object.
        the track should already contains its features
        """
        self.tracks.append(track)
        
    def append(self, *args, **kwargs):
        """
        Accepts any :class:`~biograpy.tracks` object.
        the track should already contains its features
        """
        self.add_track(*args, **kwargs)
        
    def extend(self, tracks):
        """
        Accepts a sequence of :class:`~biograpy.tracks` objects
        """
        for track in tracks:
            self.add_track(track)

    def _estimate_fig_height(self,):
        return self.drawn_lines*self.fig_width/30.

    def _draw_tracks(self, **kwargs):
        '''create an axis for each track and moves
        accordingly all the child features'''
        
        self.ax.set_axis_off()
        self.Drawn_objects = []
        self.track_axes = []
        draw_xmin = kwargs.get('xmin', None)
        draw_xmax = kwargs.get('xmax', None)
        if draw_xmin:
            self.xmin = draw_xmin
        if draw_xmax:
            self.xmax = draw_xmax
        '''estimate track height using track.drawn_lines
        find max and min x coords, and check for colorbar
        presence in at least one track'''
        cbars = False
        self.drawn_lines = 0
        Xs =[]
        track_height_user_specified = False
        for track in self.tracks:
            if track.features:#skip tracks with no features
                if track.track_height and self.fig_height: #if not fig_height is specified user track heights are ignored
                    track_height_user_specified = True
                if track_height_user_specified and not track.track_height:
                    track_height_user_specified = False #disable if some track has not a specified heigth
                    warnings.warn('All tracks need to have a specified track_height, reverting to automatic track height')
                track._sort_features(dpi = self.dpi, xoffset = self.xmin)#THIS WILL DRAW ALL THE FEATURES
                if track.draw_cb:
                    if cbars != 'label':
                        cbars = 'simple'
                        if track.cb_label:
                            cbars = 'label'
                self.drawn_lines += track.drawn_lines
                Xs.append(track.xmin)
                Xs.append(track.xmax)
        if self.xmin == None:
            self.xmin = min(Xs)
        if self.xmax == None:
            self.xmax =max(Xs)
        '''auto estimate fig_heigth and panning if needed '''
        if not self.fig_height:#automatcally set fig height basing on the total number of features
            self.fig_height=self._estimate_fig_height()
            self.vpadding = (float(self.padding)/self.dpi) / self.fig_height
            self.vtrack_padding = (float(self.track_padding)/self.dpi) / self.fig_height
            
        '''set colorbar dimension '''
        if cbars == 'label':
            cbar_extent=0.015
            cbar_axis_space = 0.05
            cbar_right_pad = 0.03
        else:
            cbar_extent=0.015
            cbar_axis_space = 0.05
            cbar_right_pad = 0.01
            
        '''arrange tracks'''
        
        axis_left_pad = self.hpadding
        default_figure_bottom_space = self.vpadding + float(self.figure_bottom_space)/self.dpi
        axis_bottom_pad = 1.0 - default_figure_bottom_space
        axis_width = 1.-2*self.hpadding
        axis_scale = None# used to persist the same scale on all the tracks
        if cbars:
            axis_width -= (cbar_extent + cbar_axis_space + cbar_right_pad)
        '''cycle trought tracks and draw them as axix object '''
        #canvas_height = 0
        for track_num, track in enumerate(self.tracks):
            if track.features:#skip tracks with no features
                '''define axis dimensions and position and create axis'''
                if track_height_user_specified:
                    axis_height = track.track_height /float(self.fig_height)
                else:
                    if axis_scale:
                        axis_height = axis_scale * track.drawn_lines
                    else:    
                        axis_height = (float(track.drawn_lines)/self.drawn_lines)  - self.vpadding/(2.*len(self.tracks)) - default_figure_bottom_space/len(self.tracks)
                        axis_scale = axis_height / float(track.drawn_lines)
                axis_bottom_pad -= (axis_height + self.vtrack_padding/2.)
                axis = matplotlib.pyplot.axes([axis_left_pad,axis_bottom_pad, axis_width, axis_height ],) 
                self.track_axes.append(axis)
                
                
                '''handle track axis display, ticks and tickslabel '''
                '''set Y lims '''
                if isinstance(track, tracks.PlotTrack):
                    if track.show_name:
                        if track.show_name == 'top':
                            axis.set_ylim(track.Ycord, track.ymax+1)
                            track_name = axis.text(self.xmin + (self.xmax * 0.01), track.ymax + .5, track.name,  horizontalalignment='left', verticalalignment='bottom', fontproperties=track.name_font_feat,)
                        elif track.show_name == 'bottom':
                            axis.set_ylim(track.Ycord-1, track.ymax)
                            track_name = axis.text(self.xmin + (self.xmax * 0.01), track.Ycord - .5, track.name,  horizontalalignment='left', verticalalignment='bottom', fontproperties=track.name_font_feat,)
                    else:
                        axis.set_ylim(track.Ycord, track.ymax)
                else:
                    if track.show_name:
                        if track.show_name == 'top':
                            axis.set_ylim(track.Ycord, track.ymax+2.5)
                            track_name = axis.text(self.xmin + (self.xmax * 0.01), track.ymax + 1.5, track.name,  horizontalalignment='left', verticalalignment='bottom', fontproperties=track.name_font_feat,)
                        elif track.show_name == 'bottom':
                            axis.set_ylim(track.Ycord-1, track.ymax + 1.5)
                            track_name = axis.text(self.xmin + (self.xmax * 0.01), track.Ycord - 0., track.name,  horizontalalignment='left', verticalalignment='bottom', fontproperties=track.name_font_feat,)
                    else: 
                        axis.set_ylim(track.Ycord, track.ymax+1.5,)
                '''set X lims'''
                axis.set_xlim(self.xmin, self.xmax)
                
                '''handle last bottom axis '''
                if (track_num+1 == len(self.tracks)) and ('force no axis' not in track.draw_axis):
                    track.draw_axis.append('bottom')
                if not self.track_padding:
                    if (track_num+1 != len(self.tracks)):
                        if 'bottom' in track.draw_axis:
                            del track.draw_axis[track.draw_axis.index('bottom')]
                        if 'top' in track.draw_axis:
                            del track.draw_axis[track.draw_axis.index('top')]
                    axis.spines["top"].set_color('none')

                '''handle axis and ticks'''
                for spine in ["right", "left", "top", "bottom"]:
                    if spine not in track.draw_axis:
                        axis.spines[spine].set_color('none')# don't draw axis 
                if ("right" not in track.draw_axis) and ("left" not in track.draw_axis):
                    axis.yaxis.set_ticks([])# dont'show ticks and labels on y
                if 'top' not in track.draw_axis:
                    axis.xaxis.set_ticks_position('bottom')# only show bottom ticks
                if 'right' not in track.draw_axis:
                    axis.yaxis.set_ticks_position('left')# only show left ticks
                
                '''handle X ticks and labels '''
                step=int(round(self.xmax/10.,1-len(str(int(self.xmax / 10.)))))
                auto_X_major_ticks = range(self.xmin, self.xmax + 1, step)
                step_min = step / 4.
                tick = auto_X_major_ticks[0]
                auto_X_minor_ticks =[]
                while tick <= self.xmax:
                    auto_X_minor_ticks.append(int(round(tick)))
                    tick+= step_min

                

                '''use sequence as X ticks '''
                if track.x_use_sequence:
                    if len(track.x_use_sequence) < self.xmax-self.xmin:
                        raise Exception('Sequence must be of the same length of X coords')
                    track.xticks_minor = []
                    track.xticklabels_minor = []
                    for i, seq in  enumerate(track.x_use_sequence):
                        if self.xmin+i <= self.xmax:
                            track.xticks_minor.append(self.xmin+i)
                            track.xticklabels_minor.append(str(seq))
                    
                
                '''major X ticks '''
                X_major_ticks_labels = None
                if track.xticks_major != None:
                    X_major_ticks = track.xticks_major
                    if (track.xticklabels_major != None) \
                       and len(track.xticklabels_major) == len(track.xticks_major):
                        X_major_ticks_labels = track.xticklabels_major
                else:
                    X_major_ticks = auto_X_major_ticks
                if 'bottom' in track.draw_axis :
                    axis.set_xticks(X_major_ticks)
                else:
                    axis.set_xticks([])
                '''major ticks labels '''
                if (track_num+1 == len(self.tracks)) or track.show_xticklabels:# last track or forced display
                    if X_major_ticks_labels == None:
                        X_major_ticks_labels = []
                        for i in X_major_ticks:
                            if isinstance(i, (float, int)):
                                X_major_ticks_labels.append(i+ self.start_position)
                            else:
                                X_major_ticks_labels.append(i)
                    axis.set_xticklabels(X_major_ticks_labels, fontsize=track.tickfontsize)
                    if track.x_use_sequence:
                        try:
                            axis.xaxis.set_tick_params(pad = 15, )
                        except: pass #not supported in matplotlib <1
                else:
                    axis.set_xticklabels([])
                
                '''minor X ticks '''
                X_minor_ticks_labels = None
                if track.xticks_minor != None:
                    X_minor_ticks = track.xticks_minor
                    if (track.xticklabels_minor != None) \
                       and len(track.xticklabels_minor) == len(track.xticks_minor):
                        X_minor_ticks_labels = track.xticklabels_minor
                else:
                    X_minor_ticks = auto_X_minor_ticks
                if 'bottom' in track.draw_axis :
                    axis.set_xticks(X_minor_ticks, minor=True)
                else:
                    axis.set_xticks([], minor=True)
                '''minor ticks labels '''
                if (track_num+1 == len(self.tracks)) or track.show_xticklabels:# last track or forced display
                    if X_minor_ticks_labels == None:
                        X_minor_ticks_labels = []
                        '''for i in X_minor_ticks:
                            if i in X_major_ticks:
                                X_minor_ticks_labels.append('')
                            else:
                                label = ''
                                if isinstance(i, (float, int)):
                                    label = str(i+ self.start_position)
                                else:
                                    label = i
                                if (len(str(i).split('.')[0])>=4) : #avoid too long minor ticks
                                    label = ''
                                    X_minor_ticks_labels = []
                                    break#no minor ticks displayed
                                X_minor_ticks_labels.append(label)'''
                    axis.set_xticklabels(X_minor_ticks_labels, fontsize=track.tickfontsize_minor, minor=True)
                else:
                    axis.set_xticklabels([], minor=True)
                    
                '''handle Y ticks and labels '''
                '''major Y ticks '''
                Y_major_ticks_labels = None
                if track.yticks_major != None:
                    Y_major_ticks = track.yticks_major
                    if (track.yticklabels_major != None) \
                       and len(track.yticklabels_major) == len(track.yticks_major):
                        Y_major_ticks_labels = track.yticklabels_major
                else:
                    Y_major_ticks = None
                if ('left' in track.draw_axis)  and track.yticks_major:
                    axis.set_yticks(Y_major_ticks)
                '''major ticks labels '''
                if Y_major_ticks and track.show_yticklabels:
                    if Y_major_ticks_labels == None:
                        Y_major_ticks_labels = []
                        for i in Y_major_ticks:
                            Y_major_ticks_labels.append(i)
                    axis.set_yticklabels(Y_major_ticks_labels, fontsize=track.tickfontsize)
                else:
                    try:
                        axis.yaxis.set_tick_params(labelsize = track.tickfontsize)
                    except: pass #not supported in matplotlib <1
                '''minor Y ticks '''
                Y_minor_ticks_labels = None
                if track.yticks_minor != None:
                    Y_minor_ticks = track.yticks_minor
                    if (track.yticklabels_minor != None) \
                       and len(track.yticklabels_minor) == len(track.yticks_minor):
                        Y_minor_ticks_labels = track.yticklabels_minor
                else:
                    Y_minor_ticks = None
                if ('left' in track.draw_axis)  and track.yticks_minor:
                    axis.set_yticks(Y_minor_ticks, minor=True)
                '''minor ticks labels '''
                if Y_minor_ticks and track.show_yticklabels:
                    if Y_minor_ticks_labels == None:
                        Y_minor_ticks_labels = []
                        for i in Y_minor_ticks:
                            if i in Y_major_ticks:
                                Y_minor_ticks_labels.append('')
                            else:
                                Y_minor_ticks_labels.append(i)
                    axis.set_yticklabels(Y_minor_ticks_labels, fontsize=track.tickfontsize_minor, minor=True)
                else:
                    try:
                        axis.yaxis.set_tick_params(which= 'minor', labelsize = track.tickfontsize)
                    except: pass #not supported in matplotlib <1
                       
                    

                '''draw grid'''
                if self.grid:
                    if (self.grid == 'major') or (self.grid == 'both'):
                        for X in auto_X_major_ticks:
                            axis.axvline(X,ls=':',c='grey',alpha=0.66, zorder = -1)
                    if (self.grid == 'minor') or (self.grid == 'both'):
                        for X in auto_X_minor_ticks:
                            axis.axvline(X,ls=':',c='grey',alpha=0.33, zorder = -1)
                
                
                '''add feature patches to track axes '''
                for feature in track.features:
                    self.Drawn_objects.append(feature)
                    for patch in feature.patches:
                        if isinstance(patch, matplotlib.lines.Line2D):
                            axis.add_line(patch)
                        elif isinstance(patch, matplotlib.patches.Patch):
                            axis.add_patch(patch)
                        else:
                            axis.add_artist(patch)
                        patch.set_transform(axis.transData)# IMPORTANT WORKAROUND!!! if not manually set, transform is not passed correctly in Line2D objects
                                                
                    for feat_name in feature.feat_name:
                        axis.add_artist(feat_name)

                if track.draw_cb:
                    cb_axis = matplotlib.pyplot.axes([axis_left_pad + axis_width + cbar_axis_space - cbar_right_pad ,axis_bottom_pad, cbar_extent, axis_height ],) 
                    if (track.min_score == None) and (track.max_score == None):
                        for feat in track.features:
                            if feat.norm != None:
                                track.norm = feat.norm
                                break

                    cb1 = matplotlib.colorbar.ColorbarBase(cb_axis, cmap=track.cm,
                                                       norm=track.norm,
                                                       alpha = track.cb_alpha,
                                                       orientation='vertical')
                    if track.cb_label:
                        cb1.set_label(track.cb_label)
                    #cb_axis.axes.set_axis_off()
                    for label in cb_axis.get_yticklabels():
                        label.set_fontsize('xx-small')
                        
                '''handle legend '''
                legend_font=FontProperties()
                legend_font.set_size('x-small')
                legend_font.set_family('serif')
                legend_font.set_weight('normal')
                axis.legend(prop = legend_font)
                
        '''set panel size and panning '''
        self.fig.set_figheight(self.fig_height)
        self.fig.set_figwidth(self.fig_width)

    def _boxes(self):
        '''must be called after Drawer.save(output)
        '''
        if not getattr(self, 'Drawn_objects', None):
            self._draw_tracks()
        trans = self.fig.get_transform() # transform should not be necessary if the plot was already plotted
        dpi = self.fig.get_dpi()
        # XXX: img_width unused ???
        # img_width = self.fig.get_figwidth() * dpi
        img_height = self.fig.get_figheight() * dpi
        for obj in self.Drawn_objects:
            xs_patches=[]
            ys_patches=[]
            for patch in obj.patches:
                if isinstance(patch,list):
                    for p in patch:
                        bbox=p.get_window_extent(None)
                        xs_patches.append(bbox.xmax)
                        xs_patches.append(bbox.xmin)
                        ys_patches.append(bbox.ymax)
                        ys_patches.append(bbox.ymin)
                else:
                    try:
                        bbox=patch.get_window_extent(None)
                        xs_patches.append(bbox.xmax)
                        xs_patches.append(bbox.xmin)
                        ys_patches.append(bbox.ymax)
                        ys_patches.append(bbox.ymin)
                    except:
                        xs_patches = ys_patches = False
            if xs_patches and ys_patches:
                xmin, ymin = trans.transform([min(xs_patches),min(ys_patches)])
                xmax, ymax = trans.transform([max(xs_patches),max(ys_patches)])
                left = xmin
                top = img_height-ymin
                right = xmax
                bottom = img_height-ymax
                yield dict(feature=obj, left=left, top=top, right=right, bottom=bottom, track=None, proceed = True)
            else:
                warnings.warn('could not find box coordinated for patch: '+str(patch) )
                yield dict(feature=obj, left=left, top=top, right=right, bottom=bottom, track=None, proceed = False)

    def _create_html_map(self, map_name = 'biograpy-map', map_id = 'biograpy-map', target = '_self', **kwargs):
        """
        returns the corresponding html map from self.Drawn_objects in self.htmlmap
        target ---> set target="_blank" on area links if needed
        """
        areas =[]
        for box in self._boxes():
            if box['proceed']:
                obj = box['feature']
                area_dict = dict(shape = 'rect', # shape: rect, circle, poly
                                 coords = '%(left)i,%(top)i,%(right)i,%(bottom)i' % box,
                                 href = obj.url or '#%s'%obj.name, #href
                                 target = target,
                                 script = obj.html_map_extend,
                                 alt = obj.name )
                area_html = '''<area shape="%(shape)s" coords="%(coords)s" href="%(href)s" target="%(target)s" alt="%(alt)s" %(script)s >''' % area_dict
                areas.append(area_html)
        self.htmlmap ='''<map name="%s" id="%s">\n %s \n</map>'''%(map_name, map_id, '\n'.join(areas))
        return



    def save(self, output, html_target='_self', xmin = None, xmax = None, **kwargs):
        '''
        
        input parameters are:
    
        =========== ============================================================
        Property    Description
        =========== ============================================================
        output      a string containing the file path or a file-like handler
        format      a format must be specified if a file handler is passed 
                    supported formats are: emf, eps, pdf, png, ps, raw, rgba, 
                    svg, svgz 
        html_target ``'_self'`` | ``'_blank'`` | ``'_parent'`` | ``'_top'`` 
                    default is ``'_self'`` . The html target value for 
                    generated hyperlinks. set to  ``'_blank'`` to open in a new
                    browser windows
        xmin        int minimum value of the X axes, use to plot just 
                    a part of the drawing, default is ``None``
        xmax        int maximum value of the X axes, use to plot just 
                    a part of the drawing, default is ``None``
        =========== ============================================================


        '''
        self._draw_tracks(xmin = xmin, xmax = xmax)
        create_html_map = False
        if kwargs.get('format',None) in ('png', 'jpg', 'jpeg'):
            create_html_map = True
        elif isinstance(output,str):
            if output[-3:] in  ('png', 'jpg', 'jpeg'):
                create_html_map = True
        if create_html_map:
            self._create_html_map(target = html_target, **kwargs)# do it just for png file
        
        matplotlib.pyplot.savefig(output, dpi=self.fig.get_dpi(), **kwargs)
        
    
    def close(self):
        '''Close to free the panel. Use it before starting a new drawing in the \
        same process. Typical usage scenario is a web server.'''
        
        matplotlib.pyplot.close()
        

