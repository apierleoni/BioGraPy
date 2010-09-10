'''
Created on 28/ott/2009

@author: andreapierleoni

quote:
That's what I like in matplotlib: no matter how hard you try,
there's always a simpler solution you're not yet aware of...


'''

import matplotlib
matplotlib.use('Agg')


class Panel(object):
    '''
::    
    
    Act as a collector of Tracks objects, and accounts for their correct arrangement and  figure saving
    Drawer MUST BE INITIALIZED before the tracks and features are created to link the drawn artists to the Panel

    self.add_track ---> add a track to the panel
    self.save ---> saves the figure to the specified output
    (self._create_html_map ---> returns an HTML map of the image and its features.)
    

    Kwargs:
        fig_width (int): (default 800)
        fig_height (int): (default None)
        fig_dpi (int): (default 80)

    >>> Panel(fig_width=800, fig_height=600, fig_dpi=300)

    gives a Panel object with a 800x600 pixel canvas at 300 dpi resolution

    
    '''


    def __init__(self,  fig_width=1500, fig_height=None, fig_dpi=80, **kwargs):
        """ """
        self.fig_width = fig_width / float(fig_dpi)
        # fig_height are handled by matplotlib in inches not pixel
        if fig_height:
            self.fig_height = fig_height / float(fig_dpi)
        else:
            self.fig_height = fig_height
        self.dpi = fig_dpi
        self.features = {}
        self.tracks = []
        #create figure object
        self.fig = matplotlib.pyplot.figure(1, dpi = self.dpi)
        self.ax = self.fig.add_subplot(111)
        self.grid = kwargs.get('grid', 'both') # can be "major", "minor", "both" or None
        self.join_tracks = kwargs.get('join_tracks', True)  

    def add_track(self, track):
        """
        Accepts any object of the  tracks object
        the track should olready contains it's features
        """
        self.tracks.append(track)

    def _draw_tracks(self, **kwargs):
        self.ax.set_axis_off()
        self.Drawn_objects = []
        self.track_axes = []
        Xmin = kwargs.get('xmin', None)
        Xmax = kwargs.get('xmax', None)
        # estimate track height using track.drawn_lines
        # and set the axes positions after they were drawn
        # order and draws tracks
        cbars = False
        self.drawn_lines = 0
        Xs =[]
        for track in self.tracks:
            if track.features:#skip tracks with no features
                track.sort_features(dpi = self.dpi)
                if track.draw_cb:
                    cbars = True
                self.drawn_lines += track.drawn_lines
                Xs.append(track.xmin)
                Xs.append(track.xmax)
        if Xmin == None:
            Xmin = min(Xs)
        if Xmax == None:
            Xmax =max(Xs)
        cbar_extent=0.01
        cbar_axis_space = 0.01
        axis_bottom_pad = 1.0
        if self.join_tracks:
            axix_vertical_pad = 0.0
        else:
            axix_vertical_pad = 1./len(self.tracks)*0.1
        axis_left_pad = 0.025
        figure_bottom_space = 1./(10. + self.drawn_lines)
        axis_width = 1.-2*axis_left_pad
        axis_scale = None# used to persist the same scale on all the tracks
        if cbars:
            axis_width -= (cbar_extent + cbar_axis_space) 
        for track_num, track in enumerate(self.tracks):
            if track.features:#skip tracks with no features
                if axis_scale:
                    axis_height = axis_scale * track.drawn_lines
                else:    
                    axis_height = (float(track.drawn_lines)/self.drawn_lines)  - axix_vertical_pad/2.
                    axis_scale = axis_height / float(track.drawn_lines)
                axis_bottom_pad -= (axis_height + axix_vertical_pad/2.) - figure_bottom_space/len(self.tracks)
                axis = matplotlib.pyplot.axes([axis_left_pad,axis_bottom_pad, axis_width, axis_height ],) 
                self.track_axes.append(axis)
                '''add feature patches to track axes '''
                for feature in track.features:
                    self.Drawn_objects.append(feature)
                    for patch in feature.patches:
                        axis.add_artist(patch)
                        patch.set_transform(axis.transData)# IMPORTANT WORKAROUND!!! if not manually set, transform is not passed correctly in Line2D objects
                    for feat_name in feature.feat_name:
                        axis.add_artist(feat_name)
                if track.draw_cb:
                    cb_axis = matplotlib.pyplot.axes([axis_left_pad + axis_width + cbar_axis_space ,axis_bottom_pad, cbar_extent, axis_height ],) 
                    #cb_axis.axes.set_axis_off()
                    for label in cb_axis.get_yticklabels():
                        label.set_fontsize('xx-small')
                # XXX: track_name unused
                # if track.show_name:
                #    track_name = axis.text(Xmax * 0.01, 1.5, track.name,  horizontalalignment='left', verticalalignment='bottom', fontproperties=track.font_feat,),
                '''handle track axis display and ticks '''
                if track.show_name:
                    axis.set_ylim(track.Ycord, 3)
                else:
                    axis.set_ylim(track.Ycord, 2)
                axis.set_xlim(Xmin, Xmax)
                if self.join_tracks:
                    if (track_num+1 == len(self.tracks)) and ('force no axis' not in track.draw_axis):
                        track.draw_axis.append('bottom')
                    else:
                        if 'bottom' in track.draw_axis:
                            del track.draw_axis[track.draw_axis.index('bottom')]
                    axis.spines["top"].set_color('none')
                else:
                    if (track_num+1 == len(self.tracks)) and ('force no axis' not in track.draw_axis):
                        track.draw_axis.append('bottom')
                for spine in ["right", "left", "top", "bottom"]:
                    if spine not in track.draw_axis:
                        axis.spines[spine].set_color('none')# don't draw axis           
                axis.yaxis.set_ticks([])# dont'show ticks and labels on y
                step=int(round(Xmax/10.,1-len(str(int(Xmax / 10.)))))
                X_major_ticks = xrange(Xmin, Xmax + 1, step)
                step_min = int(round(step / 5.))
                X_minor_ticks = xrange(Xmin, Xmax + 1, step_min)
                if 'top' not in track.draw_axis:
                    axis.xaxis.set_ticks_position('bottom')# only show bottom ticks
                if 'bottom' in track.draw_axis :
                    axis.set_xticks(X_major_ticks)
                    axis.set_xticks(X_minor_ticks,minor=True)
                else:
                    axis.set_xticks([])
                '''ticks labels '''
                if track_num+1 == len(self.tracks): # last track
                    axis.set_xticklabels(map(str,X_major_ticks), fontsize='x-small')
                else:
                    axis.yaxis.set_ticks([])
                    axis.set_xticklabels([])    
                '''draw grid'''
                if self.grid:
                    if (self.grid == 'major') or (self.grid == 'both'):
                        for X in X_major_ticks:
                            axis.axvline(X,ls=':',c='grey',alpha=0.66, zorder = -1)
                    if (self.grid == 'minor') or (self.grid == 'both'):
                        for X in X_minor_ticks:
                            axis.axvline(X,ls=':',c='grey',alpha=0.33, zorder = -1)
        '''set panel size and panning '''
        if not self.fig_height:#automatcally set fig height basing on the number of features
            self.fig_height=self.drawn_lines*self.fig_width/50
            self.fig.set_figheight(self.fig_height)
        else:
            self.fig.set_figheight(self.fig_height)
        self.fig.set_figwidth(self.fig_width)

    def boxes(self):
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
                    bbox=patch.get_window_extent(None)
                    xs_patches.append(bbox.xmax)
                    xs_patches.append(bbox.xmin)
                    ys_patches.append(bbox.ymax)
                    ys_patches.append(bbox.ymin)
            xmin, ymin = trans.transform([min(xs_patches),min(ys_patches)])
            xmax, ymax = trans.transform([max(xs_patches),max(ys_patches)])
            left = xmin
            top = img_height-ymin
            right = xmax
            bottom = img_height-ymax
            yield dict(feature=obj, left=left, top=top, right=right, bottom=bottom, track=None)

    def _create_html_map(self, map_name = 'biograpy-map', map_id = 'biograpy-map', new_win = False, **kwargs):
        """
        returns the corresponding html map from self.Drawn_objects in self.htmlmap

        new_win ---> set target="_blank" on area links
        """
        areas =[]
        if new_win:
            target = '_blank'
        else:
            target = '_self'
        for box in self.boxes():
            obj = box['feature']
            area_dict = dict(shape = 'rect', # shape: rect, circle, poly
                             coords = '%(left)i,%(top)i,%(right)i,%(bottom)i' % box,
                             href = obj.url or '#%s'%obj.name, #href
                             target = target,
                             onmouseover = obj.onmouseover,
                             alt = obj.name )#text to be displayed in pop up span?
            area_html = '''<area shape="%(shape)s" coords="%(coords)s" href="%(href)s" target="%(target)s" alt="%(alt)s">''' % area_dict
            areas.append(area_html)
        self.htmlmap ='''<map name="%s" id="%s">\n %s \n</map>'''%(map_name, map_id, '\n'.join(areas))
        return

    def imagemap(self, new_win=True, **kwargs):
        self._create_html_map(new_win = new_win, **kwargs)
        return self.htmlmap

    def save(self, output, new_win=True, **kwargs):
        '''
        output ---> a string containing the file path or a file-like handler.
        format ---> a format MUST be specified if a file handler is passed
                    supported formats are: emf, eps, pdf, png, ps, raw, rgba, svg, svgz
        '''
        self._draw_tracks()
        self._create_html_map(new_win = new_win, **kwargs)
        matplotlib.pyplot.savefig(output, dpi=self.fig.get_dpi(), **kwargs)

