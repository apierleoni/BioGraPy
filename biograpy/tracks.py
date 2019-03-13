'''
Created on 28/ott/2009

@author: andrea pierleoni

Handles a set of Graphic Representation for biological Entities.

quote:
That's what I like in matplotlib: no matter how hard you try,
there's always a simpler solution you're not yet aware of...

'''
import operator, warnings
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle, Wedge, Polygon, FancyBboxPatch, FancyArrow, FancyArrowPatch
from matplotlib.font_manager import FontProperties
from matplotlib.text import Annotation
from matplotlib.text import Text
from .features import TextSequence, PrettyTextSequence



def find_renderer(fig):

    if hasattr(fig.canvas, "get_renderer"):
        #Some backends, such as TkAgg, have the get_renderer method, which 
        #makes this easy.
        renderer = fig.canvas.get_renderer()
    else:
        #Other backends do not have the get_renderer method, so we have a work 
        #around to find the renderer.  Print the figure to a temporary file 
        #object, and then grab the renderer that was used.
        #(I stole this trick from the matplotlib backend_bases.py 
        #print_figure() method.)
        import io
        fig.canvas.print_pdf(io.BytesIO())
        renderer = fig._cachedRenderer
    return(renderer)



class BaseTrack(object):
    '''
    The :class:`BaseTrack` acts as a collector of :class:`~biograpy.features` \
    objects, and accounts for their ordering and color.
    The :class:`BaseTrack` is a functiona baseclass to use in order to generate\ 
    additional tracks.
    Track object have to be added to a :class:`~biograpy.drawer.Panel` instance\
    once they are filled with features


    ``track = BaseTrack(name = 'test')``

    gives a empty track with default values named ``'test'``

    Valid keyword arguments for :class:`BaseTrack` are:

        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        name                  Define the track `name`
        show_name             ``'top'`` | ``'bottom'`` |  ``None`` , default is
                              ``'top'`` . If a track `name` is specified  set 
                              where the name wil be displayed
        name_font_size        default ``'small'``, accepts font size parameter 
                              as in matplotlib
        name_font_family      default ``'serif'``, accepts font family 
                              parameter as in matplotlib
        name_font_weight      default ``'semibold'``, accepts font weight 
                              parameter as in matplotlib
        draw_axis             a list contanining this possible values: 
                              ``'right'`` | ``'left'`` | ``'top'`` | 
                              ``'bottom'`` | ``'force no axis'`` . default is 
                              ``['bottom']``. each specified axis position will
                              be drawn. if ``'force no axis'`` is specified no 
                              axis will be shown, this is only needed if tracks
                              are requested to be joined.
        x_use_sequence        if a string iterable is passed it will be used as
                              X tick labels  
        xticks_major          set X major ticks as in matplotlib
        xticks_minor          set X minor ticks as in matplotlib
        xticklabels_major     set X major tick labels as in matplotlib      
        xticklabels_minor     set X minor tick labels as in matplotlib                
        yticks_major          set Y major ticks as in matplotlib
        yticks_minor          set Y minor ticks as in matplotlib  
        yticklabels_major     set Y major tick labels as in matplotlib          
        yticklabels_minor     set Y minor tick labels as in matplotlib                
        show_xticklabels      ``True`` | ``False``. default is ``True``
        show_yticklabels      ``True`` | ``False``. default is ``False``
        tickfontsize_major    default ``'x-small'``, accepts font size parameter 
                              as in matplotlib
        tickfontsize_minor    default ``'xx-small'``, accepts font size 
                              parameter as in matplotlib    
        track_height          specify the the axis heigth by pixel number. 
                              default is ``0`` and will be set automatically        
        max_score             default is ``None``, if not specified will be 
                              dynamically computed basing on feature scores
        min_score             default is ``None``, if not specified will be 
                              dynamically computed basing on feature scores
        norm                  score normalizing function. default is ``None`` 
                              and will be created using maximum and minimum 
                              feature scores. could be any function taking a 
                              int or float value an returning a float between 
                              ``0.`` and ``1.``
        cm                    matplotlib colormap name to use for score-based 
                              feature color. default is ``'Accent'``
        draw_cb               ``True`` | ``False``. default is ``False``. draws
                              a colorbar inside the track to explain score-based
                              feature color
        cb_alpha              colorbar alpha value, a float between ``0.`` and 
                              ``1.``
        cb_label              colorbar label, a string
        track_lines           internal parameter used to defin track heigth in 
                              automatic mode. default is 
                              ``1*number_of feature_lines``. if forced can 
                              results in wrong feature proportions
        sort_by               ``'score'`` | ``'length'`` | ``'collapse'`` | 
                              ``None``. 
                              default is ``'collapse'`` 
                              
                              * ``'collapse'`` will arrange features in the \
                              minimum vertical space possible. 
                              * ``'score'`` will arrange features basing \
                              on score values.  
                              * ``'length'`` will arrange features basing \
                              on feature length.
                              * ``None`` will arrange features basing \
                              in the order they are passed to the track.
        sort_order            ``'top'`` | ``'bottom'``, default is ``'top'``.
                              use to reverse order
        ===================== ==================================================
'''
    
    Ycord = 0.0
    _betw_feat_space = 2.5
    
    default_cm='Accent'#deafult matplotlib colormap to use
    
    def __init__(self, *args, **kwargs):
        '''
        '''
        
        self.name = kwargs.get('name','')
        if self.name:
            self.show_name = kwargs.get('show_name', 'top')
        else:
            self.show_name = None
        self.fig = kwargs.get('fig', None)
        self.max_score = kwargs.get('max_score', None)
        self.min_score = kwargs.get('min_score', None)
        self.draw_axis = kwargs.get('draw_axis', ["bottom"] )
        self.xticks_major = kwargs.get('xticks_major', None )
        self.xticks_minor = kwargs.get('xticks_minor', None )
        self.xticklabels_major = kwargs.get('xticklabels_major', None )
        self.xticklabels_minor = kwargs.get('xticklabels_minor', None )
        self.yticks_major = kwargs.get('yticks_major', None )
        self.yticks_minor = kwargs.get('yticks_minor', None )
        self.yticklabels_major = kwargs.get('yticklabels_major', None )
        self.yticklabels_minor = kwargs.get('yticklabels_minor', None )
        self.x_use_sequence = kwargs.get('x_use_sequence', False )
        self.show_xticklabels = kwargs.get('show_xticklabels', False )
        self.show_yticklabels = kwargs.get('show_yticklabels', False )
        self.tickfontsize =  kwargs.get('tickfontsize', 'x-small' )
        self.tickfontsize =  kwargs.get('tickfontsize_major', 'x-small' )
        self.tickfontsize_minor =  kwargs.get('tickfontsize_minor', 'xx-small' )
        self.ymax = 0
        self.track_height = kwargs.get('track_height', 0)
        
        self.norm = kwargs.get('norm', None)# normalizing function, if None build one. could be any function taking a value an returning a float between 0 and 1
        if not self.norm:
            self.norm = colors.Normalize(vmin = self.min_score, vmax = self.max_score)
                
        
        self.cm = cm.get_cmap(kwargs.get('cm', self.default_cm))
        #self.colormap = cm.get_cmap(self.cm)
        self.draw_cb = kwargs.get('draw_cb', False)
        self.cb_alpha = kwargs.get('cb_alpha', 1)
        self.cb_label = kwargs.get('cb_label', None)
        
        
        self.name_font_size = kwargs.get('name_font_size', 'small' )
        self.name_font_family = kwargs.get('name_font_family', 'serif' )
        self.name_font_weight = kwargs.get('name_font_weight', 'semibold' )
        font_feat=FontProperties()
        font_feat.set_size(self.name_font_size)
        font_feat.set_family(self.name_font_family)
        font_feat.set_weight('semibold')
        self.name_font_feat = font_feat
        
        self.drawn_lines = kwargs.get('track_lines', 0 )
        if self.show_name:
            self.drawn_lines += 1
            
        self.sort_by = kwargs.get('sort_by', 'collapse' )# can also be: score, length, collapse or None. None means the features are plotted in the order they get added
        self.sort_order = kwargs.get('sort_order', 'top')# can be top or bottom
        self.features = [] # this will contains all the Graphicfeatures of the panel
        self.xmin = None
        self.xmax = None
        
        for feature in args:
            self.add_feature(feature)
        
    def add_feature(self, feature):
        '''add  :class:`~biograpy.features` object to track.'''
        def add(feature):
            self.features.append(feature)
            if self.xmin == None:
                self.xmin = feature.start
            else:
                self.xmin = min([feature.start, self.xmin])
            if self.xmax == None:
                self.xmax = feature.end
            else:
                self.xmax = max([feature.end, self.xmax])
        
        if isinstance(feature, list):
            for feat in feature:
                add(feature)
        else:
            add(feature)
            
    def append(self, feature):
        '''add  :class:`~biograpy.features` object to track.
        same as :func:`~biograpy.tracks.BaseTrack.add_feature`'''
        #check for graphicfeature istance?
        self.add_feature(feature)
        
    def extend(self, features):
        '''add a sequence of features to track'''
        #check for graphicfeature istance?
        for feature in features:
            self.add_feature(feature)
    

      
    def _collapse(self, dpi,):
        '''collapse features '''
        plt.draw()
        renderer = find_renderer(plt.gcf())
        # this is necessary to estimate the text dimensions and avoid collisions until draw
        # is not called the text object do not know its renderer and cannot call get_window_extent method.
        line_controller = []
        size_memory = {}
        draw_features = []
        last_feature_xlim = [0,0]
        while len(draw_features) < len(self.features):
            for feat_numb, feat2draw in enumerate(self.features):
                '''estimate feature lenght'''
                if feat_numb not in size_memory:
                    xs_patches = []
                    if type(feat2draw) in [TextSequence, PrettyTextSequence, Text]:
                        xs_patches.append(feat2draw.start)
                        xs_patches.append(feat2draw.end)
                    else:
                        for patch in feat2draw.patches:
                            # print("patch", patch, type(patch))
                            # print(patch.get_figure().number)
                            bbox = patch.get_window_extent(renderer=renderer)
                            # print("bbox", bbox, "bbox.xmax", bbox.xmax, "bbox.xmin", bbox.xmin)
                            xs_patches.append(bbox.xmax)
                            xs_patches.append(bbox.xmin)
                    for fname in feat2draw.feat_name:
                        # set the correct dpi to correctly estimate text size.
                        # required by Text class in matplotlib
                        # This can be a problem with ax.annotate, bounding box is too wide?
                        # We do not take into consideration text annotations
                        #bbox = fname.get_window_extent(dpi = dpi)
                        #bbox = fname.get_window_extent()
                        #xs_patches.append(bbox.xmax)
                        #xs_patches.append(bbox.xmin)
                        pass
                    size_memory[feat_numb] = {'left_margin':min(xs_patches),
                                              'right_margin':max(xs_patches)}

                ''' Check for collisions both on text and patches'''
                if feat_numb not in draw_features:
                    draw = True
                    for prev_start, prev_end in line_controller:
                        # print("line_controller",prev_start, prev_end, size_memory[feat_numb])
                        if  (prev_start <= size_memory[feat_numb]['left_margin'] <= prev_end) or \
                            (prev_start <= size_memory[feat_numb]['right_margin'] <= prev_end) or \
                            ((size_memory[feat_numb]['left_margin'] < prev_start < size_memory[feat_numb]['right_margin']) and \
                             (size_memory[feat_numb]['left_margin'] < prev_end < size_memory[feat_numb]['right_margin'])):
                            draw = False
                            break
                    if draw:
                        '''Draw if not collision '''
                        line_controller.append([size_memory[feat_numb]['left_margin'],size_memory[feat_numb]['right_margin']])
                        for patch in feat2draw.patches:
                            if isinstance(patch, Line2D):
                                current_ys = patch.get_ydata()
                                new_ys = map(operator.add, current_ys, [self.Ycord] * len(current_ys))
                                patch.set_ydata(new_ys)

                            elif isinstance(patch, FancyArrow):
                                current_xy=patch.get_xy()
                                new_xy=[]
                                for x, y in current_xy:
                                    new_xy.append([x, y + self.Ycord])
                                patch.set_xy(new_xy)
                            elif isinstance(patch, Annotation):
                                current_x, current_y = patch.xytext
                                patch.xytext = (current_x, current_y + self.Ycord)
                            else:
                                try:
                                    current_y = patch.get_y()
                                    patch.set_y(current_y + self.Ycord)
                                except AttributeError:
                                    if type(patch) is FancyArrowPatch:
                                        arrow_pos = [list(pos) for pos in patch._posA_posB]
                                        arrow_pos[0][1] = arrow_pos[0][1] + self.Ycord
                                        arrow_pos[1][1] = arrow_pos[1][1] + self.Ycord
                                        patch.set_positions(arrow_pos[0], arrow_pos[1])
                                    else:
                                        current_y = patch.get_position()[1]
                                        patch.set_y(current_y + self.Ycord)
                            
                        for iname, fname in enumerate(feat2draw.feat_name):
                            y=fname.get_position()[1]
                            feat2draw.feat_name[iname].set_y(y + self.Ycord)
                            feat2draw.Y2 = (y + self.Ycord)
                            #current_x, current_y = fname.xyann
                            #feat2draw.feat_name[iname].xytext = (current_x, current_y + self.Ycord)
                            #print("feat_name ",fname,"y",y,"new y",fname.get_position()[1])
                        draw_features.append(feat_numb)
                        
            #if len(draw_features) < len(self.features):
            
            self.Ycord-=self._betw_feat_space
            line_controller=[]
            self.drawn_lines += 1

      
    def _order_by_score(self,):
        '''order features by score '''
        feat_list = []
        for feat in self.features:
            feat_list.append([feat.score, feat])
        if self.sort_order == 'top':
            feat_list.sort()
        elif self.sort_order == 'bottom':
            feat_list.sort()
            feat_list.reverse()
        else:
            raise ValueError('Wrong Sort order option: %s'%self.sort_order )
        
        return [feat for (i,feat) in feat_list]

    
    def _order_by_length(self,):
        ''' order basing on the actual length of the feature patches in the figure '''
        
        feat_list = []
        for feat in self.features:
            xs_patches = []
            for patch in feat.patches:
                bbox=patch.get_window_extent(None)
                xs_patches.append(bbox.xmax)
                xs_patches.append(bbox.xmin)
            feat_list.append([max(xs_patches) - min(xs_patches), feat])
        if self.sort_order == 'top':
            feat_list.sort()
        elif self.sort_order == 'bottom':
            feat_list.sort()
            feat_list.reverse()
        else:
            raise ValueError('Wrong Sort order option: %s'%self.sort_order )
        
        return [feat for (i,feat) in feat_list]
        
    def _draw_ordered_features(self, feat_list=None,):
        '''draws one feature per line in the track in the order they are passed'''
        if not feat_list:
            feat_list = self.features
        for feat2draw in feat_list:
            for patch in feat2draw.patches:
                if isinstance(patch, Line2D):
                    current_ys = patch.get_ydata()
                    new_ys = map(operator.add, current_ys, [self.Ycord] * len(current_ys))
                    patch.set_ydata(new_ys)
                elif isinstance(patch, FancyArrow):
                    current_xy=patch.get_xy()
                    new_xy=[]
                    for x, y in current_xy:
                        new_xy.append([x, y + self.Ycord])
                    patch.set_xy(new_xy)
                elif isinstance(patch, Annotation):
                    current_x, current_y = patch.xytext
                    patch.xytext = (current_x, current_y + self.Ycord)
                else:
                    try: 
                        current_y = patch.get_y()
                    except AttributeError:
                        current_y = patch.get_position()[1]
                    patch.set_y(current_y + self.Ycord)
            for iname, fname in enumerate(feat2draw.feat_name):
                y = fname.get_position()[1]
                feat2draw.feat_name[iname].set_y(y + self.Ycord)
                #current_x, current_y = fname.xytext
                current_x, current_y = fname.xyann
                feat2draw.feat_name[iname].xytext = (current_x , current_y + self.Ycord)
            self.Ycord -= self._betw_feat_space
            self.drawn_lines += 1

    def _draw_features(self, fig=None, **kwargs):
        '''draw features '''
        xoffset = kwargs.get('xoffset',0)
        for feat_numb, feat2draw in enumerate(self.features):
            if feat2draw.color_by_cm:
                if feat2draw.use_score_for_color:
                    feat2draw.cm_value = feat2draw.score
                    feat2draw.fc = self.cm(feat2draw.cm_value)
                    if not feat2draw.ec:
                        feat2draw.ec = feat2draw.fc
                else:# color by feature number
                    if not feat2draw.cm_value:
                        self.norm = colors.Normalize(1,len(self.features)+1,)
                        feat2draw.cm_value = feat_numb +1
                    feat2draw.fc = self.cm(self.norm(feat2draw.cm_value))
            if False:
            # if type(feat2draw) in [PrettyTextSequence, TextSequence]:
                feat2draw.draw_feature(fig=fig)
            else:
                feat2draw.draw_feature()
            feat2draw.draw_feat_name(xoffset=xoffset)

            
    def _sort_features(self, dpi=80, **kwargs):
        ''' sort features basing on the chosen mode '''
        self._draw_features(fig=plt.gcf(), **kwargs)
        if self.sort_by =='collapse':
            self._collapse(dpi)
        elif self.sort_by =='score':
            feat_list = self._order_by_score()
            self._draw_ordered_features(feat_list,)
        elif self.sort_by =='length':
            feat_list = self._order_by_length()
            self._draw_ordered_features(feat_list, )
        else:
            self._draw_ordered_features()
            pass
    
        
        
 
            
class PlotTrack(BaseTrack):
    '''

    The :class:`PlotTrack` acts as a collector of :class:`~biograpy.features` \
    objects that encodes for plots. 
    Track object have to be added to a :class:`~biograpy.drawer.Panel` instance\
    once they are filled with features


    ``track = PlotTrack(name = 'test', ymin= 0, ymax = 10)``

    gives a empty plot track  named ``'test'`` with y axis ranging from 0 to 10.

    Valid keyword arguments for :class:`PlotTrack` are:

        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        name                  Define the track `name`
        show_name             ``'top'`` | ``'bottom'`` |  ``None`` , default is
                              ``'top'`` . If a track `name` is specified  set 
                              where the name wil be displayed
        name_font_size        default ``'small'``, accepts font size parameter 
                              as in matplotlib
        name_font_family      default ``'serif'``, accepts font family 
                              parameter as in matplotlib
        name_font_weight      default ``'semibold'``, accepts font weight 
                              parameter as in matplotlib
        draw_axis             a list contanining this possible values: 
                              ``'right'`` | ``'left'`` | ``'top'`` | 
                              ``'bottom'`` | ``'force no axis'`` . default is 
                              ``['left', 'bottom']``. each specified axis 
                              position will be drawn. 
                              if ``'force no axis'`` is specified no 
                              axis will be shown, this is only needed if tracks
                              are requested to be joined.
        ymin                  default is ``-1``. define minimum Y value to plot.
        ymax                  default is ``1``. define maximum Y value to plot.
        x_use_sequence        if a string iterable is passed it will be used as
                              X tick labels  
        xticks_major          set X major ticks as in matplotlib
        xticks_minor          set X minor ticks as in matplotlib
        xticklabels_major     set X major tick labels as in matplotlib      
        xticklabels_minor     set X minor tick labels as in matplotlib                
        yticks_major          set Y major ticks as in matplotlib
        yticks_minor          set Y minor ticks as in matplotlib  
        yticklabels_major     set Y major tick labels as in matplotlib          
        yticklabels_minor     set Y minor tick labels as in matplotlib                
        show_xticklabels      ``True`` | ``False``. default is ``True``
        show_yticklabels      ``True`` | ``False``. default is ``True``
        tickfontsize_major    default ``'x-small'``, accepts font size parameter 
                              as in matplotlib
        tickfontsize_minor    default ``'xx-small'``, accepts font size 
                              parameter as in matplotlib    
        track_height          specify the the axis heigth by pixel number. 
                              default is ``0`` and will be set automatically        
        max_score             default is ``None``, if not specified will be 
                              dynamically computed basing on feature scores
        min_score             default is ``None``, if not specified will be 
                              dynamically computed basing on feature scores
        norm                  score normalizing function. default is ``None`` 
                              and will be created using maximum and minimum 
                              feature scores. could be any function taking a 
                              int or float value and returning a float between 
                              ``0.`` and ``1.``
        cm                    matplotlib colormap name to use for score-based 
                              feature color. default is ``'Accent'``
        draw_cb               ``True`` | ``False``. default is ``False``. draws
                              a colorbar inside the track to explain score-based
                              feature color
        cb_alpha              colorbar alpha value, a float between ``0.`` and 
                              ``1.``
        cb_label              colorbar label, a string
        track_lines           default is ``4``. increase or decrease to make 
                              the track bigger or smaller
        ===================== ==================================================
    
    '''
    
    def __init__(self,*args, **kwargs):
        BaseTrack.__init__(self, *args, **kwargs)

        self.Ycord = self.ymin = kwargs.get('ymin', -1)
        self.draw_axis = kwargs.get('draw_axis', ['left', 'bottom'] )
        self.ymax = kwargs.get('ymax', 1)
        self.show_yticklabels = kwargs.get('show_yticklabels', True )
        self.drawn_lines = kwargs.get('track_lines', 4 )#  number of features to be counted do determine track height
        
                                
    def _collapse(self, dpi):
        plt.draw()
        

      
    def _order_by_score(self,):
        return self.features


    
    def _order_by_length(self,):
        return self.features
        
                
    def _draw_ordered_features(self, feat_list = None):
        return

    def _draw_features(self, **kwargs):
        xoffset = kwargs.get('xoffset',0)
        for feat_numb, feat2draw in enumerate(self.features):
            if feat2draw.color_by_cm:
                if feat2draw.use_score_for_color:
                    feat2draw.cm_value = feat2draw.score
                    feat2draw.fc = self.cm(feat2draw.cm_value)
                else:# color by feature number
                    if not feat2draw.cm_value:
                        self.norm = colors.normalize(1,len(self.features)+1,)
                        feat2draw.cm_value = feat_numb +1
                    feat2draw.fc = self.cm(self.norm(feat2draw.cm_value))
            feat2draw.draw_feature()
            feat2draw.draw_feat_name(xoffset = xoffset)

            
    def _sort_features(self, dpi = 80, **kwargs):
        self._draw_features(**kwargs)
                    
