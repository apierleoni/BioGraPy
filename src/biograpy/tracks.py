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
from matplotlib.patches import Rectangle, Circle, Wedge, Polygon, FancyBboxPatch, FancyArrow
from matplotlib.font_manager import FontProperties
from matplotlib.text import Annotation




class BaseTrack(object):
    ''' 
::
    
    base class for all tracs types
    
    track = BaseTrack ()
    track.add_feature(GraphicFeature)
'''
    
    Ycord = 0.0
    _betw_feat_space = 2
    
    default_cm='Accent'#deafult matplotlib colormap to use
    
    def __init__(self, *args, **kwargs):
        '''track name will be displayed,
        show_name = 'up', 'down', None
        '''
        
        self.name = kwargs.get('name','')
        if self.name:
            self.show_name = kwargs.get('show_name', 'top')#'top', 'bottom', None
        else:
            self.show_name = None
        self.max_score = kwargs.get('max_score', None)
        self.min_score = kwargs.get('min_score', None)
        self.draw_axis = kwargs.get('draw_axis', ["bottom"] )# list containing witch spine to draw for the track. available spines are: ["right", "left", "top", "bottom", "force no axis"] 'force no axis' option will force the bottom axis not to be shown even itf it is the last track. only needed in when tracks are requested to be joined
        self.xticks_major = kwargs.get('xticks_major', None )
        self.xticks_minor = kwargs.get('xticks_minor', None )
        self.xticklabels_major = kwargs.get('xticklabels_major', None )
        self.xticklabels_minor = kwargs.get('xticklabels_minor', None )
        self.yticks_major = kwargs.get('yticks_major', None )
        self.yticks_minor = kwargs.get('yticks_minor', None )
        self.yticklabels_major = kwargs.get('yticklabels_major', None )
        self.yticklabels_minor = kwargs.get('yticklabels_minor', None )
        self.x_use_sequence = kwargs.get('x_use_sequence', False )# if a sequence is passed it will be used to display ticklabels on X axis. must be an iterable
        self.show_xticklabels = kwargs.get('show_xticklabels', False )# force showing ticklabels on tracks different from the bottom one
        self.show_yticklabels = kwargs.get('show_yticklabels', False )
        self.tickfontsize =  kwargs.get('tickfontsize', 'x-small' )
        self.tickfontsize =  kwargs.get('tickfontsize_major', 'x-small' )
        self.tickfontsize_minor =  kwargs.get('tickfontsize_minor', 'xx-small' )
        self.ymax = 0
        self.track_height = kwargs.get('track_height', 0)# used to specify the the axis heigth by pixel number
        
        self.norm = kwargs.get('norm', None)# normalizing function, if None build one. could be any function taking a value an returning a float between 0 and 1
        if not self.norm:
            self.norm = colors.normalize(vmin = self.min_score, vmax = self.max_score)
                
        
        self.cm = cm.get_cmap(kwargs.get('cm', self.default_cm))
        #self.colormap = cm.get_cmap(self.cm)
        self.draw_cb = kwargs.get('draw_cb', False)
        self.cb_alpha = kwargs.get('cb_alpha', 1)
        self.cb_label = kwargs.get('cb_label', None)
        
        
        self.name_font_size = kwargs.get('name_font_size', 'small' )
        self.name_font_family = kwargs.get('name_font_size', 'serif' )
        self.name_font_weight = kwargs.get('name_font_weight', 'semibold' )
        font_feat=FontProperties()
        font_feat.set_size(self.name_font_size)
        font_feat.set_family(self.name_font_family)
        font_feat.set_weight('semibold')
        self.name_font_feat = font_feat
        
        self.drawn_lines = kwargs.get('track_lines', 0 )#number of features to be counted do determine track height
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
        #check for graphicfeature istance?
        self.add_feature(feature)
        
    def extend(self, features):
        #check for graphicfeature istance?
        for feature in features:
            self.add_feature(feature)
    

      
    def _collapse(self, dpi,):

        plt.draw()
        # this is necessary to estimate the text dimensions and avoid collisions until draw
        # is not called the text object do not know its renderer and cannot call get_window_extent method.
        line_controller=[]
        size_memory={}
        draw_features = []
        while len(draw_features) < len(self.features):
            for feat_numb, feat2draw in enumerate(self.features):
                '''estimate feature lenght'''
                if feat_numb not in size_memory:
                    xs_patches=[]
                    for patch in feat2draw.patches:
                        #for p in patch:
                        try:
                            bbox=patch.get_window_extent(None,)
                        except:
                            warnings.warn('could not find box coordinated for patch: '+str(patch) )
                            continue
                        xs_patches.append(bbox.xmax)
                        xs_patches.append(bbox.xmin)
                    for fname in feat2draw.feat_name:
                        # set the correct dpi to correctly estimate text size.
                        # required by Text class in matplotlib
                        bbox = fname.get_window_extent(dpi = dpi)
                        xs_patches.append(bbox.xmax)
                        xs_patches.append(bbox.xmin)
                    size_memory[feat_numb]={'left_margin' : min(xs_patches),
                                    'right_margin' : max(xs_patches)}

                ''' Check for collisions both on text and patches'''
                if feat_numb not in draw_features:
                    draw=True
                    for prev_start,prev_end in line_controller:
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
                                except AttributeError:
                                    current_y = patch.get_position()[1]
                                patch.set_y(current_y + self.Ycord)
                            
                        for iname, fname in enumerate(feat2draw.feat_name):
                            y=fname.get_position()[1]
                            feat2draw.feat_name[iname].set_y(y + self.Ycord)
                            current_x, current_y = fname.xytext
                            feat2draw.feat_name[iname].xytext = (current_x, current_y + self.Ycord)
                        draw_features.append(feat_numb)
                        
            #if len(draw_features) < len(self.features):
            self.Ycord-=self._betw_feat_space
            line_controller=[]
            self.drawn_lines += 1

      
    def _order_by_score(self,):
        
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
        '''this order basing on the actual length of the feature patches in the panel '''
        
        
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
        
    def _draw_ordered_features(self, feat_list = None,):
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
                y=fname.get_position()[1]
                feat2draw.feat_name[iname].set_y(y + self.Ycord)
                current_x, current_y = fname.xytext
                feat2draw.feat_name[iname].xytext = (current_x , current_y + self.Ycord)
            self.Ycord-=self._betw_feat_space
            self.drawn_lines += 1

    def _draw_features(self, **kwargs):
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
                        self.norm = colors.normalize(1,len(self.features)+1,)
                        feat2draw.cm_value = feat_numb +1
                    feat2draw.fc = self.cm(self.norm(feat2draw.cm_value))
            feat2draw.draw_feature()
            feat2draw.draw_feat_name(xoffset = xoffset)

            
    def sort_features(self, dpi = 80, **kwargs):
        
        self._draw_features(**kwargs)
        if self.sort_by =='collapse':
            self._collapse(dpi, )
        elif self.sort_by =='score':
            feat_list = self._order_by_score()
            self._draw_ordered_features(feat_list,)
        elif self.sort_by =='length':
            feat_list = self._order_by_length()
            self._draw_ordered_features(feat_list, )
        else:
            self._draw_ordered_features()
    
        
        
 
            
class PlotTrack(BaseTrack):
    '''
::    
    
    Track to be used for plotting  PlotFeatures
    multiple PlotFeatures in the same track shares the same axes'''
    
    def __init__(self,*args, **kwargs):
        BaseTrack.__init__(self, *args, **kwargs)

        self.Ycord = self.ymin = kwargs.get('ymin', -1)
        self.draw_axis = kwargs.get('draw_axis', ['left', 'bottom'] )
        self.ymax = kwargs.get('ymax', 1)
        self.show_yticklabels = kwargs.get('show_yticklabels', True )# force showing ticklabels on y axis, default is True. will only be showed if ticks are defined
        self.drawn_lines = kwargs.get('track_lines', 4 )#abs(self.ymax - self.ymin)  number of features to be counted do determine track height
        
                                
    def _collapse(self, dpi):
        plt.draw()
        

      
    def _order_by_score(self,):
        return self.features


    
    def _order_by_length(self,):
        return self.features
        
                
    def _draw_ordered_features(self, feat_list = None):
        return

    def _draw_features(self):
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
            feat2draw.draw_feat_name()

            
    def sort_features(self, dpi = 80, **kwargs):
        self._draw_features(**kwargs)
                    