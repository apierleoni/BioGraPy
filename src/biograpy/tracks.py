
'''

track = BaseTrack ()

track.add_feature(GraphicFeature)



'''
import operator
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
    
    base class for all tracs types'''
    
    Ycord = 0.0
    _betw_feat_space = 2
    
    default_cm='jet'#deafult matplotlib colormap to use
    
    def __init__(self, *args, **kwargs):
        '''track name will be displayed,
        show_name = 'up', 'down', None
        '''
        
        self.name = kwargs.get('name','')
        if not self.name:
            self.show_name = None
        else:
            self.show_name = kwargs.get('show_name', 'up')#'up', 'down', None
        self.max_score = kwargs.get('max_score', None)
        self.min_score = kwargs.get('min_score', None)
        self.draw_axis = kwargs.get('draw_axis', ['bottom'] )# list containing witch spine to draw for the track. available spines are: ["right", "left", "top", "bottom", "force no axis"]
        
        self.norm = kwargs.get('norm', None)# normalizing function, if None build one. could be any function taking a value an returning a float between 0 and 1
        if not self.norm:
                self.norm = colors.normalize(vmin = self.min_score, vmax = self.max_score)
                
        
        self.cm = cm.get_cmap(kwargs.get('cm', self.default_cm))
        #self.colormap = cm.get_cmap(self.cm)
        self.draw_cb = kwargs.get('draw_cb', False)
        
        
        font_feat=FontProperties()
        font_feat.set_size('small')
        font_feat.set_family('serif')
        font_feat.set_weight('semibold')
        self.font_feat = font_feat
        
        self.drawn_lines = 0
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
    

      
    def _collapse(self, dpi):

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
                        bbox=patch.get_window_extent(None,)
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
                            else:
                                current_y = patch.get_y()
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
        
    def _draw_ordered_features(self, feat_list = None):
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
                else:
                    current_y = patch.get_y()
                    patch.set_y(current_y + self.Ycord)
            for iname, fname in enumerate(feat2draw.feat_name):
                y=fname.get_position()[1]
                feat2draw.feat_name[iname].set_y(y + self.Ycord)
                current_x, current_y = fname.xytext
                feat2draw.feat_name[iname].xytext = (current_x, current_y + self.Ycord)
            self.Ycord-=self._betw_feat_space
            self.drawn_lines += 1

    def _draw_features(self):
        for feat_numb, feat2draw in enumerate(self.features):
            if feat2draw.color_by_cm:
                if feat2draw.use_score_for_color:
                    feat2draw.cm_value = feat2draw.score
                else:# color by feature number
                    if not feat2draw.cm_value:
                        self.norm = colors.normalize(1,len(self.features)+1,)
                        feat2draw.cm_value = feat_numb +1
                feat2draw.fc = self.cm(self.norm(feat2draw.cm_value))
            feat2draw.draw_feature()
            feat2draw.draw_feat_name()

            
    def sort_features(self, dpi = 80):
        self._draw_features()
        if self.sort_by =='collapse':
            self._collapse(dpi)
        elif self.sort_by =='score':
            feat_list = self._order_by_score()
            self._draw_ordered_features(feat_list)
        elif self.sort_by =='length':
            feat_list = self._order_by_length()
            self._draw_ordered_features(feat_list)
        else:
            self._draw_ordered_features()
    
        
        
class ReferenceTrack(BaseTrack):
    '''
::
    
    Track to be used for reference objects'''
    
    def __init__(self, **kwargs):
        BaseTrack.__init__(self, **kwargs)
        if self.name:
            self.show_name = None
            
            
            
class PlotTrack(BaseTrack):
    '''
::    
    
    Track to be used for plotting  PlotFeatures
    multiple PlotFeatures in the same track shares the same axes'''
    
    def __init__(self, **kwargs):
        BaseTrack.__init__(self, **kwargs)
        if self.name:
            self.show_name = 'down'