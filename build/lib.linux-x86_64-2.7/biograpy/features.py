'''
Created on 28/ott/2009

@author: andrea pierleoni

Handles a Graphic Representation for biological Entities.
Can read info from BioPython SeqFeature object

quote:
That's what I like in matplotlib: no matter how hard you try,
there's always a simpler solution you're not yet aware of...
'''

import operator
import random
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import Subplot
from matplotlib.font_manager import FontProperties
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.patches import Rectangle,Circle, Wedge, Polygon,FancyBboxPatch, FancyArrow, FancyArrowPatch, ArrowStyle
from matplotlib.colors import Normalize,LogNorm,ListedColormap,LinearSegmentedColormap
from Bio.SeqFeature import ExactPosition, CompoundLocation


defaultFontSize = 11
defaultFontSizeScalingFactor = 110    # Works well for monospace font
textFeatureHeight = 1.1
arrowHeadLengthRelative = 0.01



class BaseGraphicFeature(object):
    ''' 

    Base class for all `GraphicFeature` types and is not functional.
    This class template misses the `draw_feature()` method that has to be
    implemented in the full functional classes.
    To create a new functional `GraphicFeature` class this conditions 
    must me met:

    * `self.draw_feature()` must return one or a list of matplotlib Artists in \
      `self.patches`
    * `self.start` must return the lower extent of the feature in the \ 
      coordinate system
    * `self.end` must return the higher extent of the feature in the \
      coordinate system
    * `self.__init__` method must be expanded accordingly to the new input
    * `self.draw_feat_name()` method can be overridden, but must set a list of \
      matplotlib Annotate objects in `self.feat_name`
                     
    Valid keyword arguments for :class:`BaseGraphicFeature` are:

        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        name                  a string text that will be drawn under the \
                              feature.
        type                  define the GrafhFeature type, is needed by  
                              `Panel` for feature grouping, default is ``''``
                              can be ``'plot'`` for plot-like features.
        score                 a float value that can be used to color the 
                              feature accordingly using a colormap
        height                feature y axis extension. default is ``1``
        cm                    matplotlib colormap to use, custom colormap can be
                              easily created but ar not supported right now. 
                              A list of default colormaps can be found at: 
                              www.scipy.org_
        color_by_cm           ``True`` | ``False``, if ``True`` set the color by
                              using the colormap. default is ``True``
        cm_value              float value between ``0.`` and ``1.`` to use for 
                              chosing the color from the colormap. default is 
                              ``None``
        use_score_for_color   ``True`` | ``False``. if ``True`` use the score 
                              value to pick a color from the colormap and will 
                              ovverride ``color_by_cm`` as ``True``.
                              if scores are not between 0 and 1 they can be 
                              normalized by using Normalize and LogNorm from 
                              ``matplotlib.colors``. 
        ec                    define edgecolor as in matplotlib,  overrides 
                              settings from  color_by_cm and use_score_for_color
        fc                    define facecolor as in matplotlib,  overrides 
                              settings from  color_by_cm and use_score_for_color
        lw                    defines edge line width, as in matplotlib
        alpha                 defines facecolor aplha value as in matplotlib
        boxstyle              default is ``'square, pad=0.'``, accepts any 
                              matplotlib ``FancyBboxPatch`` input string. 
                              take a look at matplotlib_ docs.
        norm                  normalization function to be used for in feature
                              colormap based coloring. tipically used in plot 
                              like features.default is ``None`` 
                              could be any function taking a  int or float value
                              and returning a float between  ``0.`` and ``1.``
        url                   set url fot the feature to be used in htmlmaps \
                              or SVG <a> elements.
        html_map_extend       included text will be added in html area tab, 
                              use to add onmouseover events, etc... 
                              must be a valid html attrib.
        ===================== ==================================================

.. _www.scipy.org: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

.. _matplotlib: http://matplotlib.sourceforge.net/api/artist_api.html?highlight\
=fancybboxpatch#matplotlib.patches.FancyBboxPatch

            '''
    
    default_cm = 'winter'    #deafult matplotlib colormap to use
    
    def __init__(self, **kwargs):
        '''
        init
        '''
        self.name = kwargs.get('name','')
        self.type = kwargs.get('type','')
        self.score = kwargs.get('score',0.)
        self.height = kwargs.get('height',1.)
        # html map options
        self.url = kwargs.get('url','')
        self.html_map_extend = kwargs.get('html_map_extend','')
        #feature style options
        self.cm = cm.get_cmap(kwargs.get('cm', self.default_cm))
        self.color_by_cm = kwargs.get('color_by_cm',True)
        self.cm_value = kwargs.get('cm_value',None)
        self.use_score_for_color = kwargs.get('use_score_for_color',False)
        if self.use_score_for_color:
            self.color_by_cm = True
        self.fc = kwargs.get('fc', self.cm(self.cm_value or 0))
        self.ec = kwargs.get('ec', None)
        self.lw = kwargs.get('lw', 0.5)
        self.ls = kwargs.get('ls','-')
        self.alpha = kwargs.get('alpha',.8)
        self.boxstyle = kwargs.get('boxstyle','square, pad=0.')
        self.Y = 0.0
        self.Y2 = 0.0
        self.patches = []   # all the patches must be returned inside this list
        self.feat_name = []   # all the patches labels must be returned inside this list
        self.norm = kwargs.get('norm', None)
        

    def find_arrow_dimensions(self, start, end, verbose=0):
        # verbose = 2
        # The following mess is needed because FancyArrowPatch only takes arrow head widths and lengths
        # in points (1/72 inches), the method is meant to draw arrows independently of axes sizes, same as a text annotation.
        # Because we want to set the arrow size relative to the axis, we have to transform between
        # coordinates systems: data, axis, display pixels and inches and points.
        # What a mess!!!

        # Adjust the arrow head length relative to the axis width
        axis_to_data = plt.gca().transAxes + plt.gca().transData.inverted()
        headLengthData = (axis_to_data.transform((arrowHeadLengthRelative, 0))[0] -
                          axis_to_data.transform((0, 0))[0])
        if verbose >= 2: print("headLengthData1", headLengthData)
        if headLengthData > end - start:
            headLengthData = 0.75*(end - start)
        if verbose >= 2: print("(end - start)", (end - start), "automaticHeadLengthData2", headLengthData)

        # Adjust body and head widths as fixed values in data coordinates
        bodyWidthData = self.height*1.0
        headWidthData = self.height*1.5

        # Convert all measures to x and y display inches and points
        data_to_display = plt.gca().transData
        xDataUnitInDisplayInches = (data_to_display.transform((1, 0))[0] -
                                    data_to_display.transform((0, 0))[0]) / plt.gcf().get_dpi()
        yDataUnitInDisplayInches = (data_to_display.transform((0, 1))[1] -
                                    data_to_display.transform((0, 0))[1]) / plt.gcf().get_dpi()
        xDataUnitInDisplayPt = xDataUnitInDisplayInches*72
        yDataUnitInDisplayPt = yDataUnitInDisplayInches*72
        if verbose >= 2:
            print("xDataUnitInDisplayInches", xDataUnitInDisplayInches)
            print("yDataUnitInDisplayInches", yDataUnitInDisplayInches)
            print("xDataUnitInDisplayPt", xDataUnitInDisplayPt)
            print("yDataUnitInDisplayPt", yDataUnitInDisplayPt)

        headLengthInches = headLengthData*xDataUnitInDisplayInches
        bodyWidthInches = bodyWidthData*yDataUnitInDisplayInches
        headWidthInches = headWidthData*yDataUnitInDisplayInches
        if verbose >= 2:
            print("headLengthInches", headLengthInches)
            print("headWidthData", headWidthData, "headWidthInches", headWidthInches)

        headLengthPoints = headLengthInches*72
        bodyWidthPoints = bodyWidthInches*72
        headWidthPoints = headWidthInches*72

        return (headLengthData, bodyWidthData, headWidthData,
                headLengthInches, bodyWidthInches, headWidthInches,
                headLengthPoints, bodyWidthPoints, headWidthPoints)



    def draw_feat_name(self,**kwargs):
        '''calling self.draw_feat_name() agailn will overwrite the previous \
        feature name stored in self.feat_name 
        Valid keyword arguments for :class:`BaseTrack` are:

        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        set_size              as in matplotlib, default is `'x-small'`
        set_family            as in matplotlib, default is `'serif'`
        set_weight            as in matplotlib, default is `'normal'`
        va                    as in matplotlib, default is `'top'`
        ha                    as in matplotlib, default is `'left'`  
        xoffset               to be used to move the text annotation along x axis
        ===================== ==================================================

        '''
        self.feat_name = []
        font_feat = FontProperties()
        set_size = kwargs.get('set_size','x-small')
        ax = kwargs.get('ax',plt.gca())
        set_family = kwargs.get('set_family','serif')
        set_weight = kwargs.get('set_weight','normal')
        va = kwargs.get('va','top')
        ha = kwargs.get('ha','left')
        xoffset = kwargs.get('xoffset',0)
        font_feat.set_size(set_size)
        font_feat.set_family(set_family)
        font_feat.set_weight(set_weight)
        text_x = self.start
        if (self.start < xoffset) and (xoffset <= self.end):
            text_x += xoffset
        #print("drawing feature name", self.name, (text_x, self.Y), (text_x, self.Y - self.height/5.), self.Y2)
        #self.feat_name = [ax.annotate(self.name, xy = (text_x, self.Y2), xytext = (text_x, self.Y2 - self.height/5.), fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va)]
        #for fname in self.feat_name:
        #    fname.set_visible(False)
        self.feat_name = [Text(x=text_x, y=self.Y2 - 0.12*self.height, text=self.name,
                          fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va)]

        
class Simple(BaseGraphicFeature):
    '''
    Handle a feature not derived from a SeqFeature. 
    Just need a `start` and `end` position.
    Minimum definition are start and end positions.
    Usage eg.
    
    ``features.Simple(start,end,name='factor 7',**kwargs)``

    '''
    def __init__(self,start,end,**kwargs):
    
        BaseGraphicFeature.__init__(self, **kwargs)
        self.start = start
        self.end = end


    def draw_feature(self):
        feat_draw = FancyBboxPatch((self.start,self.Y),
                                   width=(self.end-self.start),
                                   height=self.height,
                                   boxstyle=self.boxstyle,
                                   lw=self.lw,
                                   ec=self.ec,
                                   fc=self.fc,
                                   alpha=self.alpha,
                                   url=self.url,
                                   mutation_aspect=0.3)
        self.patches.append(feat_draw)


class GenericSeqFeature(BaseGraphicFeature):
    '''

    Handle a feature  derived from a `Biopython.SeqFeature.SeqFeature` as a \
    simple rectangle.
    
    Requires a `SeqFeature` object in input, `start` and `end` will be \
    automatically detected from seqfeature.
    
    Usage eg.

    ::
    
        from Bio.SeqFeature import SeqFeature, FeatureLocation 
        feat = SeqFeature (FeatureLocation(10,124))
        features.GenericSeqFeature(feat,name='factor 7', score=0.2, **kwargs)
    
    '''
    def __init__(self,feature,**kwargs):
        '''

        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        self.start = kwargs.get('start',min([feature.location.start.position,feature.location.end.position]))
        self.end = kwargs.get('end',max([feature.location.start.position,feature.location.end.position]))
        self.type = kwargs.get('type',feature.type)
        if 'score' in feature.qualifiers:
            self.score = kwargs.get('score',feature.qualifiers['score'])


    def draw_feature(self):
        feat_draw=FancyBboxPatch((self.start,self.Y), width=(self.end-self.start),
            height=self.height, boxstyle=self.boxstyle, lw=self.lw,
            ec=self.ec, fc=self.fc,alpha=self.alpha, mutation_aspect = 0.3, url = self.url,)
        self.patches.append(feat_draw)


class GeneSeqFeature(BaseGraphicFeature):
    '''
    
    Draws a Gene Feature as an arrow.
    Requires a SeqFeature with ``SeqFeature.type = 'gene'`` and optionally a 
    list of exons SeqFeatures that can be drawn over the gene patch.
    
    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        head_length           defines matplotlib  FancyArrow head_length keyword
                              argument to set up gene arrow head length
        ===================== ==================================================
    
    Usage eg.
    
    ::
    
        from Bio.SeqFeature import SeqFeature, FeatureLocation 
        genefeat = SeqFeature (FeatureLocation(103,1053), type = 'gene', strand=1,)
        features.GeneSeqFeature(genefeat,name='factor 7', **kwargs)
    
    
    '''
    def __init__(self,feature,exons=[],**kwargs):

        BaseGraphicFeature.__init__(self,**kwargs)
        self.height = kwargs.get('height',1.)
        self.start = kwargs.get('start',min([feature.location.start.position,feature.location.end.position]))
        self.end = kwargs.get('end',max([feature.location.start.position,feature.location.end.position]))
        default_head_length = 35
        if abs(self.end - self.start <= 50):
            default_head_length = (self.end-self.start)/10.
        self.head_length = kwargs.get('head_length', default_head_length)
        self.type = kwargs.get('type','gene')
        self.feature = feature
        self.exons = exons


    def draw_feature(self):

        (headLengthData, bodyWidthData, headWidthData,
         headLengthInches, bodyWidthInches, headWidthInches,
         headLengthPoints, bodyWidthPoints, headWidthPoints) = self.find_arrow_dimensions(self.start, self.end)

        self.patches = []
        if self.feature.strand == 1:
            arrow_start = self.start
            arrow_direction = self.end - self.start
        elif self.feature.strand == -1:
            arrow_start = self.end
            arrow_direction = self.start - self.end
        else:
            raise ValueError('Gene feature must have strand equal to 1 or -1')
        
        # feat_draw = FancyArrowPatch(
        #     (arrow_start, self.Y + 0.5),
        #     (arrow_start + arrow_direction, self.Y + 0.5),
        #     arrowstyle=ArrowStyle('simple',
        #                           head_length=headLengthInches,
        #                           head_width=headWidthInches,
        #                           tail_width=bodyWidthInches),
        #     ec=self.ec, fc=self.fc, alpha=self.alpha, lw=self.lw,
        #     mutation_scale=1)

        # feat_draw = FancyArrowPatch(
        #     (arrow_start, self.Y + 0.5),
        #     (arrow_start + arrow_direction, self.Y + 0.5),
        #     arrowstyle=ArrowStyle('simple',
        #                           head_length=headLengthPoints,
        #                           head_width=headWidthPoints,
        #                           tail_width=bodyWidthPoints),
        #     ec=self.ec, fc=self.fc, alpha=self.alpha, lw=self.lw,
        #     mutation_scale=1)

        feat_draw = FancyArrow(
            x=arrow_start, y=self.Y + self.height/2., dx=arrow_direction, dy=0,
            ec=self.ec, fc=self.fc, alpha=self.alpha, lw=self.lw,
            width=bodyWidthData,
            head_length=headLengthData,
            head_width=headWidthData,
            length_includes_head=True, head_starts_at_zero=False)

        self.patches.append(feat_draw)

        # for exon in self.exons:
        #     feat_draw = FancyBboxPatch(
        #         (int(exon.location.start.position),self.Y),
        #         width=(int(exon.location.end.position)-int(exon.location.start.position)),
        #         height=body_width/2., boxstyle=self.boxstyle,lw=0, ec=self.ec,
        #         fc=self.fc,alpha=self.alpha+0.1, url=self.url,)
        #     self.patches.append(feat_draw)


class TextSequence(BaseGraphicFeature):
    '''
    
    Draws a biological sequence as text at the given positions   
    
    Requires a `sequence` iterable tipically a string.
    
    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        start                 feature start position. default is ``0``
        end                   feature end position. default is 
                              ``self.start +  len(self.sequence)``
        set_size              as in matplotlib, default is ``'xx-small'``
        set_family            as in matplotlib, default is ``'monospace'``     
        set_weight            as in matplotlib, default is ``'normal'``     
        va                    as in matplotlib, default is ``'bottom'`` 
        ha                    as in matplotlib, default is ``'center'``
        ===================== ==================================================
    
    Usage eg.
    
    ::
    
        features.TextSequence('ATCATGATACTGNCAT',start = 103, **kwargs)
    
    
    '''
    def __init__(self,sequence,**kwargs):

        BaseGraphicFeature.__init__(self,**kwargs)
        self.height = kwargs.get('height', textFeatureHeight)
        self.sequence = sequence
        self.start = kwargs.get('start',0)
        self.end = kwargs.get('end', self.start + len(self.sequence))
        self.customFontSize = 'set_size' in kwargs
        self.set_size = kwargs.get('set_size','medium')
        self.set_family = kwargs.get('set_family','monospace')
        self.set_weight = kwargs.get('set_weight','normal')
        self.va = kwargs.get('va','bottom')
        self.ha = kwargs.get('ha','center')
        self.fontProp = FontProperties()
        self.fontProp.set_size(self.set_size)
        self.fontProp.set_family(self.set_family)
        self.fontProp.set_weight(self.set_weight)
        
    def draw_feature(self):
        if not self.customFontSize:
            # Adjust default text size such that we uses the biggest font size possible to fill all the space
            # between two positions. The idea is to make the font size proportional to the distance
            # in inches between position i and i+1 (distance 1 in axis coordinates).
            width1AxisUnitInDisplayPixels = plt.gca().transData.transform((1, 0))[0] - plt.gca().transData.transform((0, 0))[0]
            width1AxisUnitInDisplayInches = width1AxisUnitInDisplayPixels / plt.gcf().get_dpi()
            self.fontProp.set_size(width1AxisUnitInDisplayInches*defaultFontSizeScalingFactor)

        self.patches = []
        for i, char in enumerate(self.sequence):
            self.patches.append(Text(x=self.start + i,
                                     y=0,
                                     text=char,
                                     fontproperties=self.fontProp,
                                     horizontalalignment=self.ha,
                                     verticalalignment=self.va,
                                     url=self.url))



class PrettyTextSequence(BaseGraphicFeature):
    """
    characterList should be a list of (character, attributeDictionary).
    """

    def __init__(self, sequence, highlightList=None, **kwargs):

        BaseGraphicFeature.__init__(self,**kwargs)
        self.sequence = sequence
        self.highlightList = highlightList
        self.height = kwargs.get('height', textFeatureHeight)
        self.start = kwargs.get('start', 0)
        self.end = kwargs.get('end', self.start + len(self.sequence))
        self.hasCustomGeneralFontSize = 'set_size' in kwargs
        self.set_size = kwargs.get('set_size','medium')
        self.set_family = kwargs.get('set_family','monospace')
        self.set_weight = kwargs.get('set_weight','normal')
        self.va = kwargs.get('va','bottom')
        self.ha = kwargs.get('ha','center')
        self.fontProp = FontProperties()
        self.fontProp.set_size(self.set_size)
        self.fontProp.set_family(self.set_family)
        self.fontProp.set_weight(self.set_weight)
        
    def draw_feature(self):

        if not self.hasCustomGeneralFontSize:
            # Adjust default text size such that we uses the biggest font size possible to fill all the space
            # between two positions. The idea is to make the font size proportional to the distance
            # in inches between position i and i+1 (distance 1 in axis coordinates).
            width1AxisUnitInDisplayPixels = plt.gca().transData.transform((1, 0))[0] - plt.gca().transData.transform((0, 0))[0]
            width1AxisUnitInDisplayInches = width1AxisUnitInDisplayPixels / plt.gcf().get_dpi()
            automaticFontSize = width1AxisUnitInDisplayInches*defaultFontSizeScalingFactor
            self.fontProp.set_size(automaticFontSize)

        self.patches = []

        # Transform the sequence in a list of (char, propArgsDict)
        if type(self.sequence) is str:
            self.sequence = list(self.sequence)

        for i, charTuple in enumerate(self.sequence):
            if type(charTuple) is str:
                # We only have a character without any custom font properties
                char = charTuple
                propArgsDict = {'fontproperties':self.fontProp.copy()}

            # Add default font properties if needed
            if 'fontproperties' not in propArgsDict:
                propArgsDict['fontproperties'] = self.fontProp.copy()

            # Update char tuple
            self.sequence[i] = (char, propArgsDict)

        # Apply highlighting
        if self.highlightList is not None:
            for highlightRegion in self.highlightList:
                bc = highlightRegion.get('background_color')
                if bc is None:
                    bc = highlightRegion.get('bc')
                bc_alpha = highlightRegion.get('background_color_alpha')
                if bc_alpha is None:
                    bc_alpha = highlightRegion.get('bc_alpha')
                fc = highlightRegion.get('foreground_color')
                if fc is None:
                    fc = highlightRegion.get('fc')
                weight = highlightRegion.get('weight')

                for i in range(highlightRegion['start'], highlightRegion['end']):
                    char, propArgsDict = self.sequence[i]
                    if bc is not None:
                        bc_alpha = 1.0 if bc_alpha is None else bc_alpha
                        propArgsDict['bbox'] = dict(facecolor=bc, alpha=bc_alpha, edgecolor='none', lw=0, pad=0.1)
                    if fc is not None:
                        propArgsDict['color'] = fc
                    if weight is not None:
                        propArgsDict['fontproperties'].set_weight(weight)
                    # print("highlight", i, fc, bc, bc_alpha, weight, "propArgsDict", propArgsDict)
            
        for i, charTuple in enumerate(self.sequence):
            (char, propArgsDict) = charTuple

            # Adjust the font size automatically, unless explicitly defined as custom size in the font properties dictionary
            # Note that here we override any font size that would be defined in the FontProperties object passed in 
            # `propArgsDict['fontproperties']`. This is because there is no way to know if the FontProperties object uses
            #  the default font size or a user defined font size. We make it explicit by requiring the custom font size
            # as a separate option in the `propArgsDict`.
            customFontSize = propArgsDict.get('size')
            if customFontSize is not None:
                propArgsDict['fontproperties'].set_size(customFontSize)
            else:
                propArgsDict['fontproperties'].set_size(automaticFontSize)

            # If the font size is too small, we do not draw any feature
            if propArgsDict['fontproperties'].get_size() > 3:
                feat_draw = Text(x=self.start + i,
                                 y=0,
                                 text=char,
                                 horizontalalignment=self.ha,
                                 verticalalignment=self.va,
                                 **propArgsDict)
                    # # Here we plot directly the text on the current figure, however we will have to remove from it afterwards
                    # feat_draw = fig.text(self.start + i, 0, char,
                    #                            horizontalalignment=self.ha,
                    #                            verticalalignment=self.va,
                    #                            **propArgsDict)
                self.patches.append(feat_draw)



class PlotFeature(BaseGraphicFeature):
    ''' 
    
    Draws a lineplot of continuous values for a feature such as an 
    hydrophobicity plot, or prediction confidency values
    must be used within a :class:`~biograpy.tracks.PlotTrack`.

    Y values are necessary
    If no X values are provided they will be created starting from position 1 
    for each Y value provided
    
    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        y                     list of Y numerical values  
        x                     list of X numerical values  
        style                 line plotting style as in matplotlib. default is
                              ``'-'``
        label                 set label to be used in `PlotTrack` legend. if no
                              label is defined for all the features the legend
                              will not be drawn. 
        ===================== ==================================================

    Usage eg.
    
    ::
    
        features.PlotFeature([1.3, 2.4, -1,2],  **kwargs)
    
    '''
    def __init__(self,  y , x = [], style = '-' , **kwargs):

        BaseGraphicFeature.__init__(self,**kwargs)
        self.feat_type = 'plot' #to be checked in plot track
        self.label = kwargs.get('label', '')

        if  not isinstance(y[0], (int,float)):
            raise ValueError('PlotFeature objects only accepts int or float data')
        self.x = x
        self.y = y
        self.style = style
        if not self.x:
            self.x = range(1,len(self.y)+1)
        self.start = min(self.x)
        self.end = max(self.x)

            
    def draw_feature(self):
        if self.y:
            plotted_data = plt.plot(self.x,
                                    self.y,
                                    self.style,
                                    label = self.label,
                                    lw = self.lw,
                                    ls =self.ls,
                                    #fc = self.fc,
                                    #ec = self.ec,
                                    alpha = self.alpha,
                                    url = self.url,)
            self.patches.extend(plotted_data)
            
        

class BarPlotFeature(BaseGraphicFeature):
    '''
    
    
    Draws a bar plot values for a feature such as prediction confidence values
    must be used within a PlotTrak to have axix.
   
    
   
    Draws a barplot of continuous values for a feature such as 
    prediction confidency values must be used within a 
    :class:`~biograpy.tracks.PlotTrack`.

    Y values are necessary
    If no X values are provided they will be created starting from position 1 
    for each Y value provided
    
    ``color_by_cm = True`` will color each point based on the y values using
    the supplied colormap ``cm``.
    
    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        y                     list of Y numerical values  
        x                     list of X numerical values  
        style                 line plotting style as in matplotlib. default is
                              ``'-'``
        label                 set label to be used in `PlotTrack` legend. if no
                              label is defined for all the features the legend
                              will not be drawn. 
        width                 as in matplotlib, default is ``.8``
        bottom                as in matplotlib, default is ``0``
        align                 as in matplotlib, default is ``'center'``
        xerr                  as in matplotlib, default is ``None``
        yerr                  as in matplotlib, default is ``None``
        ecolor                as in matplotlib, default is ``'k'``
        capsize               as in matplotlib, default is ``1``
        ls                    as in matplotlib, default is ``'solid'``
        ===================== ==================================================

    Usage eg.
    
    ::
    
        feat = features.PlotFeature([1.3, 2.4, -1,2],  **kwargs)


    '''
    def __init__(self,  y , x = [], style = '-' , **kwargs):

        BaseGraphicFeature.__init__(self,**kwargs)
        self.feat_type = 'plot' #to be checked in plot track
        self.label = kwargs.get('label', '')
        self.width = kwargs.get('width', .8)
        self.bottom = kwargs.get('bottom', 0)
        self.align = kwargs.get('align', 'center')
        self.xerr = kwargs.get('xerr', None)
        self.yerr = kwargs.get('yerr', None)
        self.ecolor = kwargs.get('ecolor', 'k')
        self.capsize = kwargs.get('capsize', 1)
        self.ls = kwargs.get('ls', 'solid')#['solid' | 'dashed' | 'dashdot' | 'dotted']
        if  not isinstance(y[0], (int,float)):
            raise ValueError('PlotFeature objects only accepts int or float data')
        self.x = x
        self.y = y
        self.style = style
        if not self.x:
            self.x = range(1,len(self.y)+1)
        self.start = min(self.x)
        self.end = max(self.x)

            
    def draw_feature(self):
        if self.y:
            bars = plt.bar(self.x,
                           self.y,
                           width = self.width,
                           bottom = self.bottom, 
                           label = self.label,
                           lw = self.lw,
                           ls =self.ls,
                           fc = self.fc,
                           ec = self.ec,
                           alpha = self.alpha,
                           xerr = self.xerr,
                           yerr = self.yerr,
                           ecolor = self.ecolor,
                           capsize = self.capsize,
                           align = self.align,
                           url = self.url,
                           )
                        
            if self.color_by_cm:
                self.norm = Normalize(min(self.y), max(self.y))
                for bar, value in zip(bars, self.y):
                    color = self.cm(self.norm(value))
                    bar.set_color(color)
                    
            self.patches.extend(bars)
            
        
        


class SegmentedSeqFeature(BaseGraphicFeature):
    '''
    
    Draws a SeqFeature carrying ``'join'`` subfeatures
    
    requires a SeqFeature object carryng ``'joins'`` in SeqFeature.subfeatures, 
    `start` and `end` can be supplied or automatically detected from SeqFeature.
    
    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        start                 feature start position. default is the minimum 
                              value of all SubFeature.location postions
        end                   feature end position. default is the maximum 
                              value of all SubFeature.location postions  
        type                  feature type
        ec                    default is set to ``'k'``
        ===================== ==================================================
    
    Usage eg.
    
    ::
    
        feat = SegmentedSeqFeature(feature,name='factor 7',s)
    '''

    def __init__(self,feature,**kwargs):
        '''
        
        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        self.start=kwargs.get('start',min([feature.location.start.position,feature.location.end.position]))
        self.end=kwargs.get('end',max([feature.location.start.position,feature.location.end.position]))
        self.type=kwargs.get('type',feature.type)
        self.ec=kwargs.get('ec','k')
        if 'score' in feature.qualifiers:
            self.score=kwargs.get('score',feature.qualifiers['score'])
        self.feature=feature


    def draw_feature(self):
        if len(self.feature.sub_features):
            junction_start=False
            if self.feature.strand==-1:
                self.feature.sub_features.reverse()
            for sub_feature in self.feature.sub_features:
                if sub_feature.location_operator=='join':
                    feat_draw=FancyBboxPatch((int(sub_feature.location.start.position),self.Y), width=(int(sub_feature.location.end.position)-int(sub_feature.location.start.position)), height=self.height, boxstyle=self.boxstyle,lw=self.lw, ec=self.ec, fc=self.fc, alpha=self.alpha, url = self.url,)
                    self.patches.append(feat_draw)
                    if junction_start:
                        junction_end=float(sub_feature.location.start.position)
                        junction_middle=float((junction_start+junction_end)/2.)
                        Yends=self.Y+self.height/2.
                        Ymiddle=self.Y+self.height
                        join=plt.plot([junction_start,junction_middle,junction_end],[Yends,Ymiddle,Yends], lw=1 ,ls='-', c=self.ec, alpha=0.5, url = self.url,)
                        self.patches.extend(join)
                    junction_start=float(sub_feature.location.end.position)
        else:#if SegmentedSeqFeature is called whihout subfeatures returns a GenericSeqFeature
            feat_draw=FancyBboxPatch((self.start,self.Y), width=(self.end-self.start), height=self.height, boxstyle=self.boxstyle, lw=self.lw, ec=self.ec, fc=self.fc,alpha=self.alpha,)
            self.patches.append(feat_draw)


class CoupledmRNAandCDS(BaseGraphicFeature):
    '''

    Draws an mRNa at a small alpha value, and the corresponding CDS on top of 
    it at `alpha = 1.`
    Requires an mRNA SeqFeature object and a CDS seqfeature object. 'joins' will 
    be threated as in :class:`~biograpy.features.SegmentedSeqFeature`

    
    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        mRNA_start            mRNA feature start position. default is the 
                              minimum value of all SubFeature.location postions
        mRNA_end              mRNA feature end position. default is the maximum 
                              value of all SubFeature.location postions  
        CDS_start             CDS feature start position. default is the 
                              minimum value of all SubFeature.location postions
        CDS_end               CDS feature end position. default is the maximum 
                              value of all SubFeature.location postions  
        type                  feature type
        ec                    default is set to ``'k'``
        ===================== ==================================================
    
    Usage eg.
    
    ::
    
        feat = CoupledmRNAandCDS(mRNA_feature, CDS_feature, name='factor 7',s)

    
    '''
    def __init__(self, mRNA, CDS, **kwargs):
        '''
        
        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        self.mRNA_start = kwargs.get('mRNA_start',min([mRNA.location.start.position,mRNA.location.end.position]))
        self.mRNA_end = kwargs.get('mRNA_end',max([mRNA.location.start.position,mRNA.location.end.position]))
        self.CDS_start = kwargs.get('CDS_start',min([CDS.location.start.position,CDS.location.end.position]))
        self.CDS_end = kwargs.get('CDS_end',max([CDS.location.start.position,CDS.location.end.position]))
        self.start = self.mRNA_start
        self.end = self.mRNA_end
        self.type = kwargs.get('type','mRNA')
        if 'score' in mRNA.qualifiers:
            self.score = kwargs.get('score',mRNA.qualifiers['score'])
        self.mRNA = mRNA
        self.CDS = CDS
        self.ec = kwargs.get('ec','k')
        # Arrow head length
        default_head_length = 35
        if abs(self.end - self.start <= 50):
            default_head_length = (self.end-self.start)/10.
        self.head_length = kwargs.get('head_length', default_head_length)


    def draw_feature(self):

        #draw mRNA
        if type(self.mRNA.location) is CompoundLocation:
            junction_start = False
            if self.mRNA.strand == -1:
                self.mRNA.sub_features.reverse()
            for sub_feature in self.mRNA.sub_features:
                if sub_feature.location_operator == 'join':
                    feat_draw=FancyBboxPatch((int(sub_feature.location.start.position),self.Y), width=(int(sub_feature.location.end.position)-int(sub_feature.location.start.position)), height=self.height, boxstyle=self.boxstyle,lw=self.lw, ec=self.ec, fc=self.fc,alpha=self.alpha-0.7, url = self.url,)
                    self.patches.append(feat_draw)
                    if junction_start:
                        junction_end=float(sub_feature.location.start.position)
                        junction_middle=float((junction_start+junction_end)/2.)
                        Yends=self.Y+self.height/2.
                        Ymiddle=self.Y+self.height
                        join=plt.plot([junction_start,junction_middle,junction_end],[Yends,Ymiddle,Yends], lw=0.5 ,ls='-', c=self.ec,alpha=0.5)
                        self.patches.extend(join)
                    junction_start=float(sub_feature.location.end.position)
        else:
            feat_draw = FancyBboxPatch((self.mRNA_start, self.Y), width=(self.mRNA_end-self.mRNA_start),
                                       height=self.height, boxstyle=self.boxstyle, lw=self.lw, ec=self.ec,
                                       fc=self.fc,alpha=self.alpha-0.5, url=self.url,)
            self.patches.append(feat_draw)

        # draw CDS
        if type(self.CDS.location) is CompoundLocation:
            # Note: this is probably not correct for the newer Biopython CompoundLocation class...
            junction_start = False
            if self.CDS.strand == -1:
                self.CDS.sub_features.reverse()
            for sub_feature in self.CDS.sub_features:
                if sub_feature.location_operator == 'join':
                    feat_draw = FancyBboxPatch((int(sub_feature.location.start.position),self.Y), width=(int(sub_feature.location.end.position)-int(sub_feature.location.start.position)), height=self.height, boxstyle=self.boxstyle,lw=0, ec=self.ec, fc=self.fc,alpha=self.alpha,)
                    self.patches.append(feat_draw)
        else:
            (headLengthData, bodyWidthData, headWidthData,
             headLengthInches, bodyWidthInches, headWidthInches,
             headLengthPoints, bodyWidthPoints, headWidthPoints) = self.find_arrow_dimensions(self.CDS_start, self.CDS_end)

            if self.CDS.strand == 1:
                arrow_start = self.CDS_start
                arrow_direction = self.CDS_end-self.CDS_start
            elif self.CDS.strand == -1:
                arrow_start = self.CDS_end
                arrow_direction = self.CDS_start-self.CDS_end
            else:
                raise ValueError('Gene feature must have strand equal to 1 or -1')
            shape = 'full'   # if we prefer a full arrow instead of half arrow head.

            feat_draw = FancyArrow(
                x=arrow_start, y=self.Y + self.height/2., dx=arrow_direction, dy=0,
                ec=self.ec, fc=self.fc, alpha=self.alpha, lw=self.lw,
                width=bodyWidthData,
                head_length=headLengthData,
                head_width=headWidthData,
                length_includes_head=True, head_starts_at_zero=False)

            self.patches.append(feat_draw)

            

            
            
class SinglePositionFeature(BaseGraphicFeature):
    '''
 
    
    Draws features localized just to one position

    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        marker                marker symbol, accepts all matplotlib marker types
                              default is circle marker ``'o'``
        markersize            marker size as in matplotlib. default  is ``4``.
        ===================== ==================================================
        
    Available matplotlib marker types:
        
            
        * ``'.'``     point marker
        * ``','``     pixel marker
        * ``'o'``     circle marker
        * ``'v'``     triangle_down marker
        * ``'^'``     triangle_up marker
        * ``'<'``     triangle_left marker
        * ``'>'``     triangle_right marker
        * ``'1'``     tri_down marker
        * ``'2'``     tri_up marker
        * ``'3'``     tri_left marker
        * ``'4'``     tri_right marker
        * ``'s'``     square marker
        * ``'p'``     pentagon marker
        * ``'*'``     star marker
        * ``'h'``     hexagon1 marker
        * ``'H'``     hexagon2 marker
        * ``'+'``     plus marker
        * ``'x'``     x marker
        * ``'D'``     diamond marker
        * ``'d'``     thin_diamond marker
        * ``'|'``     vline marker
        * ``'_'``     hline marker
    
    Usage eg.
    
    ::
    
        feat = SinglePositionFeature(feature, name='N-glycosilation',s)


    
    '''
    def __init__(self,feature,**kwargs):
        '''
        
        '''
        BaseGraphicFeature.__init__(self,**kwargs)
       
        self.start = self.end = kwargs.get('start',min([feature.location.start.position,feature.location.end.position]))
        self.type=kwargs.get('type',feature.type)
        if 'score' in feature.qualifiers:
            self.score=kwargs.get('score',feature.qualifiers['score'])

        self.marker=kwargs.get('marker','o')
        self.markersize=kwargs.get('markersize',4)
        

        
    def draw_feature(self):
        feat_draw=plt.plot(self.start, self.Y, marker=self.marker, markerfacecolor=self.fc, markeredgecolor='k', markersize=self.markersize, alpha=self.alpha, url = self.url,)
        self.patches=feat_draw

        
        
class TMFeature(BaseGraphicFeature):
    '''  
    
    Draws  a line plot of proteins contaning transmembrane regions with 
    cytoplasmic/non cytoplasmic sides.
    
    A list of `transmembrane` SeqFeatures is required. additional `cytoplasmic`
    and `non cytoplasmic` SeqFeatures can be passed. If only `transmembrane` 
    SeqFeatures are given, the remaining space is filled with an horizontal 
    line.

    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        TM                    list of SeqFeatures reporting Transmembrane 
                              regions. required to draw something.
        cyto                  list of SeqFeatures reporting cytoplasmic regions.
        non_cyto              list of SeqFeatures reporting non cytoplasmic 
                              regions. 
        type                  feature type
        TM_ec                 transmembrane region edgecolor, as in matplotlib.
        TM_fc                 transmembrane region facecolor, as in matplotlib.
        TM_label              text that will appear under the transmembrane 
                              regions. default is ``'TM'``
        cyto_color            cytoplasmic line color default is ``self.ec``
        cyto_linestyle        cytoplasmic line style default is solid ``'-'``
        cyto_label            text that will appear under the cytoplasmic 
                              regions. default is ``'cyto'``
        non_cyto_color        non cytoplasmic line color default is ``self.ec``
        non_cyto_linestyle    non cytoplasmic line style default is solid 
                              ``'-'``
        non_cyto_label        text that will appear under the non cytoplasmic 
                              regions. default is ``'non cyto'``
        fill_color            automatic fill line color default is ``self.ec``         
        fill_linestyle        automatic fill line style default is solid ``'-'``
        connection_lw         linewidth for cytoplasmic, non cytoplasming and
                              automatic filled lines. default is ``'2``
        do_not_fill           ``True`` | ``False``. if ``True`` will not 
                              automatically fill the space between transmembrane
                              regions with a line. default is ``False``
        ===================== ==================================================

        
        Available linestyles for lines:
        
            * ``'-'``     solid line style
            * ``'--'``    dashed line style
            * ``'-.'``    dash-dot line style
            * ``':'``     dotted line style
    
    Usage eg.
            
    ::
    
        feat = TMFeature(TM_features, cyto_features, non_features, 
                         non_cyto_label= 'extracellular', )
    
    '''
    def __init__(self, **kwargs):
        '''
        
        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        
        self.TM = kwargs.get('TM',[])
        self.cyto = kwargs.get('cyto',[])
        self.non_cyto = kwargs.get('non_cyto',[])
        self.type=kwargs.get('type',self.TM[0].type)
        self.TM_ec = kwargs.get('TM_ec',self.ec)
        self.TM_fc = kwargs.get('TM_fc',self.fc)
        self.TM_label = kwargs.get('TM_label','TM')
        self.cyto_color = kwargs.get('cyto_color',self.ec)
        self.cyto_linestyle = kwargs.get('cyto_linestyle','-')
        self.cyto_label = kwargs.get('cyto_label','cyto')
        self.non_cyto_color = kwargs.get('non_cyto_color',self.ec)
        self.non_cyto_linestyle = kwargs.get('non_cyto_linestyle','-')
        self.non_cyto_label = kwargs.get('non_cyto_label','non cyto')
        self.fill_color = kwargs.get('fill_color',self.cyto_color)
        self.fill_linestyle = kwargs.get('fill_linestyle',self.cyto_linestyle)
        self.connection_lw = kwargs.get('connection_lw',2)
        self.do_not_fill = kwargs.get('do_not_fill',False)
        
        Xs = []
        for feature in self.TM:
            Xs.append(feature.location.start.position)
            Xs.append(feature.location.end.position)
        if self.cyto:
            for feature in self.cyto:
                Xs.append(feature.location.start.position)
                Xs.append(feature.location.end.position)
        if self.non_cyto:
            for feature in self.non_cyto:
                Xs.append(feature.location.start.position)
                Xs.append(feature.location.end.position)
        self.start = min(Xs)
        self.end = max(Xs)
        
        if (self.cyto and self.non_cyto) or self.do_not_fill:
            self.fill = False 
        else:
            self.fill = True 
            
        #label positions
        self.TM_starts = []
        self.cyto_starts = []
        self.non_cyto_starts = []
        
        

    
    
    def draw_feature(self):
        
        def draw_orizontal_line(start, end, y, linestyle, color, url):
            return plt.plot((start,end), (y,y), linestyle=linestyle, color=color, 
                            lw=self.connection_lw, aa=False, alpha=self.alpha,
                            url = url)
        
        TMs=[]
        for feature in self.TM:
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            type = feature.type
            self.TM_starts.append((start,end))
            if 'score' in feature.qualifiers:#not used yet
                score=kwargs.get('score',feature.qualifiers['score'])
            feat_draw=FancyBboxPatch((start,self.Y), width=(end-start), height=self.height, boxstyle=self.boxstyle, lw=self.lw, ec=self.TM_ec, fc=self.TM_fc,alpha=self.alpha, url = self.url,)
            self.patches.append(feat_draw)
        if self.fill:
            start = 1
            for TM_start,TM_end in sorted(self.TM_starts):
                end = TM_start# - 1
                feat_draw=draw_orizontal_line(start, end, self.Y + self.height/2, self.fill_linestyle, self.fill_color, url = self.url,)
                self.patches.extend(feat_draw)#use extend, since plot returns a list of Line2D objects
                start = TM_end #+ 1
        'draw cytplasmic regions'
        for feature in self.cyto:
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            self.cyto_starts.append((start,end))
            feat_draw=draw_orizontal_line(start, end, self.Y + self.height/2, self.cyto_linestyle, self.cyto_color, url = self.url,)
            self.patches.extend(feat_draw)
        'draw non cytplasmic regions'
        for feature in self.non_cyto:
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            self.non_cyto_starts.append((start,end))
            feat_draw=draw_orizontal_line(start, end, self.Y + self.height/2, self.non_cyto_linestyle, self.non_cyto_color, url = self.url,)
            self.patches.extend(feat_draw)
                 
                
    def draw_feat_name(self,**kwargs):
        '''draws the feature name and all the TM, cyto a non cyto labels'''
        self.feat_name=[]
        font_feat=FontProperties()
        set_size=kwargs.get('set_size','xx-small')
        set_family=kwargs.get('set_family','serif')
        set_weight=kwargs.get('set_weight','normal')
        va=kwargs.get('va','top')
        ha=kwargs.get('ha','center')
        font_feat.set_size(set_size)
        font_feat.set_family(set_family)       
        font_feat.set_weight(set_weight)  
        
        for start, end in self.TM_starts:
            #self.feat_name.append(plt.text(start+((end-start)/2), self.Y-self.height/5., self.TM_label, fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va))
            self.feat_name.append(plt.annotate(self.TM_label, xy = (start + (end - start) / 2, self.Y), xytext = (start + (end - start) / 2, self.Y - self.height/5.), fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va))
        for start, end in self.cyto_starts:
            #self.feat_name.append(plt.text(start+((end-start)/2), self.Y-self.height/5., self.cyto_label, fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va))
            self.feat_name.append(plt.annotate(self.cyto_label, xy = (start + (end - start) / 2, self.Y), xytext = (start + (end - start) / 2, self.Y - self.height/5.), fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va))
        for start, end in self.non_cyto_starts:
            #self.feat_name.append(plt.text(start+((end-start)/2), self.Y-self.height/5., self.non_cyto_label, fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va))
            self.feat_name.append(plt.annotate(self.non_cyto_label, xy = (start + (end - start) / 2, self.Y), xytext = (start + (end - start) / 2, self.Y - self.height/5.), fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va))

            
            
class SecStructFeature(BaseGraphicFeature):
    '''
    
    Draws secondary structures. requires at least a list of SeqFeatures 
    indicating beta strands (`betas`), alpha helices (`alphah`), and random
    coil (`coil`).

    Beta strands are drawn as arrows, alpha helices as rectangles and coils as
    horizontal lines. 


    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        betas_ec              beta strands edgecolor, as in matplotlib.
        betas_fc              beta strands facecolor, as in matplotlib.
        alphah_ec             alpha helices edgecolor, as in matplotlib.
        alphah_fc             alpha helices facecolor, as in matplotlib.
        coil_color            random coil line color default is ``self.ec``
        coil_linestyle        random coil line style default is solid 
        lw                    linewidth for coil lines . default is ``'0.5'``
        type                  feature type. default is ``''``
        filter_struct_length  if != 0 smooths secondary structures, by 
                              displaying only those longher than the given 
                              value. Force coil generation, and ignore supplied 
                              coil regions. default is ``'0'``
        ===================== ==================================================
        
        Available linestyles for random coils:
        
            * ``'-'``     solid line style
            * ``'--'``    dashed line style
            * ``'-.'``    dash-dot line style
            * ``':'``     dotted line style

    Usage eg.
            
    ::
    
        feat = SecStructFeature(betas = betas_features, 
                                alphah = alphah_features)
    
    
    '''
    def __init__(self, betas = [], alphah = [], coil = [], **kwargs):
        '''

        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        
        self.betas = betas
        self.alphah = alphah
        self.coil = coil
        self.type=kwargs.get('type','')
        
        self.betas_ec = kwargs.get('betas_ec','k')
        self.betas_fc = kwargs.get('betas_fc','yellow')
        self.alphah_ec = kwargs.get('alphah_ec','k')
        self.alphah_fc = kwargs.get('alphah_fc','magenta')
        self.coil_color = kwargs.get('coil_color','k')
        self.coil_linestyle = kwargs.get('coil_linestyle','-')
        self.lw = kwargs.get('lw',2)
        
        self.filter_struct_length = kwargs.get('filter_struct_length', 0)#if != 0 smooths secondary structures, by displaying only those longher than the given value. self.filter_struct_length force coil geneation, and ignore supplied coils regions
        
        self.structured_regions =[]

        Xs = []
        if self.betas:
            for feature in self.betas:
                Xs.append(feature.location.start.position)
                Xs.append(feature.location.end.position)
        if self.alphah:
            for feature in self.alphah:
                Xs.append(feature.location.start.position)
                Xs.append(feature.location.end.position)
        if self.coil:
            for feature in self.coil:
                Xs.append(feature.location.start.position)
                Xs.append(feature.location.end.position)
        self.start = min(Xs)
        self.end = max(Xs)


    def draw_feature(self):
        def draw_orizontal_line(start, end, y, linestyle, color, url):
            return plt.plot((start,end), (y,y), linestyle=linestyle, color=color, 
                                lw=1, aa=False, alpha=self.alpha, url = url)
        'draw beta strands'
        for feature in self.betas:
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            if (end-start) >= self.filter_struct_length: #skip short secondary structures
                type = feature.type
                self.structured_regions.append((start,end))
                if 'score' in feature.qualifiers:#not used yet
                    score=kwargs.get('score',feature.qualifiers['score'])

                feat_draw = FancyArrow(start, self.Y+self.height/2., dx=end-start, dy=0, ec=self.betas_ec, fc=self.betas_fc, alpha=self.alpha, 
                                       width=self.height/2., head_length=(end-start)*.33, head_width=self.height, 
                                       lw=self.lw, length_includes_head=True,  head_starts_at_zero=False, url = self.url,)
                self.patches.append(feat_draw)


        'draw alpha helix'
        for feature in self.alphah:
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            if (end-start) >= self.filter_struct_length: #skip short secondary structures
                type = feature.type
                self.structured_regions.append((start,end))
                if 'score' in feature.qualifiers:#not used yet
                    score=kwargs.get('score',feature.qualifiers['score'])

                feat_draw = FancyBboxPatch((start,self.Y), width=(end-start), height=self.height, boxstyle=self.boxstyle, 
                                           lw=self.lw, ec=self.alphah_ec, fc=self.alphah_fc,alpha=self.alpha, url = self.url,)
                self.patches.append(feat_draw)
        
        if (not self.coil) or (self.filter_struct_length):
            start = 1
            for region_start,region_end in sorted(self.structured_regions):
                end = region_start #- 1
                feat_draw=draw_orizontal_line(start, end, self.Y + self.height/2, self.coil_linestyle, self.coil_color, url = self.url,)
                self.patches.extend(feat_draw)#use extend, since plot returns a list of Line2D objects
                start = region_end #+ 1
        else:
             for feature in self.coil:
                start = min([feature.location.start.position,feature.location.end.position])
                end = max([feature.location.start.position,feature.location.end.position])
                type = feature.type
                if 'score' in feature.qualifiers:#not used yet
                    score=kwargs.get('score',feature.qualifiers['score'])
                feat_draw=draw_orizontal_line(start, end, self.Y + self.height/2, self.coil_linestyle, self.coil_color, url = self.url,)
                self.patches.extend(feat_draw)#use extend, since plot returns a list of Line2D objects
                
                
class DomainFeature(BaseGraphicFeature):
    '''
    
    Draws one ore more domains (encoded as seqfeatures) with `FancyBboxPatch`.
    Can optionally draw a line representing the complete sequence.
    Accepts a list of SeqFeature objects.
    
    Additional valid attributes:
        ===================== ==================================================
        Key                   Description
        ===================== ==================================================
        seq_line              sequence line start and end positions tuple. 
                              default is (0,0)
        name                  inherited by 
                              :class:`~biograpy.features.BaseGraphicFeature`,
                              accepts multiple values
        boxstyle              inherited by 
                              :class:`~biograpy.features.BaseGraphicFeature`,
                              accepts multiple values
        ec                    inherited by 
                              :class:`~biograpy.features.BaseGraphicFeature`, 
                              accepts multiple values
        fc                    inherited by 
                              :class:`~biograpy.features.BaseGraphicFeature`, 
                              accepts multiple values
        alpha                 inherited by
                              :class:`~biograpy.features.BaseGraphicFeature`,
                              accepts multiple values
        zorder                use zorder to display overlapping domains. default
                              is ``2``
        ===================== ==================================================

    
    Available boxstyles_:
    
    * ``'larrow'``
    * ``'rarrow'``
    * ``'round'``
    * ``'round4'`` 
    * ``'roundtooth'`` 
    * ``'sawtooth'``
    * ``'square'``
    
    with examples_.
    
    Usage eg.
    
    ::
    
        feat = features.DomainFeature([seqfeature,seqfeature2], 
               name = ['test domain 1', 'test domain 2'], 
               boxstyle = ['sawtooth, tooth_size=1.4',  
                           'rarrow, pad = 0.1'],  
               seq_line = (1, 766),)

.. _boxstyles: http://matplotlib.sourceforge.net/api/artist_api.html?highlight\
=fancybboxpatch#matplotlib.patches.FancyBboxPatch

.. _examples: http://matplotlib.sourceforge.net/examples/pylab_examples/\
fancybox_demo2.html

    '''
    def __init__(self, domains, **kwargs):
        '''

        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        if isinstance(domains, list):
            self.domains = domains
        else:
            raise ValueError('needs a list of SeqFeature')
        self.name=kwargs.get('name','')
        self.seq_line_start, self.seq_line_end  = kwargs.get('seq_line',(0,0)) # this requires a tuple of start and stop positions
        self.zorder = kwargs.get('zorder',2)
        Xs = []
        for feature in self.domains:
            Xs.append(feature.location.start.position)
            Xs.append(feature.location.end.position)
        
        if self.seq_line_start:
            Xs.append(self.seq_line_start)
        if self.seq_line_end:
            Xs.append(self.seq_line_end)
        self.start = min(Xs)
        self.end = max(Xs)
        
        
            
    def draw_feature(self):
        if self.seq_line_start:
            line_y = self.Y + self.height/2
            feat_draw = plt.plot((self.seq_line_start, self.seq_line_end),
                                 (line_y, line_y), 
                                 linestyle = '-', 
                                 color = '#CCCCCC', 
                                 lw = 2,
                                 aa = False, 
                                 alpha = 1.0,
                                 zorder = 1,
                                 url = self.url,)
            self.patches.extend(feat_draw)
            
        for i,feature in enumerate(self.domains):
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            width = float(end - start)
            if isinstance(self.boxstyle, list):
                boxstyle = self.boxstyle[i]
            else:
                boxstyle = self.boxstyle
            if isinstance(self.ec, list):
                ec = self.ec[i]
            else:
                ec = self.ec
            if isinstance(self.ec, list):
                ec = self.ec[i]
            else:
                ec = self.ec
            if isinstance(self.fc, list):
                fc = self.fc[i]
            else:
                fc = self.fc
            if isinstance(self.alpha, list):
                alpha = self.alpha[i]
            else:
                alpha = self.alpha
            if isinstance(self.zorder, list):
                zorder = self.zorder[i]
            else:
                zorder = self.zorder
            feat_draw=FancyBboxPatch((start,self.Y), 
                                     width=(end-start),
                                     height=self.height, 
                                     boxstyle=boxstyle, 
                                     lw=self.lw,
                                     ec=ec, 
                                     fc=fc,
                                     alpha=alpha,
                                     mutation_aspect = 0.3,
                                     zorder = zorder,
                                     url = self.url,)
            self.patches.append(feat_draw)

            
    def draw_feat_name(self,**kwargs):
        self.feat_name=[]
        font_feat=FontProperties()
        set_size=kwargs.get('set_size','xx-small')
        set_family=kwargs.get('set_family','serif')
        set_weight=kwargs.get('set_weight','bold')
        va=kwargs.get('va','center')
        ha=kwargs.get('ha','center')
        font_feat.set_size(set_size)
        font_feat.set_family(set_family)       
        font_feat.set_weight(set_weight)  
        
        for i,feature in enumerate(self.domains):
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            width = float(end - start)
            if isinstance(self.name, list):
                name = self.name[i]
            else:
                name = self.name
            self.feat_name.append(plt.annotate(name, xy = (start + (end - start) / 2, self.Y), xytext = (start + (end - start) / 2, self.Y + self.height/2.), fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va))

