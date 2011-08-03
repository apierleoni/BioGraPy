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
from matplotlib.patches import Rectangle,Circle, Wedge, Polygon,FancyBboxPatch, FancyArrow
from matplotlib.colors import Normalize,LogNorm,ListedColormap,LinearSegmentedColormap


class BaseGraphicFeature(object):
    ''' 
::

    base class for all GraphicFeature types
    this class template misses the "draw_feature" method that has to be
    implemented in the full functional classes

    self.draw_feature() must return one or more matplotlib patches
    self.start  must return the lower extent of the feature in the coordinate system
    self.end must return the higher extent of the feature in the coordinate system

    "__init__" method must be expanded accordingly to GraphicFeature input
    "draw_feat_name" method can be overridden, but must set a list of
                     plt.text object in self.feat_name
                     
    Optional values:
            name --> text that will be drawn under the GraphicFeature
            type --> define the GrafhFeature type, is needed by Drawer for feature grouping
            score --> this value can be used to color the feature accordingly
            height ---> GraphicFeature y axis extension
            cm ---> matplotlib colormap to use, custom colormap can be easily created but ar not supported right now. a list of default colormaps can be found at: http://www.astro.princeton.edu/~msshin/science/code/matplotlib_cm/
            color_by_cm ---> if True set the color by using the colormap (True by default)
            cm_value ---> value between 0 and 1 to use for chosing the color from the colormap (default random)
            use_score_for_color ---> if True use the score value to pick a color from the colormap. if scores are not between 0 and 1 they can be normalized by using Normalize and LogNorm from matplotlib.colors 
            ec ---> define edgecolor,  overrides color_by_cm and use_score_for_color
            fc ---> define facecolor,  overrides color_by_cm and use_score_for_color
            lw ---> defines edge line width
            alpha ---> defines facecolor aplha value
            boxstyle ---> FancyBboxPatch input string. take a look at http://matplotlib.sourceforge.net/api/artist_api.html?highlight=fancybboxpatch#matplotlib.patches.FancyBboxPatch
    '''
    
    default_cm='winter'#deafult matplotlib colormap to use
    
    def __init__(self, **kwargs):
        '''
        init
        '''
        self.name=kwargs.get('name','')
        self.type=kwargs.get('type','')
        self.score=kwargs.get('score',0.)
        self.height=kwargs.get('height',1.)
        # html map options
        self.url =  kwargs.get('url','')
        self.html_map_extend = kwargs.get('html_map_extend','') #included text will be added in html area tab, use onmouseover, etc... must be valid html attrib
        #feature style options
        self.cm = cm.get_cmap(kwargs.get('cm', self.default_cm))
        self.color_by_cm = kwargs.get('color_by_cm',True)
        self.cm_value = kwargs.get('cm_value',None)
        self.use_score_for_color=kwargs.get('use_score_for_color',False)
        if self.use_score_for_color:
            self.color_by_cm = True
        self.fc = kwargs.get('fc', self.cm(self.cm_value or 0))
        self.ec = kwargs.get('ec', None)
        self.lw = kwargs.get('lw',0.5)
        self.ls = kwargs.get('ls','-')
        self.alpha=kwargs.get('alpha',.8)
        self.boxstyle=kwargs.get('boxstyle','square, pad=0.')
        self.Y=0.0
        self.patches = [] # all the patches must be returned inside this list
        self.feat_name= [] # all the patches labels must be returned inside this list
        self.norm = kwargs.get('norm', None)
        

    def draw_feat_name(self,**kwargs):
        '''recalling self.draw_feat_name() will overwrite the previous feature name stored in self.feat_name '''
        self.feat_name = []
        font_feat = FontProperties()
        set_size = kwargs.get('set_size','x-small')
        set_family = kwargs.get('set_family','serif')
        set_weight = kwargs.get('set_weight','normal')
        va=kwargs.get('va','top')
        ha=kwargs.get('ha','left')
        xoffset = kwargs.get('xoffset',0)# to be used to move the text annotation along x axis
        font_feat.set_size(set_size)
        font_feat.set_family(set_family)
        font_feat.set_weight(set_weight)
        text_x = self.start
        if (self.start < xoffset) and (xoffset <= self.end):
            text_x+= xoffset
        self.feat_name = [plt.annotate(self.name, xy = (text_x, self.Y), xytext = (text_x, self.Y - self.height/5.), fontproperties=font_feat, horizontalalignment=ha, verticalalignment=va)]

        
class Simple(BaseGraphicFeature):
    '''
::    
    
    Handle a feature not derived from a SeqFeature. Just need a Start and Stop position
    Minimum definition are start and end positions.
        Eg.
        GraphicFeature.Simple(start,end,name='factor 7',**kwargs)
        Optional values:
        
            name --> text that will be drawn under the GraphicFeature
            type --> define the GrafhFeature type, is needed by Drawer for feature grouping
            score --> this value can be used to color the feature accordingly
            height ---> GraphicFeature y axis extension
            cm ---> matplotlib colormap to use, custom colormap can be easily created but ar not supported right now. a list of default colormaps can be found at: http://www.astro.princeton.edu/~msshin/science/code/matplotlib_cm/
            color_by_cm ---> if True set the color by using the colormap (True by default)
            cm_value ---> value between 0 and 1 to use for chosing the color from the colormap (default random)
            use_score_for_color ---> if True use the score value to pick a color from the colormap. if scores are not between 0 and 1 they can be normalized by using Normalize and LogNorm from matplotlib.colors 
            ec ---> define edgecolor,  overrides color_by_cm and use_score_for_color
            ec ---> define facecolor,  overrides color_by_cm and use_score_for_color
            lw ---> defines edge line width
            alpha ---> defines facecolor aplha value
            boxstyle ---> FancyBboxPatch input string. take a look at http://matplotlib.sourceforge.net/api/artist_api.html?highlight=fancybboxpatch#matplotlib.patches.FancyBboxPatch
        
        returns
        
            self.patches ---> list of drawn patches
            self.feat_name ---> matplotlib Text object reporting feature name, drawn under the patches
            self.type ---> feature type used by Drawer
        

    '''
    def __init__(self,start,end,**kwargs):
        '''
        init
        '''
        BaseGraphicFeature.__init__(self, **kwargs)
        self.start = start
        self.end = end



    def draw_feature(self):        
        feat_draw=FancyBboxPatch((self.start,self.Y), 
                                  width=(self.end-self.start),
                                  height=self.height, 
                                  boxstyle=self.boxstyle, 
                                  lw=self.lw,
                                  ec=self.ec,
                                  fc=self.fc, 
                                  alpha=self.alpha, 
                                  url = self.url, 
                                  mutation_aspect = 0.3)
        self.patches.append(feat_draw)


class GenericSeqFeature(BaseGraphicFeature):
    '''
::    
    
    Draws any seqfeature as a simple rectangular
    
    unlike GraphicFeature.Simple requires a SeqFeature object in input, while start and end can be automatically detected from seqfeature.
    GraphicFeature.Generic(feature,name='factor 7',start=10,end=143,score=0.2)
    '''
    def __init__(self,feature,**kwargs):
        '''

        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        self.start=kwargs.get('start',min([feature.location.start.position,feature.location.end.position]))
        self.end=kwargs.get('end',max([feature.location.start.position,feature.location.end.position]))
        self.type=kwargs.get('type',feature.type)
        if 'score' in feature.qualifiers:
            self.score=kwargs.get('score',feature.qualifiers['score'])


    def draw_feature(self):
        feat_draw=FancyBboxPatch((self.start,self.Y), width=(self.end-self.start),
            height=self.height, boxstyle=self.boxstyle, lw=self.lw,
            ec=self.ec, fc=self.fc,alpha=self.alpha, mutation_aspect = 0.3)
        self.patches.append(feat_draw)


class GeneSeqFeature(BaseGraphicFeature):
    '''
::    
    
    Draws a Gene Feature
    requires a SeqFeature (type='gene') and optionally a list of exons SeqFeatures that can be drawn over the gene patch
    
    '''
    def __init__(self,feature,exons=[],**kwargs):
        '''
        
        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        self.height=kwargs.get('height',2.)
        self.start=kwargs.get('start',min([feature.location.start.position,feature.location.end.position]))
        self.end=kwargs.get('end',max([feature.location.start.position,feature.location.end.position]))
        default_head_length = 10
        if abs(self.end -  self.start<= 50):
            default_head_length = (self.end-self.start)/10.
        self.head_length=kwargs.get('head_length',default_head_length)
        self.type=kwargs.get('type','Gene')
        self.feature=feature
        self.exons=exons


    def draw_feature(self):
        self.patches=[]
        if self.feature.strand==1:
            arrow_start=self.start
            arrow_direction=self.end-self.start
            shape='right'
            body_width=self.height*.6667
            head_width=self.height
        elif self.feature.strand==-1:
            arrow_start=self.end
            arrow_direction=-self.end-self.start
            shape='left'
            body_width=self.height*.6667
            head_width=self.height
        else:
            raise ValueError('Gene feature must have strand equal to 1 or -1')
        feat_draw=FancyArrow(arrow_start, self.Y, dx=arrow_direction, dy=0, ec=self.ec,
            fc=self.fc,alpha=self.alpha, width=body_width, head_length = self.head_length,
            head_width=head_width,lw=self.lw,length_includes_head=True,
            shape=shape, head_starts_at_zero=False)
        self.patches.append(feat_draw)
        for exon in self.exons:
            feat_draw=FancyBboxPatch((int(exon.location.start.position),self.Y),
                width=(int(exon.location.end.position)-int(exon.location.start.position)),
                height=body_width/2., boxstyle=self.boxstyle,lw=0, ec=self.ec,
                fc=self.fc,alpha=self.alpha+0.1,)
            self.patches.append(feat_draw)


class TextSequence(BaseGraphicFeature):
    '''
::    
    
    Draws a text biological sequence    
    '''
    def __init__(self,sequence,**kwargs):
        '''
        
        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        self.sequence = sequence
        self.start = kwargs.get('start',0)
        self.end = kwargs.get('end',len(self.sequence))
        self.set_size = kwargs.get('set_size','xx-small')
        self.set_family = kwargs.get('set_family','monospace')
        self.set_weight = kwargs.get('set_weight','normal')
        self.va = kwargs.get('va','bottom')
        self.ha = kwargs.get('ha','left')
        self.font_feat=FontProperties()
        self.font_feat.set_size(self.set_size)
        self.font_feat.set_family(self.set_family)       
        self.font_feat.set_weight(self.set_weight) 
        
    def draw_feature(self):
        self.patches=[]
        for i, char in enumerate(self.sequence):
            '''self.patches.append(plt.annotate(char, 
                                             xy = (self.start + i, 0), 
                                             xytext = (self.start + i, 0-self.height/5.), 
                                             fontproperties = self.font_feat, 
                                             horizontalalignment = self.ha, 
                                             verticalalignment = self.va ))'''
            self.patches.append(Text(x = self.start + i,
                                     y= 0,
                                     text = char, 
                                     fontproperties = self.font_feat, 
                                     horizontalalignment = self.ha, 
                                     verticalalignment = self.va ))
         
        

class PlotFeature(BaseGraphicFeature):
    '''
::    
    
    Draws a plot of contnuous value for a feature such as an hydrophobicity scale
    must be used within a PlotTrak to have axix.
    
**TO DO**


import pylab
from Bio import SeqIO
for subfigure in [1,2]:
    filename = "SRR001666_%i.fastq" % subfigure
    pylab.subplot(1, 2, subfigure)
    for i,record in enumerate(SeqIO.parse(filename, "fastq")):
        if i >= 50 : break #trick!
        pylab.plot(record.letter_annotations["phred_quality"])
    pylab.ylim(0,45)
    pylab.ylabel("PHRED quality score")
    pylab.xlabel("Position")
pylab.savefig("SRR001666.png")
print "Done"

    from Bio.SeqUtils import ProtParamData, ProtParam
    analysis = ProtParam.ProteinAnalysis('AKLSMCSKLAPERRALPRLALAPRALLRAPPAMMKMSKMKMKCVVCAAS')
    KDdata = analysis.protein_scale(ProtParamData.kd, 9, 0.4)
    
    input data  needs to be:
    [1,2,3], [1,2,3], 'go-', 
      y          x     style     
     
    y is necessary
    if no x is provided x will be created starting from position 1

    
    '''
    def __init__(self,  y , x = [], style = '-' , **kwargs):
        '''
        TO DO
        '''
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
                                    alpha = self.alpha)
            self.patches.extend(plotted_data)
            
        

class BarPlotFeature(BaseGraphicFeature):
    '''
::    
    
    Draws a bar plot values for a feature such as prediction confidence values
    must be used within a PlotTrak to have axix.
    

    '''
    def __init__(self,  y , x = [], style = '-' , **kwargs):
        '''
        TO DO
        
        color_by_cm will color based on the y values
        '''
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
                           )
                        
            if self.color_by_cm:
                self.norm = Normalize(min(self.y), max(self.y))
                for bar, value in zip(bars, self.y):
                    color = self.cm(self.norm(value))
                    bar.set_color(color)
                    
            self.patches.extend(bars)
            
        
        


class SegmentedSeqFeature(BaseGraphicFeature):
    '''
::    
    
    Draws a SeqFeature carrying 'join' subfeatures
    
    requires a SeqFeature object carryng 'joins' in SeqFeature.subfeatures, 
    start and end can be automatically detected from seqfeature.
    SegmentedSeqFeature(feature,name='factor 7',s)
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
                    feat_draw=FancyBboxPatch((int(sub_feature.location.start.position),self.Y), width=(int(sub_feature.location.end.position)-int(sub_feature.location.start.position)), height=self.height, boxstyle=self.boxstyle,lw=self.lw, ec=self.ec, fc=self.fc, alpha=self.alpha,)
                    self.patches.append(feat_draw)
                    if junction_start:
                        junction_end=float(sub_feature.location.start.position)
                        junction_middle=float((junction_start+junction_end)/2.)
                        Yends=self.Y+self.height/2.
                        Ymiddle=self.Y+self.height
                        join=plt.plot([junction_start,junction_middle,junction_end],[Yends,Ymiddle,Yends], lw=1 ,ls='-', c=self.ec, alpha=0.5)
                        self.patches.extend(join)
                    junction_start=float(sub_feature.location.end.position)
        else:#if SegmentedSeqFeature is called whihout subfeatures returns a GenericSeqFeature
            feat_draw=FancyBboxPatch((self.start,self.Y), width=(self.end-self.start), height=self.height, boxstyle=self.boxstyle, lw=self.lw, ec=self.ec, fc=self.fc,alpha=self.alpha,)
            self.patches.append(feat_draw)


class CoupledmRNAandCDS(BaseGraphicFeature):
    '''
::

    Draws mRNA with joins at a lower alpha and the corresponding CDS on top of it
    requires an mRNA seqfeature object and a CDS seqfeature object
    
    
    '''
    def __init__(self,mRNA,CDS,**kwargs):
        '''
        
        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        self.mRNA_start=kwargs.get('mRNA_start',min([mRNA.location.start.position,mRNA.location.end.position]))
        self.mRNA_end=kwargs.get('mRNA_end',max([mRNA.location.start.position,mRNA.location.end.position]))
        self.CDS_start=kwargs.get('CDS_start',min([CDS.location.start.position,CDS.location.end.position]))
        self.CDS_end=kwargs.get('CDS_end',max([CDS.location.start.position,CDS.location.end.position]))
        self.start = self.mRNA_start
        self.end = self.mRNA_end
        self.type=kwargs.get('type','mRNA')
        if 'score' in mRNA.qualifiers:
            self.score=kwargs.get('score',mRNA.qualifiers['score'])
        self.mRNA=mRNA
        self.CDS=CDS
        self.ec=kwargs.get('ec','k')


    def draw_feature(self):
        #draw mRNA
        if len(self.mRNA.sub_features):
            junction_start=False
            if self.mRNA.strand==-1:
                self.mRNA.sub_features.reverse()
            for sub_feature in self.mRNA.sub_features:
                if sub_feature.location_operator=='join':
                    feat_draw=FancyBboxPatch((int(sub_feature.location.start.position),self.Y), width=(int(sub_feature.location.end.position)-int(sub_feature.location.start.position)), height=self.height, boxstyle=self.boxstyle,lw=self.lw, ec=self.ec, fc=self.fc,alpha=self.alpha-0.5,)
                    self.patches.append(feat_draw)
                    if junction_start:
                        junction_end=float(sub_feature.location.start.position)
                        junction_middle=float((junction_start+junction_end)/2.)
                        Yends=self.Y+self.height/2.
                        Ymiddle=self.Y+self.height
                        join=plt.plot([junction_start,junction_middle,junction_end],[Yends,Ymiddle,Yends], lw=0.5 ,ls='-', c=self.ec,alpha=0.5)
                        self.patches.extend(join)
                    junction_start=float(sub_feature.location.end.position)
        else:#if SegmentedSeqFeature is called whihout subfeatures returns a GenericSeqFeature
            feat_draw=FancyBboxPatch((self.mRNA_start,self.Y), width=(self.mRNA_end-self.mRNA_start), height=self.height, boxstyle=self.boxstyle, lw=self.lw, ec=self.ec, fc=self.fc,alpha=self.alpha,)
            self.patches.append(feat_draw)
        #draw CDS
        if len(self.CDS.sub_features):
            junction_start=False
            if self.CDS.strand==-1:
                self.CDS.sub_features.reverse()
            for sub_feature in self.CDS.sub_features:
                if sub_feature.location_operator=='join':
                    feat_draw=FancyBboxPatch((int(sub_feature.location.start.position),self.Y), width=(int(sub_feature.location.end.position)-int(sub_feature.location.start.position)), height=self.height, boxstyle=self.boxstyle,lw=0, ec=self.ec, fc=self.fc,alpha=self.alpha,)
                    self.patches.append(feat_draw)
        else:#if SegmentedSeqFeature is called whihout subfeatures returns a GenericSeqFeature
            feat_draw=FancyBboxPatch((self.CDS_start,self.Y), width=(self.CDS_end-self.CDS_start), height=self.height, boxstyle=self.boxstyle, lw=self.lw, ec=self.ec, fc=self.fc,alpha=self.alpha,)
            self.patches.append(feat_draw)
        
            

            
            
class SinglePositionFeature(BaseGraphicFeature):
    '''
::    
    
    Draws features spanning just one position
    
        marker = 'o'
        markersize = 4
        
        
        available marker types:
        
            
        '.'     point marker
        ','     pixel marker
        'o'     circle marker
        'v'     triangle_down marker
        '^'     triangle_up marker
        '<'     triangle_left marker
        '>'     triangle_right marker
        '1'     tri_down marker
        '2'     tri_up marker
        '3'     tri_left marker
        '4'     tri_right marker
        's'     square marker
        'p'     pentagon marker
        '*'     star marker
        'h'     hexagon1 marker
        'H'     hexagon2 marker
        '+'     plus marker
        'x'     x marker
        'D'     diamond marker
        'd'     thin_diamond marker
        '|'     vline marker
        '_'     hline marker
        
    
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
        feat_draw=plt.plot(self.start, self.Y, marker=self.marker, markerfacecolor=self.fc, markeredgecolor='k', markersize=self.markersize, alpha=self.alpha,)
        self.patches=feat_draw

        
        
class TMFeature(BaseGraphicFeature):
    '''  
::    
    
    Draws TransMembrane helix with predicted cytoplasmic/non cytoplasmic sides

    TM = list of seqfeatures reporting Transmembrane regions 
        cyto =  list of seqfeatures reporting cytoplasmic regions
        non_cyto = list of seqfeatures reporting non cytoplasmic regions
        
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
        self.non_cyto_linestyle = kwargs.get('non_cyto_linestyle',':')
        self.non_cyto_label = kwargs.get('non_cyto_label','not cyto')
        self.fill_color = kwargs.get('fill_color',self.cyto_color)
        self.fill_linestyle = kwargs.get('fill_linestyle',self.cyto_linestyle)
        
        available linestyles:
        
        '-'     solid line style
        '--'     dashed line style
        '-.'     dash-dot line style
        ':'     dotted line style
        

    
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
        self.non_cyto_linestyle = kwargs.get('non_cyto_linestyle',':')
        self.non_cyto_label = kwargs.get('non_cyto_label','non cyto')
        self.fill_color = kwargs.get('fill_color',self.cyto_color)
        self.fill_linestyle = kwargs.get('fill_linestyle',self.cyto_linestyle)
        self.connection_lw = kwargs.get('connection_lw',2)
        
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
        
        if self.cyto and self.non_cyto:
            self.fill = False 
        else:
            self.fill = True #draw line between TM
            
        #label positions
        self.TM_starts = []
        self.cyto_starts = []
        self.non_cyto_starts = []
        
        

    
    
    def draw_feature(self):
        'draw TM regions'
        
        def drow_orizontal_line(start, end, y, linestyle, color):
            return plt.plot((start,end), (y,y), linestyle=linestyle, color=color, 
                            lw=self.connection_lw, aa=False, alpha=self.alpha,)
        
        TMs=[]
        for feature in self.TM:
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            type = feature.type
            self.TM_starts.append((start,end))
            if 'score' in feature.qualifiers:#not used yet
                score=kwargs.get('score',feature.qualifiers['score'])
            feat_draw=FancyBboxPatch((start,self.Y), width=(end-start), height=self.height, boxstyle=self.boxstyle, lw=self.lw, ec=self.TM_ec, fc=self.TM_fc,alpha=self.alpha,)
            self.patches.append(feat_draw)
        if self.fill:
            start = 1
            for TM_start,TM_end in sorted(self.TM_starts):
                end = TM_start# - 1
                feat_draw=drow_orizontal_line(start, end, self.Y + self.height/2, self.fill_linestyle, self.fill_color)
                self.patches.extend(feat_draw)#use extend, since plot returns a list of Line2D objects
                start = TM_end #+ 1
        'draw cytplasmic regions'
        for feature in self.cyto:
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            self.cyto_starts.append((start,end))
            feat_draw=drow_orizontal_line(start, end, self.Y + self.height/2, self.cyto_linestyle, self.cyto_color)
            self.patches.extend(feat_draw)
        'draw non cytplasmic regions'
        for feature in self.non_cyto:
            start = min([feature.location.start.position,feature.location.end.position])
            end = max([feature.location.start.position,feature.location.end.position])
            self.non_cyto_starts.append((start,end))
            feat_draw=drow_orizontal_line(start, end, self.Y + self.height/2, self.non_cyto_linestyle, self.non_cyto_color)
            self.patches.extend(feat_draw)
                 
                
    def draw_feat_name(self,**kwargs):
        '''recalling self.draw_feat_name() will overwrite the previous feature name stored in self.feat_name '''
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
::    
    
    Draws secondary structures
    
        self.betas = kwargs.get('betas',[])
        self.alphah = kwargs.get('alphah',[])
        self.coil = kwargs.get('coil',[])
        self.type=kwargs.get('type','')
        
        self.betas_ec = kwargs.get('betas_ec','k')
        self.betas_fc = kwargs.get('betas_fc','yellow')
        self.alphah_ec = kwargs.get('alphah_ec','k')
        self.alphah_fc = kwargs.get('alphah_fc','magenta')
        self.coil_color = kwargs.get('coil_color','k')
        self.coil_linestyle = kwargs.get('coil_linestyle','-')
        self.lw = kwargs.get('lw',0.5)
        
        self.filter_struct_length = kwargs.get('filter_struct_length',3)#if != 0 smooths secondary structures, by displaying only those longher than the given value. self.filter_struct_length force coil geneation, and ignore supplied coils regions
        
    
    '''
    def __init__(self,**kwargs):
        '''

        '''
        BaseGraphicFeature.__init__(self,**kwargs)
        
        self.betas = kwargs.get('betas',[])
        self.alphah = kwargs.get('alphah',[])
        self.coil = kwargs.get('coil',[])
        self.type=kwargs.get('type','')
        
        self.betas_ec = kwargs.get('betas_ec','k')
        self.betas_fc = kwargs.get('betas_fc','yellow')
        self.alphah_ec = kwargs.get('alphah_ec','k')
        self.alphah_fc = kwargs.get('alphah_fc','magenta')
        self.coil_color = kwargs.get('coil_color','k')
        self.coil_linestyle = kwargs.get('coil_linestyle','-')
        self.lw = kwargs.get('lw',2)
        
        self.filter_struct_length = kwargs.get('filter_struct_length',3)#if != 0 smooths secondary structures, by displaying only those longher than the given value. self.filter_struct_length force coil geneation, and ignore supplied coils regions
        
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
        def drow_orizontal_line(start, end, y, linestyle, color):
            return plt.plot((start,end), (y,y), linestyle=linestyle, color=color, 
                                lw=1, aa=False, alpha=self.alpha,)
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
                                       lw=self.lw, length_includes_head=True,  head_starts_at_zero=False)
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
                                           lw=self.lw, ec=self.alphah_ec, fc=self.alphah_fc,alpha=self.alpha,)
                self.patches.append(feat_draw)
        
        if (not self.coil) or (self.filter_struct_length):
            start = 1
            for region_start,region_end in sorted(self.structured_regions):
                end = region_start #- 1
                feat_draw=drow_orizontal_line(start, end, self.Y + self.height/2, self.coil_linestyle, self.coil_color)
                self.patches.extend(feat_draw)#use extend, since plot returns a list of Line2D objects
                start = region_end #+ 1
        else:
             for feature in self.coil:
                start = min([feature.location.start.position,feature.location.end.position])
                end = max([feature.location.start.position,feature.location.end.position])
                type = feature.type
                if 'score' in feature.qualifiers:#not used yet
                    score=kwargs.get('score',feature.qualifiers['score'])
                feat_draw=drow_orizontal_line(start, end, self.Y + self.height/2, self.coil_linestyle, self.coil_color)
                self.patches.extend(feat_draw)#use extend, since plot returns a list of Line2D objects
                
                
class DomainFeature(BaseGraphicFeature):
    '''
::    
    
    Draws domains(encoded as seqfeatures) with FancyBboxPatch
    can optionally drow a a line representing the complete sequence
    
    domains = accepts a list of SeqFeature objects
    self.boxstyle is inherited by  BaseGraphicFeature ---> in this class accepts one or a list of matplotlib FancyBboxPatch boxstyle
    
    available boxstyle: (http://matplotlib.sourceforge.net/api/artist_api.html?highlight=fancybboxpatch#matplotlib.patches.FancyBboxPatch)
    larrow 
    rarrow 
    round 
    round4 
    roundtooth 
    sawtooth
    square
    
    
    style examples at:
    http://matplotlib.sourceforge.net/examples/pylab_examples/fancybox_demo2.html
    
    self.seq_line_start, self.seq_line_end = kwargs.get('self.seq_line',(0,0))
    
    self.boxstyle, self.ec, self.fc, self.alpha, self.zorder can be a single value (applied to all domain features) 
    or a list of values (each one corresponding to a  domain features)

    use zorder to display overlapping domains
        
**DRAFT**

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
        '''draw domains '''
        if self.seq_line_start:
            line_y = self.Y + self.height/2
            feat_draw = plt.plot((self.seq_line_start, self.seq_line_end),
                                 (line_y, line_y), 
                                 linestyle = '-', 
                                 color = '#CCCCCC', 
                                 lw = 2,
                                 aa = False, 
                                 alpha = 1.0,
                                 zorder = 1)
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
                                     zorder = zorder)
            self.patches.append(feat_draw)

            
    def draw_feat_name(self,**kwargs):
        '''recalling self.draw_feat_name() will overwrite the previous feature name stored in self.feat_name '''
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

