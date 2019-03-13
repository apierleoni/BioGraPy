BioGraPy - Biological Graphical Library in Python
=================================================

Fork that adds support to plot the tracks below an existing figure, and adds some new features to draw such as text sequence with automatic font size and highlighted text.

```python
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np

import biograpy
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq


fig1, ax = plt.subplots(figsize=(10,3), dpi=300)

ax.scatter(x=list(range(0, 1200)), y=np.random.random(1200))
ax.set_xlim([0,1200])
ax.set_ylim([0,2])

# We use the package biograpy to create Panel, tracks, and features
panel = biograpy.Panel(fig1)
track = biograpy.tracks.BaseTrack()

# Simple feature drawn as a rectangle
track.add_feature(biograpy.features.Simple(name='Simple1', start = 50, end = 300))

# Gene feature drawn as an arrow
genefeat = SeqFeature(FeatureLocation(100,500), type='gene', strand=1)
track.add_feature(biograpy.features.GeneSeqFeature(genefeat, name='GeneSeqFeature1'))

# Simple feature with color
track.add_feature(biograpy.features.Simple(name='Simple_colored1', start = 500, end = 820, color_by_cm=False, fc='r'))

# Gene feature with color
track.add_feature(biograpy.features.GeneSeqFeature(genefeat,name='GeneSeqFeature2_colored',fc='r'))

# Very short gene feature (to test the arrow head automatic sizing)
genefeat = SeqFeature(FeatureLocation(800,810), type='gene', strand=1)
track.add_feature(biograpy.features.GeneSeqFeature(genefeat, name='GeneSeqFeature_very_short'))

# mRNA with one CDS drawn as a shaded rectangle with an arrow on top
CDS_feature = SeqFeature(FeatureLocation(180,1000), type='CDS', strand=-1)
mRNA_feature = SeqFeature(FeatureLocation(180,1100), type='mRNA', strand=-1)
mRNAandCDSfeat = biograpy.features.CoupledmRNAandCDS(mRNA_feature, CDS_feature,
                                                     name='CoupledmRNAandCDS1_strandm', ec='k')
track.add_feature(mRNAandCDSfeat)

# mRNA with CDS, custom color
CDS_feature = SeqFeature( FeatureLocation(200,600), type='CDS', strand=1)
mRNA_feature = SeqFeature( FeatureLocation(100,800), type='mRNA', strand=1)
mRNAandCDSfeat = biograpy.features.CoupledmRNAandCDS(mRNA_feature, CDS_feature,
                                                     name='CoupledmRNAandCDS1_strandp_colored', ec='r', fc='r')
track.add_feature(mRNAandCDSfeat)

# Very short mRNA with CDS (to test the arrow head automatic sizing)
CDS_feature = SeqFeature( FeatureLocation(250,260), type='CDS', strand=1)
mRNA_feature = SeqFeature( FeatureLocation(220,300), type='mRNA', strand=1)
mRNAandCDSfeat = biograpy.features.CoupledmRNAandCDS(mRNA_feature, CDS_feature,
                                                     name='CoupledmRNAandCDS1_very_short')
track.add_feature(mRNAandCDSfeat)

# Gene feature on minus strand
geneSeqfeat = SeqFeature( FeatureLocation(900,200), type = 'gene', strand=-1)
genefeat = biograpy.features.GeneSeqFeature(geneSeqfeat, name='GeneSeqFeature3_strandm')
track.add_feature(genefeat)

track2 = biograpy.tracks.BaseTrack(biograpy.features.Simple(name='Simple3_track2', start= 0, end = 80))

panel.add_track(track)
panel.add_track(track2)
panel._draw_tracks()
panel.save('biograpy_test1.png')
```

![example1](/examples/biograpy_test1.png)

```python
fig1, ax = plt.subplots(figsize=(11,3), dpi=400)

ax.scatter(x=list(range(0, 100)), y=np.random.random(100))
ax.set_xlim([0,100])
ax.set_ylim([0,2])

# Create panel and track
panel = biograpy.Panel(fig1)
track = biograpy.tracks.BaseTrack()

# Text sequence feature
track.add_feature(biograpy.features.TextSequence(
        'MKKVIVIGVNHAGTSFIRTLLSKSKDFQVNAYDRNTNISFLGCGIALAVSGVVKNTEDLFYSTPEELKAMGANVFMAHDVVGLDLDKKQVIVKDL',
        start=10, name="TextSequence1"))

"""
Pretty text sequence feature allows to pass a text sequence together with a list of regions
to highlight. The font size is calculated automatically to fit the scale of the x axis.
"""

prettyTextFeat = biograpy.features.PrettyTextSequence(
    'KSKDFQVNAYDRNMKKVIVIGVNHAGTSFIRTLLSKSKDFQVNAYDRNTNISFLGCGIALAVSGVVKNTEDLFYSTPEELKAMGANVFMAHDVVGLDLDKKQVIVKDLATGKETVDHY',
    highlightList=[{'start':0, 'end':10, 'background_color':'red', 'background_color_alpha':0.4, 'weight':900},
                   {'start':20, 'end':25, 'foreground_color':'blue', 'weight':100}
                  ],
    start=20, name='PrettyTextSequence1')
track.add_feature(prettyTextFeat)

panel.add_track(track)
panel.save('biograpy_test2.png')
panel._draw_tracks()
```

![example1](/examples/biograpy_test2.png)