from plone4bio.graphics import Panel, tracks, features

panel=Panel(fig_width=2000, fig_dpi=100, join_tracks = False, grid = 'both')#initialize a drawer



test_track = tracks.BaseTrack(name = 'test',sort_by= None)


feat1=features.Simple(name='feat1',start=100,end=356,fc='r', alpha =0.5, color_by_cm = False)
feat2=features.Simple(name='feat2',start=30,end=856,fc='pink', alpha =0.5, color_by_cm = False)
feat3=features.Simple(name='feat3',start=400,end=1356,fc='y', alpha =0.5,color_by_cm = False)
feat4=features.Simple(name='feat4',start=800,end=1256,fc='g', alpha =0.5,color_by_cm = False)
feat5=features.Simple(name='feat5',start= 1367,end=1806,fc='b', alpha =0.5,color_by_cm = False)

test_track.append(feat1)
test_track.append(feat2)
test_track.append(feat3)
test_track.append(feat4)
test_track.append(feat5)

panel.add_track(test_track)


test_track2 = tracks.BaseTrack(features.Simple(name='feat1',start=100,end=356, alpha = 0.75, color_by_cm = True),
                               features.Simple(name='feat2',start=30,end=856, alpha = 0.75,  color_by_cm = True),
                               features.Simple(name='feat3',start=400,end=1356, alpha = 0.75, color_by_cm = True),
                               features.Simple(name='feat4',start=800,end=1256, alpha = 0.75, color_by_cm = True),
                               features.Simple(name='feat5',start= 1367,end=1806, alpha = 0.75, color_by_cm = True),
                               sort_by = 'length', 
                               draw_cb = True)



panel.add_track(test_track2)

test_track3 = tracks.BaseTrack(name = 'test3', cm = 'jet', draw_cb = True)

feat1=features.Simple(name='feat1',start=100,end=356,alpha = 0.75, color_by_cm = True)
feat2=features.Simple(name='feat2',start=30,end=856,alpha = 0.75, color_by_cm = True)
feat3=features.Simple(name='feat3',start=400,end=1356,alpha = 0.75, color_by_cm = True)
feat4=features.Simple(name='feat4',start=800,end=1256,alpha = 0.75, color_by_cm = True)
feat5=features.Simple(name='feat5',start= 1367,end=1806,alpha = 0.75, color_by_cm = True)


test_track3.extend([feat1, feat2, feat3, feat4, feat5])

panel.add_track(test_track3)

'''SinglePositionFeature'''
from Bio import SeqFeature

test_track4 = tracks.BaseTrack(name = 'test single positions',
                               sort_by = 'collapse', 
                               cm = 'Pastel2',
                               draw_cb = True)
for pos in [50, 60, 100, 300, 500, 699, 1200, 1201, 1215,1280, 1400, 1600]:
    feat = SeqFeature.SeqFeature()
    feat.location = SeqFeature.FeatureLocation(pos, pos)
    test_track4.append(features.SinglePositionFeature(feat, name = str(pos), alpha = 0.75, color_by_cm = True))
panel.add_track(test_track4)


panel.save('../docs/images/test.png')

#print panel.htmlmap




''' draw seqrecord '''

from plone4bio.graphics import Panel, tracks, features
from plone4bio.graphics.seqrecord import SeqRecordDrawer
from Bio import SeqIO

seqrec = SeqIO.read(open('test_uniprot.xml'),'uniprot')
grseqrec = SeqRecordDrawer(seqrec)
grseqrec.save('../docs/images/test_uniprot.png')




from plone4bio.graphics import Panel, tracks, features
from plone4bio.graphics.seqrecord import SeqRecordDrawer
from Bio import SeqIO

seqrec = SeqIO.read(open('factor7.embl'),'embl')
grseqrec = SeqRecordDrawer(seqrec)
grseqrec.save('../docs/images/factor7.png')


