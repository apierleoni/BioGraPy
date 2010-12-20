Simple test

    >>> import tempfile
    >>> import os
    >>> import Image
    >>> from Bio import SeqIO

    >>> from biograpy import Panel, tracks, features

    >>> panel=Panel(2000,fig_width=1000,fig_dpi=100)#initialize a drawer
    >>> panel.draw_ref_obj(10,950,name='reference',height=0.02)#draw a reference object
    >>> test_track = tracks.BaseTrack('test') 

Create 5 graphicfeatures

    >>> feat1=features.Simple(name='feat1',start=100,end=756,fc='r',aplha=0.5,height=0.015)
    >>> feat2=features.Simple(name='feat2',start=300,end=1056,fc='pink',aplha=0.5,height=0.015)
    >>> feat3=features.Simple(name='feat3',start=600,end=1356,fc='y',aplha=0.5,height=0.015)
    >>> feat4=features.Simple(name='feat4',start=800,end=1356,fc='g',aplha=0.5,height=0.015)
    >>> feat5=features.Simple(name='feat5',start= 1357,end=1806,fc='b',aplha=0.5,height=0.015)

Add the features to the track::

    >>> test_track.add_feature(feat1)
    >>> test_track.add_feature(feat2)
    >>> test_track.add_feature(feat3)
    >>> test_track.add_feature(feat4)
    >>> test_track.add_feature(feat5)

Add the track to the panel::

    >>> panel.add_track(test_track)

Save the drawn image in PDF format::

    >>> (fhandle, fname) = tempfile.mkstemp(suffix='.pdf')
    >>> os.close(fhandle)
    >>> panel.save(fname)
    >>> open(fname).read(4)
    '%PDF'
    >>> os.unlink(fname)




from SeqRecord ...

>>> from biograpy import Panel, tracks, features
>>> from biograpy.seqrecord import SeqRecordDrawer
>>> from Bio import SeqIO
>>> seqrec = SeqIO.read(open('plone4bio/graphics/tests/test_uniprot.xml'),'uniprot')
>>> grseqrec = SeqRecordDrawer(seqrec)
>>> grseqrec.save('plone4bio/graphics/tests/test_uniprot.png')
