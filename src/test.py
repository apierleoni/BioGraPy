from cStringIO import StringIO
from biograpy.seqrecord import SeqRecordDrawer
from Bio import SeqIO
import cProfile

def test():
    seqrecord = SeqIO.parse(open("biograpy/tests/P31946.gbk","r"), "genbank").next()
    imgdata=StringIO()
    SeqRecordDrawer(seqrecord, fig_width=1500).save(imgdata, format='PNG')
    open("P31946.png", "w").write(imgdata.getvalue())
    
def test2():
    seqrecord = SeqIO.parse(open("biograpy/tests/test_uniprot.xml","r"), "uniprot").next()
    imgdata=StringIO()
    SeqRecordDrawer(seqrecord, fig_width=1500).save(imgdata, format='PNG')
    open("biograpy/tests/uniprot.png", "w").write(imgdata.getvalue())

cProfile.run("test2()", "test.profile")

import pstats
stat = pstats.Stats('test.profile')
#stat.sort_stats('time').print_stats(10)
stat.sort_stats('time').print_stats('drawer', 10)
#stat.sort_stats('calls').print_stats()
