
from Bio import SeqIO 
import urllib 

def test_uniprot():
    from biograpy.seqrecord import SeqRecordDrawer
    test_entry = 'S4A7_HUMAN'
    seqrec = SeqIO.read(urllib.urlopen('http://www.ebi.ac.uk/Tools/webservices/rest/dbfetch/uniprotkb/%s/uniprotxml'%test_entry), 'uniprot')
    handler = SeqRecordDrawer(seqrec,  fig_width=1000)
    fh = tempfile.TemporaryFile()
    handler.save(fh + '.png')
    fh.seek(0)
    img = Image.open(fh)
    self.assertEqual(img.format, 'PNG')

    
    
if __name__ == '__main__':
    test_uniprot()