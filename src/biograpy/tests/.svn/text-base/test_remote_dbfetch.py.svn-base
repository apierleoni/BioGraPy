
from Bio import SeqIO 
import urllib 

def test_uniprot():
    from plone4bio.graphics import SeqRecordDrawer
    test_entry = 'S4A7_HUMAN'
    seqrec = SeqIO.read(urllib.urlopen('http://www.ebi.ac.uk/Tools/webservices/rest/dbfetch/uniprotkb/%s/uniprotxml'%test_entry), 'uniprot')
    handler = SeqRecordDrawer(seqrec,  fig_width=1000)
    handler.save(test_entry + '.png')
    
    
if __name__ == '__main__':
    test_uniprot()