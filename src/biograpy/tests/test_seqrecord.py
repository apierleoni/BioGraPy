import unittest
import os
import tempfile
import Image
from Bio import SeqIO
from biograpy import Panel, tracks, features, seqrecord

class TestSeqrecord(unittest.TestCase):
    def setUp(self):
        self.emblpath = os.path.sep.join([os.path.dirname(__file__), "factor7.embl"])
        self.seqr=SeqIO.parse(open(self.emblpath), 'embl').next()

    def test_embl(self):
        handler = seqrecord.SeqRecordDrawer(self.seqr, fig_width=1000)
        (fhandle, fname) = tempfile.mkstemp(suffix='.png')
        os.close(fhandle)
        handler.save(fname)
        img = Image.open(fname)
        self.assertEqual(img.size, (1250, 587))
        self.assertEqual(img.format, 'PNG')
        os.unlink(fname)

    def test_embl_500x500(self):
        handler = seqrecord.SeqRecordDrawer(self.seqr, fig_width=500, fig_height=500)
        (fhandle, fname) = tempfile.mkstemp(suffix='.png')
        os.close(fhandle)
        handler.save(fname)
        img = Image.open(fname)
        self.assertEqual(img.size, (625, 625))
        self.assertEqual(img.format, 'PNG')
        os.unlink(fname)

    def test_patches(self):
        # TODO: test
        pass

    def test_keeporder(self):
        # TODO: test
        pass

def test_suite():
    return unittest.makeSuite(TestSeqrecord)

if __name__ == '__main__':
    unittest.main(defaultTest='test_suite')

