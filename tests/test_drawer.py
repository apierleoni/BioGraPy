import unittest
import tempfile
import Image
from biograpy import Panel, tracks, features

class TestDrawer(unittest.TestCase):
    def test_simple_features(self):
        panel=Panel(fig_width=1000)#initialize a drawer
        test_track = tracks.BaseTrack(features.Simple(name='feat1',start=100,end=756,),
            features.Simple(name='feat2',start=300,end=1056,),
            features.Simple(name='feat3',start=600,end=1356,),
            features.Simple(name='feat4',start=800,end=1356,),
            features.Simple(name='feat5',start= 1357,end=1806,),
            name = 'test')
        panel.add_track(test_track)
        fh = tempfile.TemporaryFile()
        panel.save(fh, format='png')
        fh.seek(0)
        img = Image.open(fh)
        #self.assertEqual(img.size, (1000, 220))
        self.assertEqual(img.format, 'PNG')

def test_suite():
    return unittest.makeSuite(TestDrawer)

if __name__ == '__main__':
    unittest.main(defaultTest='test_suite')

