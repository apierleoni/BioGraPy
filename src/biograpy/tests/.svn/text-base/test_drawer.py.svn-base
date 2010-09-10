import unittest
import tempfile
import Image
from plone4bio.graphics import drawer
from plone4bio.graphics import features

class TestDrawer(unittest.TestCase):
    def test_simple_features(self):
        panel = drawer.Panel(2000,fig_width=1000,fig_dpi=100)#initialize a drawer
        panel.draw_ref_obj(10,950,name='reference',height=0.02)#draw a reference object
        panel.add_feature(features.Simple(name='feat1',start=100,end=756,fc='r',aplha=0.5,height=0.015))
        panel.add_feature(features.Simple(name='feat2',start=300,end=1056,fc='pink',aplha=0.5,height=0.015))
        panel.add_feature(features.Simple(name='feat3',start=600,end=1356,fc='y',aplha=0.5,height=0.015))
        panel.add_feature(features.Simple(name='feat4  ',start=800,end=1356,fc='g',aplha=0.5,height=0.015))
        panel.add_feature(features.Simple(name='feat5',start= 1357,end=1806,fc='b',aplha=0.5,height=0.015))
        fh = tempfile.TemporaryFile()
        panel.save(fh, format='png')
        fh.seek(0)
        img = Image.open(fh)
        self.assertEqual(img.size, (1000, 220))
        self.assertEqual(img.format, 'PNG')

def test_suite():
    return unittest.makeSuite(TestDrawer)

if __name__ == '__main__':
    unittest.main(defaultTest='test_suite')

