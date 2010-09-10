import unittest
import doctest

def test_suite():
    return unittest.TestSuite([
        doctest.DocFileSuite('README.txt', package='plone4bio.graphics'),
        ])

if __name__ == '__main__':
    unittest.main(defaultTest='test_suite')

