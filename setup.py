from setuptools import setup, find_packages
import os

version = '1.0a'

setup(name='biograpy',
      version=version,
      description="",
      long_description=open("README.txt").read() + "\n" +
                       open(os.path.join("src", "biograpy",  "docs", "HISTORY.txt")).read() + "\n" +
                       open(os.path.join("src", "biograpy",  "docs", "AUTHORS.txt")).read(),
      # Get more strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
        "Programming Language :: Python",
        ],
      keywords='',
      author='Andrea Pierleoni',
      author_email='andrea@biocomp.unibo.it',
      url='',
      license='LGPL',
      packages=find_packages('src', exclude=['ez_setup']),
      package_dir = {'':'src'},
      namespace_packages=['biograpy'],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'setuptools',
          # -*- Extra requirements: -*-
          'biopython',
          'matplotlib>=1.0',
          'numpy>=1.1',
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
