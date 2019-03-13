from setuptools import setup, find_packages
import os

version = '1.1'

setup(name='biograpy',
      version=version,
      description="",
      long_description=open("README.md").read() + "\n" +
                       open(os.path.join("docs", "HISTORY.txt")).read() + "\n" +
                       open(os.path.join("docs", "AUTHORS.txt")).read(),
      # Get more strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
        "Programming Language :: Python",
        ],
      keywords='',
      author='Andrea Pierleoni',
      author_email='apierleoni.dev@gmail.com',
      url='',
      license='LGPL',
	packages =['biograpy'],
#      packages=find_packages(exclude=['ez_setup', 'docs', 'tests*']),
#      package_dir={'':'biograpy'},
      namespace_packages=['biograpy'],
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          'setuptools',
          # -*- Extra requirements: -*-
          'biopython',
          'matplotlib>=2.0.0',
          'numpy>=1.1',
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
