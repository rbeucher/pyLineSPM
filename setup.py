from setuptools import setup, find_packages
from os import path

MAJOR = 0
MINOR = 1
MICRO = 0
ISRELEASED = False
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='pyLineSPM',
    setup_requires=[
        'setuptools>=18.0',
        ],
    version=VERSION,
    description='Erosion on a line',
    long_description=long_description,
    url='https://github.com/rbeucher/pyLineSPM.git',
    author='Romain Beucher',
    author_email='romain.beucher@anu.edu.au',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',

        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Build Tools',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords=['erosion', "surface processes", "SPM"],
    install_requires=requirements,
)
