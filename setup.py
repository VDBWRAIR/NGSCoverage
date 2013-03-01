import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description. It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def scripts( ):
    return [os.path.join( 'bin', f ) for f in os.listdir( 'bin' )]

# The next three lines are modified from Biopython
__version__ = "Undefined"
for line in open('coverage/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())
        break

setup(
    name = "coverage",
    version = __version__,
    author = "Tyghe Vallard",
    author_email = "vallardt@gmail.com",
    description = ("Coverage for NGS sequencing runs"),
    keywords = "biopython walter reed research python coverage gaps",
    url = "https://github.com/VDBWRAIR/coverage",
    packages = ['coverage'],
    scripts = scripts(),
    data_files = [
    ],
    install_requires = [
        "numpy >=1.6",
        "biopython >=1.59",
        "wrairlib >= 0.0.5",
        'matplotlib',
    ],
    long_description=read('README.md'),
)
