from distutils.core import setup
from codecs import open
import os

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the relevant file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='vcfiterator',
    packages=['vcfiterator'],
    version='0.1-pre4',
    description='A package for iterating of .vcf files, parsing it to data structures. Especially useful for converting vcf files to JSON.',
    long_description=long_description,
    author='Svein Tore Koksrud Seljebotn',
    author_email='sveint@gmail.com',
    url='https://github.com/sveint/vcfiterator',
    license='MIT',
    keywords=['parsing', 'iterator', 'vcf', 'json'],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
    ],

)
