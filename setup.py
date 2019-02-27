from distutils.core import setup

setup(
    name='vcfiterator',
    packages=['vcfiterator'],
    version='0.1.3',
    description='A package for iterating of .vcf files, parsing it to data structures. Especially useful for converting vcf files to JSON.',
    author='Svein Tore Koksrud Seljebotn',
    author_email='sveint@gmail.com',
    url='https://github.com/sveint/vcfiterator',
    license='MIT',
    keywords=['parsing', 'iterator', 'vcf', 'json'],
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
    ],

)
