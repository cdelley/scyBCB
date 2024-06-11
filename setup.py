from setuptools import setup

setup(
    name='scyBCB',
    version='0.5.0',    
    description="""A Python package to analyze single-cell sequencing experiments (DNA, Protein, DAb-seq, etc.) with custom barcode beads""",
    url='https://github.com/cdelley/scyBCB',
    author='Cyrille L. Delley',
    author_email='cyrille.delley@ucsf.edu',
    license='MIT',
    packages=['scyBCB'],
    install_requires=[
        'numpy', 
        'Levenshtein',
    ],

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3.12',
    ],
)

