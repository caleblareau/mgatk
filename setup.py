"""
mgatk: a mitochondrial genome analysis toolkit
"""
from setuptools import find_packages, setup

dependencies = ['click', 'pysam', 'pytest', 'snakemake', 'biopython', 'numpy', 'pandas', 'optparse-pretty', 'regex', 'ruamel.yaml']

setup(
    name='mgatk',
    version='0.6.1',
    url='https://github.com/caleblareau/mgatk',
    license='MIT',
    author='Caleb Lareau',
    author_email='caleb.lareau@gmail.com',
    description='Mitochondrial genome analysis toolkit.',
    long_description=__doc__,
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=dependencies,
    entry_points={
        'console_scripts': [
            'mgatk = mgatk.cli:main',
            'mgatk-del-find = mgatk.delete.clifind:main',
            'mgatk-del = mgatk.delete.clidel:main'
        ],
    },
    classifiers=[
        # As from http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'Development Status :: 1 - Planning',
        # 'Development Status :: 2 - Pre-Alpha',
        # 'Development Status :: 3 - Alpha',
        # 'Development Status :: 4 - Beta',
         'Development Status :: 5 - Production/Stable',
        # 'Development Status :: 6 - Mature',
        # 'Development Status :: 7 - Inactive',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
)
