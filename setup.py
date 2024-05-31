#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = [ ]

setup(
    author="Yasser Alemán Gómez",
    author_email='yasseraleman@protonmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="The thalamic nuclei are involved in many neurodegenerative diseases and therefore, their identification is of key importance in numerous clinical treatments. Automated segmentation of thalamic subparts is currently achieved by exploring diffusion-weighted magnetic resonance imaging (DW-MRI), but in absence of such data, atlas-based segmentation can be used as an alternative. Currently, there is a limited number of available digital atlases of the thalamus. Moreover, all atlases are created using a few subjects only, thus are prone to errors due to the inter-subject variability of the thalamic morphology. In this work, we present a probabilistic atlas of anatomical subparts of the thalamus built upon a relatively large dataset where the individual thalamic parcellation was done by employing a recently proposed automatic diffusion-based clustering method. Our analyses, comparing the segmentation performance between the atlas-based and the clustering method, demonstrate the ability of the provided atlas to substitute the automated diffusion-based subdivision in the individual space when the DW-MRI is not available.",
    entry_points={
        'console_scripts': [
            'mialthalnucparc=mialthalnucparc.cli:main',
        ],
    },
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='mialthalnucparc',
    name='mialthalnucparc',
    packages=find_packages(include=['mialthalnucparc', 'mialthalnucparc.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/yaleman/mialthalnucparc',
    version='0.1.0',
    zip_safe=False,
)
