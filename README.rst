=============================================================================================
Atlas-based thalamic nuclei parcellation using the atlas developed by Najdenovska et al. 2018
=============================================================================================


.. image:: https://img.shields.io/pypi/v/mialthalnucparc.svg
        :target: https://pypi.python.org/pypi/mialthalnucparc

.. image:: https://img.shields.io/travis/yaleman/mialthalnucparc.svg
        :target: https://travis-ci.com/yaleman/mialthalnucparc

.. image:: https://readthedocs.org/projects/mialthalnucparc/badge/?version=latest
        :target: https://mialthalnucparc.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




The thalamic nuclei are involved in many neurodegenerative diseases and therefore, their identification is of key importance in numerous clinical treatments. Automated segmentation of thalamic subparts is currently achieved by exploring diffusion-weighted magnetic resonance imaging (DW-MRI), but in absence of such data, atlas-based segmentation can be used as an alternative. Currently, there is a limited number of available digital atlases of the thalamus. Moreover, all atlases are created using a few subjects only, thus are prone to errors due to the inter-subject variability of the thalamic morphology. In this work, we present a probabilistic atlas of anatomical subparts of the thalamus built upon a relatively large dataset where the individual thalamic parcellation was done by employing a recently proposed automatic diffusion-based clustering method. Our analyses, comparing the segmentation performance between the atlas-based and the clustering method, demonstrate the ability of the provided atlas to substitute the automated diffusion-based subdivision in the individual space when the DW-MRI is not available.


* Free software: Apache Software License 2.0
* Documentation: https://mialthalnucparc.readthedocs.io.


Features
--------

* TODO

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
