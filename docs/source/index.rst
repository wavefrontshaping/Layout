
Introduction
============

Brief description
-----------------

SLMlayout is designed to help generating amplitude and/or phase masks using Spatial Light Modulators (SLMs)
with in mind applications of wavefrontshaping in complex media.
It allows creating patterns which consist of elementary segments on which we want to control the optical field.
This pattern are to be displayed with phase only SLMs or amplitude only Digital Micromirror Devices (DMDs).
A set of advanced features are available to use DMDs for phase or amplitude and phase modulation 
using spatial filtering in the Fourier plane.

SLMlayout is compatible with the `ALP4lib <https://github.com/wavefrontshaping/ALP4lib/>`_ 
module to control Vialux DMDs 
but can also be used for other types of SLMs or DMDs. 
In particular, it allows generating bitplanes in C format data that can be sent to the DMD with 
fast transfer time (see Vialux documentation). The corresponding function are coded in Cython 
for optimized computation time.


Who and why
-----------

This module was created and is maintained by Sébastien M. Popoff and improved by Maxime W. Matthès.
It was developped at the `Langevin Institute <https://www.institut-langevin.espci.fr/>`_
in the context of the ANR Project Molotof (ANR-16-CE25-0008-01 MOLOTOF).


Website
-------

Related discussions and codes are located on our website `wavefrontshaping.net <wavefrontshaping.net>`_.

Citing the code
---------------

If you find this tool usefull, please consider citing our paper:

`M. Matthès, P. del Hougne, J. de Rosny, G. Lerosey, and S. Popoff, 
"Optical complex media as universal reconfigurable linear operators," 
Optica *6*, 465-472 (2019). <https://doi.org/10.1364/OPTICA.6.000465>`_

Content
-------

.. toctree::
    :maxdepth: 2
    :caption: Contents
    
    introduction
    installation
    usage







Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
