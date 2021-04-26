Typical usage
=============


Amplitude and phase modulation with DMDs
----------------------------------------

By nature, Digital Micromirror Devices (**DMDs**) and binary amplitude
modulators.

They are sometimes preferred to liquid crystal phase modulators for their speed, 
few kHz compared to 10 to 100Hz typically.

However, phase modulation is often required for wavefront shaping manipulation 
as it is the best way to control interference.
To do so, one popular technique is the `Lee Hologram approach <https://www.wavefrontshaping.net/post/id/16>`_.
It consists in:
* introducing a low pass filter in the Fourier plane of the DMD,
* encoding the phase of the field in the spatial phase of amplitude fringes.


.. figure::  /_static/images/setup.png
   :align:   center

   Typical setup for amplitude and phase modulation using DMDs.


.. code-block:: python3

    import SLMlayout   
    layout = Layout.Hexagons(radius = 350, 
                         hexSize = 60,
                         resolution = [1200,1920], 
                         center = [1200//2,1920//2], 
                         gap = 3)





.. automodule:: SLMlayout.core
    :members:

.. autoclass:: SLMlayout.core.Layout
    :members:

.. automethod:: SLMlayout.core.Layout.getMaskFromImage

.. automodule:: SLMlayout.hexagons
    :members:
    :undoc-members:
    :show-inheritance:


    
.. toctree::
   :maxdepth: 2
   :caption: Contents: