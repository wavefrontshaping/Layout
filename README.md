# SLMlayout
A library to create layouts of various shapes for wavefrontshaping using DMDs or SLMs.
It is specifically designed to work with the [ALP4lib](https://github.com/wavefrontshaping/ALP4lib) module to control Vialux DMDs but can also be used for other types of SLMs or DMDs.
In particular, it allows generating **bitplanes** in C format data that can be sent to the DMD with fast transfer time (see Vialux documentation). The corresponding function are coded in Cython for optimized computation speed.

## Citing the code

If you find this tool usefull, please consider citing our paper:
[M. Matthès, P. del Hougne, J. de Rosny, G. Lerosey, and S. Popoff, "Optical complex media as universal reconfigurable linear operators," Optica *6*, 465-472 (2019).](https://doi.org/10.1364/OPTICA.6.000465)


## Installation

Simply install it with

```shell
pip install SLMlayout 
```

Alternatively, you can download the files and execute the following command

```shell
python setup.py install
```

## Presentation and example

Please read our short presentation of the module on [https://wavefrontshaping.net]() and see our Jupyter notebook example codes [here](https://github.com/wavefrontshaping/WFS.net/tree/master/Layout).
