########
tesscube
########

.. image:: https://github.com/tessgi/tesscube/actions/workflows/python-app.yml/badge.svg
    :target: https://github.com/tessgi/tesscube/actions/workflows/python-app.yml
    :alt: Test status

.. image:: https://badge.fury.io/py/tesscube.svg
    :target: https://badge.fury.io/py/tesscube
    :alt: PyPI version

.. image:: https://img.shields.io/badge/documentation-live-blue.svg
    :target: https://tessgi.github.io/tesscube/
    :alt: Documentation

.. <!-- intro content start -->

``tesscube`` is a package designed to help you obtain TESS data by cutting it out of the FFI cubes at the `Barbara A. Mikulski Archive for Space Telescopes (MAST) <https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html>`_.

``tesscube`` works with data that is available in the cloud, and will return TESS data in formats similar to the NASA TESS official mission products. You do not need any login credentials to use `tesscube`, and can use this tool by pip installing either on your local machine or in the cloud.

.. note:: Data formats
  `tesscube` will return `astropy` native objects, including `astropy.fits.HDUList` objects, and `astropy.wcs.WCS` objects. You can work with these directly, or load them into `lightkurve`.

.. <!-- intro content end -->

.. <!-- quickstart content start -->

Quickstart
==========

The easiest way to install ``tesscube`` and all of its dependencies is to use the ``pip`` command.

To install ``tesscube``, run the following command in a terminal window:

.. code-block:: console

  $ python -m pip install tesscube --upgrade

The ``--upgrade`` flag is optional, but recommended if you already
have ``tesscube`` installed and want to upgrade to the latest version.

Depending on the specific Python environment, you may need to replace ``python``
with the correct Python interpreter, e.g., ``python3``.

Load an FFI cube
----------------

You can work with an FFI cube by loading it using a sector, camera, and CCD number.

.. code-block:: python

  from tesscube import TESSCube
  cube = TESSCube(sector=1, camera=1, ccd=4)


Obtain an FFI Using `tesscube`
------------------------------

You can obtain an FFI image by indexing into a cube

.. code-block:: python

  from tesscube import TESSCube
  cube = TESSCube(sector=1, camera=1, ccd=4)
  ffi = cube[300]

This will return an `astropy.fits.HDUList`


Obtain a TPF
------------

You can obtain a TPF in two ways, either you can either pass a pixel position

.. code-block:: python

  from tesscube import TESSCube
  from astropy.coordinates import SkyCoord
  corner = (1282, 1750)
  cube = TESSCube(sector=1, camera=1, ccd=4)
  tpf = cube.get_tpf(corner, shape=(10, 11))

Or you can pass an astropy SkyCoord object containing the RA and Dec of the target

.. code-block:: python

  from tesscube import TESSCube
  from astropy.coordinates import SkyCoord
  coord = SkyCoord.from_name("AU Mic")
  cube = TESSCube(sector=1, camera=1, ccd=4)
  tpf = cube.get_tpf(coord, shape=(10, 11))

Alternatively, you can index into the cube like so:

.. code-block:: python

  from tesscube import TESSCube
  cube = TESSCube(sector=1, camera=1, ccd=4)
  tpf = cube[:, 401:410, 503:510]

Both will return an `astropy.fits.HDUList`, with a file format similar to the official mission products.

Obtain a time-coadded TPF
-------------------------

TESS data can be coadded in time to increase signal to noise at the expense of time resolution. You can obtain a lower time resolution by either passing in a `frame_bin` parameter, which will downsample the resultant TPF,

.. code-block:: python

  from tesscube import TESSCube
  from astropy.coordinates import SkyCoord
  corner = (1282, 1750)
  cube = TESSCube(sector=1, camera=1, ccd=4)
  tpf = cube.get_tpf(corner, shape=(10, 11), frame_bin=10)

Or you can slice the cube, which will return a downsampled TPF


.. code-block:: python
  
  from tesscube import TESSCube
  cube = TESSCube(sector=1, camera=1, ccd=4)
  tpf = cube[::10, 401:410, 503:510]

Both will return an `astropy.fits.HDUList`, with a file format similar to the official mission products, with the time resolution reduced by a factor of 10.

Obtain a TICA TPF
-----------------

TICA data are FFIs produced by MIT with a quick turn around time that are rapidly delivered to MIT to enable quick processing of data for transient events. You can create a TPF out of TICA products, but this functionality is potentially brittle.

TICA FFIs do not always conform to the same standard and sometimes have larger headers. This means they are hard to read in by tesscube, which assumes that there are always the same numbers of bytes. You will have to specify for a given sector how many header blocks each TICA cube has. Usually this is 6, but sometimes it is 7, depending on the sector.

.. code-block:: python

  from tesscube import TESSCubes
  cube = TESSCube(sector=1, camera=1, ccd=4, tica=True, nhdr_blocks=6)

Once this object is initialized you should be able to work with it to extract a TPF

.. code-block:: python

  cube.get_tpf()

TICA FFIs do not have errors. When functions should return errors, these will be filled with zeros.

Currently there is not a way to index into the cube object and get the TICA FFI, because these are not stored in the AWS cloud. 

.. <!-- quickstart content end -->

.. <!-- Contributing content start -->

Contributing
============

``tesscube``  is an open-source, community driven package. 
We welcome users to contribute and develop new features for lksearch.  

For further information, please see the `Lightkurve Community guidelines <https://docs.lightkurve.org/development/contributing.html>`_.

.. <!-- Contributing content end -->

.. <!-- Citing content start -->

Citing
======

If you find ``tesscube`` useful in your research, please cite it and give us a GitHub star!

If you use Lightkurve for work or research presented in a publication, we request the following acknowledgment or citation:

`This research made use of Lightkurve, a Python package for Kepler and TESS data analysis (Lightkurve Collaboration, 2018).`

See full citation instuctions, including dependencies, in the `Lightkurve documentation <https://docs.lightkurve.org/about/citing.html>`_. 

.. <!-- Citing content end -->

.. <!-- Contact content start -->

Contact
=======
``tesscube`` is an open source community project created by the `TESS Science Support Center`_. 
The best way to contact us is to `open an issue`_ or to e-mail tesshelp@bigbang.gsfc.nasa.gov.
  
Please include a self-contained example that fully demonstrates your problem or question.
  .. _`TESS Science Support Center`: https://heasarc.gsfc.nasa.gov/docs/tess/
  .. _`open an issue`: https://github.com/tessgi/tesscube/issues/new


.. <!-- Contact content end -->

.. <!-- Changelog content start -->

Changelog:
==========

  - Fixed barycentric timing correction with lkspacecraft package
  - Added from_name method
  - Added ability to use "TICA" FFIs. This is experimental and might be buggy.
  - Patch removes the un-needed `fitsio` dependency
  - Initial v1.0.0 release of `tesscube`.
  

.. <!-- Changelog content end -->
