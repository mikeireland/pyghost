************
Introduction
************

.. highlight:: python

PyGHOST is a pure python (i.e. no IRAF) analysis and simulation package for
GHOST data. The hope that is that it will become neat enough for easy 
incoporation into the Gemini recipe system.

 * reduction - the main module containing the reduction code.
 * simulation - the main module containing the simulation code.
 * doc - this documentation.
 * cal - calibration data

Basic Code Use
==============

As of May 2015, this code can only easily be used by running/copying the example lines of
code down the bottom of ghostsim.py.

Simulation Algorithms
=====================

The simulator is not a complete optical simulator - e.g. it doesn't directly deal with 
lenses and only propagates light according to first principles in collimated beams. The
two key equations in use are the grating equation, which in the plane orthogonal to the grating lines is:

.. math::

    n\lambda = d (\sin(\theta_i) - \sin(\theta_o)),

where :math: `n` is the order, :math: `\lambda` the wavelength, :math: `\theta_i` the input angle and :\math `theta_o` the output angle. In three dimensions, this becomes the following vector equation:

.. math::

    \mathbf{\hat{v}} \cdot \mathbf{\hat{s}} = \mathbf{\hat{u}} \cdot \mathbf{\hat{s}} + \frac{n\lambda}{d}

with :math:`\mathbf{s}` a unit vector perpendicular to the grating lines in the plane of the grating. Snell's law is also included in a similar way:

.. math::

    n_o \sin(\theta_o) &= n_i \sin(\theta_i) \\
    \mathbf{\hat{v}} &= \mathbf{\hat{n}} \cos(\theta_o) + \mathbf{\hat{p}} \sin(theta_i),
    
where :math: `\mathbf{\hat{n}}` is the surface normal, and `\mathbf{\hat{p}}` is a unit vector in both in the plane of the surface and in the plane defined by the input vector and surface normal.
