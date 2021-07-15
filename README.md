# hyperbolic

Hyperbolic geometry and rigidbody physics library.

This is a header library of functions related to the geometry and physics of hyperbolic space using the hyperboloid model.  The only pre-requisite is the glm vector library.

* glm_lorentz.h - consists basic functions related to the Lorentz group of orientation and time preserving symmetries of the 3+1 dimensional Riemannian space.  It is written in the style and builds on the glm library.

* vector_types.h - the choice of vector types used in the following files, see this for changing between float/double precision.

* imaginary_extension.h - adds an imaginary component to a type, used internally by glm_lorentz.h

* hyperbolic.h - basic functions for reasoning with and manipulating points, lines and planes in hyperbolic space.

* lorentz_lie_algebra.h - functions related to the lie algebra of the Lorentz group, elements of which represent the velocity of a frame in hyperbolic space.

* rigid_body.h - a representation of a rigid body in hyperbolic space, along with functions for simulating its mass distribution and physics.

## Things you should know

The hyperboloid representation of hyperbolic space is convenient and efficient, however numerical issues occur quickly as one moves away from the center (0, 0, 0, 1).
So this representation is only appropriate for representing points / frames relative to another point / frame of interest.  For example floating point precision is sufficient for rendering objects visible from the frame.

## Potential future features

* a point class for representing point in hyperbolic space far from the origin