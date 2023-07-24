# Gordon-Wixom surface

A *C1* transfinite interpolation surface based on a
[classic paper](https://doi.org/10.1137/0711072).
Uses my [geometry library](https://github.com/salvipeter/libgeom/).

In this variation the function to be interpolated maps into 3D,
so it can be used for hole filling.
The domain is assumed to be a regular polygon;
boundary constraints are given as pairs of B-spline curves.

The test program reads a file of the following format:
```
<number of sides>
<1st outer B-spline>
<1st inner B-spline>
<2nd outer B-spline>
<2nd inner B-spline>
...
```
where a B-spline is described as:
```
<degree>
<number of knots> <knot values: u1 u2 ...>
<number of control points>
<control points 1: x1 y1 z1>
<control points 2: x2 y2 z2>
```
