#pragma once

#include <geometry.hh>

// Based on
//   W.J. Gordon, J.A. Wixom: Pseudo-harmonic interpolation on convex domains
//   SIAM Journal on Numerical Analysis 11(5):909-933, 1974.
//   https://doi.org/10.1137/0711072

// But:
// - the domain is always a regular n-sided polygon
// - the boundary constraints are given as pairs of B-spline curves
//   (normal derivatives are evaluated as `inner(t) - outer(t)`)

// The boundary curves are assumed to be
// - parameterized in [0,1]
// - correctly oriented (curve at 1 = next curve at 0)

// As in the original paper, the boundary constraint may be discontinuous.

class GordonWixom {
public:
  using CurveVector = std::vector<std::shared_ptr<Geometry::BSCurve>>;
  GordonWixom(const CurveVector &outer, const CurveVector &inner);
  Geometry::Point2DVector domain() const;
  Geometry::Point3D eval(const Geometry::Point2D &uv) const;
  Geometry::TriMesh eval(size_t resolution) const;
private:
  Geometry::Point3D hermitePoint(size_t i1, double t1, double cos1,
                                 size_t i2, double t2, double cos2,
                                 double ratio) const;
  Geometry::Point2DVector parameters(size_t resolution) const;
  Geometry::TriMesh meshTopology(size_t resolution) const;

  size_t n_;
  CurveVector outer_, inner_;
  Geometry::Point2DVector domain_;
};
