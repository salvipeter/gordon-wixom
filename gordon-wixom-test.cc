#include <fstream>

#include "gordon-wixom.hh"

using namespace Geometry;

Point3D readPoint(std::istream &is) {
  Point3D p;
  is >> p[0] >> p[1] >> p[2];
  return p;
}

BSCurve readBSpline(std::istream &is) {
  size_t degree, n_knots, n_cpts;
  DoubleVector knots;
  PointVector cpts;
  is >> degree;
  is >> n_knots;
  knots.resize(n_knots);
  for (size_t i = 0; i < n_knots; ++i)
    is >> knots[i];
  is >> n_cpts;
  cpts.resize(n_cpts);
  for (size_t i = 0; i < n_cpts; ++i)
    cpts[i] = readPoint(is);
  return { degree, knots, cpts };
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " <model.gwx> [resolution]" << std::endl;
    return 1;
  }

  size_t resolution = 15;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  std::ifstream f(argv[1]);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t n;
  f >> n;

  GordonWixom::CurveVector outer, inner;
  for (size_t i = 0; i < n; ++i) {
    outer.push_back(std::make_shared<BSCurve>(readBSpline(f)));
    inner.push_back(std::make_shared<BSCurve>(readBSpline(f)));
  }

  f.close();

  GordonWixom surface(outer, inner);
  surface.eval(resolution).writeOBJ("test.obj");
}
