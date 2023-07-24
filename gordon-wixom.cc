#include <cmath>
#include <functional>
#include <optional>

#include "gordon-wixom.hh"

using namespace Geometry;

GordonWixom::GordonWixom(const CurveVector &outer, const CurveVector &inner)
  : n_(outer.size()), outer_(outer), inner_(inner)
{
  // Setup Domain
  if (n_ == 4)
    domain_ = { {1,1}, {-1,1}, {-1,-1}, {1,-1} };
  else {
    double alpha = 2.0 * M_PI / n_;
    for (size_t i = 0; i < n_; ++i)
      domain_.emplace_back(std::cos(alpha * i), std::sin(alpha * i));
  }
}


// Numerical integration

static Point3D trapezoid(const std::function<Point3D(double)> &f, double a, double b,
                         size_t n, const Point3D &s = { 0.0, 0.0, 0.0 }) {
  if (n == 1)
    return (f(a) + f(b)) / 2 * (b - a);
  size_t k = 1 << (n - 2);
  double h = (b - a) / k;
  double start = a + h / 2;
  Point3D sum(0.0, 0.0, 0.0);
  for (size_t i = 0; i < k; ++i)
    sum += f(start + h * i);
  return (s + sum * h) / 2;
}

[[maybe_unused]]
static Point3D integral(const std::function<Point3D(double)> &f, double a, double b,
                        size_t max_iterations, double tolerance) {
  Point3D s, s_prev, st_prev(0.0, 0.0, 0.0);
  for (size_t i = 1; i <= max_iterations; ++i) {
    Point3D st = trapezoid(f, a, b, i, st_prev);
    s = (st * 4 - st_prev) / 3;
    if (i > 5 && ((s - s_prev).norm() < tolerance * s_prev.norm() ||
                  (s.norm() == 0.0 && s_prev.norm() == 0.0)))
      break;
    s_prev = s;
    st_prev = st;
  }
  return s;
}

[[maybe_unused]]
static Point3D quadrature(const std::function<Point3D(double)> &f,
                          double a, double b, size_t n) {
  const static double gauss[] = {-0.861136312, 0.347854845,
                                 -0.339981044, 0.652145155,
                                 0.339981044, 0.652145155,
                                 0.861136312, 0.347854845};
  Point3D sum(0.0, 0.0, 0.0);
  for (size_t j = 0; j < n; ++j) {
    double a1 = a + (b - a) * j / n;
    double b1 = a1 + (b - a) / n;
    for (size_t i = 0; i < 8; i += 2) {
      double u = ((b1 - a1) * gauss[i] + a1 + b1) * 0.5;
      sum += f(u) * gauss[i+1] * (b1 - a1) * 0.5;
    }
  }
  return sum;
}


// Domain & mesh generation

Point2DVector GordonWixom::domain() const {
  return domain_;
}

static size_t meshSize(size_t n, size_t resolution) {
  if (n == 3)
    return (resolution + 1) * (resolution + 2) / 2;
  if (n == 4)
    return (resolution + 1) * (resolution + 1);
  return 1 + n * resolution * (resolution + 1) / 2;
}

Point2DVector GordonWixom::parameters(size_t resolution) const {
  size_t size = meshSize(n_, resolution);
  Point2DVector result;
  result.reserve(size);

  if (n_ == 3) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = domain_[0] * u + domain_[2] * (1 - u);
      auto q = domain_[1] * u + domain_[2] * (1 - u);
      for (size_t k = 0; k <= j; ++k) {
        double v = j == 0 ? 1.0 : (double)k / j;
        result.push_back(p * (1 - v) + q * v);
      }
    }
  } else if (n_ == 4) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      auto p = domain_[0] * (1 - u) + domain_[1] * u;
      auto q = domain_[3] * (1 - u) + domain_[2] * u;
      for (size_t k = 0; k <= resolution; ++k) {
        double v = (double)k / resolution;
        result.push_back(p * (1 - v) + q * v);
      }
    }
  } else { // n_ > 4
    Point2D center(0.0, 0.0);
    result.push_back(center);
    for (size_t j = 1; j <= resolution; ++j) {
      double u = (double)j / (double)resolution;
      for (size_t k = 0; k < n_; ++k)
        for (size_t i = 0; i < j; ++i) {
          double v = (double)i / (double)j;
          Point2D ep = domain_[(k+n_-1)%n_] * (1.0 - v) + domain_[k] * v;
          Point2D p = center * (1.0 - u) + ep * u;
          result.push_back(p);
        }
    }
  }
  return result;
}

TriMesh GordonWixom::meshTopology(size_t resolution) const {
  TriMesh mesh;
  mesh.resizePoints(meshSize(n_, resolution));

  if (n_ == 3) {
    size_t prev = 0, current = 1;
    for (size_t i = 0; i < resolution; ++i) {
      for (size_t j = 0; j < i; ++j) {
        mesh.addTriangle(current + j, current + j + 1, prev + j);
        mesh.addTriangle(current + j + 1, prev + j + 1, prev + j);
      }
      mesh.addTriangle(current + i, current + i + 1, prev + i);
      prev = current;
      current += i + 2;
    }
  } else if (n_ == 4) {
    for (size_t i = 0; i < resolution; ++i)
      for (size_t j = 0; j < resolution; ++j) {
        size_t index = i * (resolution + 1) + j;
        mesh.addTriangle(index, index + resolution + 1, index + 1);
        mesh.addTriangle(index + 1, index + resolution + 1, index + resolution + 2);
      }
  } else { // n_ > 4
    size_t inner_start = 0, outer_vert = 1;
    for (size_t layer = 1; layer <= resolution; ++layer) {
      size_t inner_vert = inner_start, outer_start = outer_vert;
      for (size_t side = 0; side < n_; ++side) {
        size_t vert = 0;
        while(true) {
          size_t next_vert = (side == n_ - 1 && vert == layer - 1) ? outer_start : (outer_vert + 1);
          mesh.addTriangle(inner_vert, outer_vert, next_vert);
          ++outer_vert;
          if (++vert == layer)
            break;
          size_t inner_next = (side == n_ - 1 && vert == layer - 1) ? inner_start : (inner_vert + 1);
          mesh.addTriangle(inner_vert, next_vert, inner_next);
          inner_vert = inner_next;
        }
      }
      inner_start = outer_start;
    }
  }
  return mesh;
}


// Evaluation

static std::optional<double> intersect(const Point2D &p, double theta,
                                       const Point2D &a, const Point2D &b) {
  size_t x = 0, y = 1;
  Vector2D v(std::cos(theta), std::sin(theta));
  if (std::abs(v[0]) < std::abs(v[1]))
    std::swap(x, y);
  double r = v[y] / v[x];
  double denom = (b[y] - a[y]) - r * (b[x] - a[x]);
  if (denom != 0.0) {
    double s = ((p[y] - a[y]) - r * (p[x] - a[x])) / denom;
    if (s >= 0 && s <= 1)
      return { s };
  }
  return {};
}

Point3D GordonWixom::hermitePoint(size_t i1, double t1, double cos1,
                                  size_t i2, double t2, double cos2,
                                  double ratio) const {
  VectorVector der;
  auto p1 = outer_[i1]->eval(t1, 1, der);
  auto du1 = der[1];
  auto dv1 = inner_[i1]->eval(t1) - p1;
  auto d1 = du1 * cos1 + dv1 * (1.0 - std::abs(cos1));
  auto p2 = outer_[i2]->eval(t2, 1, der);
  auto du2 = der[1];
  auto dv2 = inner_[i2]->eval(t2) - p2;
  auto d2 = du2 * cos2 + dv2 * (1.0 - std::abs(cos2));
  BSCurve hermite({ p1, p1 + d1 / 3, p2 + d2 / 3, p2 });
  return hermite.eval(ratio);
}

static double cosangle(const Point2D &a, const Point2D &b, const Point2D &c) {
  auto safenorm = [](const Point2D &u) {
    double n = u.norm();
    return n > 0 ? u / n : u;
  };
  return safenorm(a - b) * safenorm(c - b);
}

// Actually checks only that it is on the line,
// but since we are on a convex domain, it is good enough.
static bool onEdge(const Point2D &p, const Point2D &a, const Point2D &b) {
  auto v = (b - a).normalize();
  return ((a - p) - v * ((a - p) * v)).norm() < epsilon;
}

Point3D GordonWixom::eval(const Point2D &uv) const {
  for (size_t i = 0; i < n_; ++i)
    if (onEdge(uv, domain_[(i+n_-1)%n_], domain_[i]))
      return outer_[i]->eval((uv - domain_[(i+n_-1)%n_]).norm() /
                             (domain_[i] - domain_[(i+n_-1)%n_]).norm());
  auto f = [&](double theta) {
    std::vector<size_t> indices;
    std::vector<Point2D> endpoints;
    std::vector<double> params;
    for (size_t i = 0; i < n_; ++i) {
      auto s = intersect(uv, theta, domain_[(i+n_-1)%n_], domain_[i]);
      if (!s)
        continue;
      indices.push_back(i);
      endpoints.push_back(domain_[(i+n_-1)%n_] * (1 - *s) + domain_[i] * *s);
      params.push_back(*s);
    }
    if (indices.size() < 2)
      throw std::runtime_error("intersection error");
    if (indices.size() > 2) {
      // kutykurutty
      if (indices[1] == indices[0] + 1) {
        indices.erase(indices.begin());
        params.erase(params.begin());
        endpoints.erase(endpoints.begin());
      }
    }
    return hermitePoint(indices[0], params[0], cosangle(uv, endpoints[0], domain_[indices[0]]),
                        indices[1], params[1], cosangle(uv, endpoints[1], domain_[indices[1]]),
                        (uv - endpoints[0]).norm() / (endpoints[1] - endpoints[0]).norm());
  };
  // return integral(f, 0.0, 2.0 * M_PI, 10, 1e-8) / (2.0 * M_PI);
  return quadrature(f, 0.0, 2.0 * M_PI, 100) / (2.0 * M_PI);
}

TriMesh GordonWixom::eval(size_t resolution) const {
  TriMesh mesh = meshTopology(resolution);
  Point2DVector uvs = parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [&](const Point2D &uv) { return eval(uv); });
  mesh.setPoints(points);
  return mesh;
}
