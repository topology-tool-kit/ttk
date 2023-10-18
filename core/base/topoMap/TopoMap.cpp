#include <TopoMap.h>

ttk::TopoMap::TopoMap() {
  // inherited from Debug: prefix will be printed at the beginning of every msg
  this->setDebugMsgPrefix("TopoMap");
}
ttk::TopoMap::~TopoMap() = default;

#if defined(TTK_ENABLE_QHULL) && defined(Qhull_FOUND)

#include <libqhullcpp/Qhull.h>
bool computeConvexHull_aux(const std::vector<double> &coords,
                           std::vector<size_t> &res,
                           std::string &errMsg) {
  size_t nbPoint = coords.size() / 2;
  char qHullFooStr[1] = ""; // We use no Qhull options.
  orgQhull::Qhull qhull;
  try {
    qhull.runQhull(qHullFooStr, 2, nbPoint, coords.data(), qHullFooStr);
  } catch(orgQhull::QhullError &e) {
    errMsg = "Error with qHull module: " + std::string(e.what());
    return false;
  }

  // Qhull gives us the coordinates of the points in the convex hull. Here we
  // retrive the indices of this points in the list we provided. We will also
  // compute the barycenter of the points in the convex hull.
  for(const auto &u : qhull.vertexList()) {
    const orgQhull::QhullPoint &qhullPt = u.point();
    auto coordsCur = qhullPt.coordinates();
    for(size_t j = 0; j < coords.size() / 2; j++) {
      if(fabs(coords[2 * j] - coordsCur[0])
           + fabs(coords[2 * j + 1] - coordsCur[1])
         < EpsilonDBL) {
        res.push_back(j);
        break;
      }
    }
  }

  if(res.size() != qhull.vertexList().size()) {
    errMsg = "Error : could not retrieve all vertices in the convex hull.";
    return false;
  }
  return true;
}

#else // Using boost.

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/geometries/polygon.hpp>
namespace bg = boost::geometry;
BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

using Point = boost::tuple<double, double>;
using Polygon = boost::geometry::model::polygon<Point>;
using Mpoints = boost::geometry::model::multi_point<Point>;

bool computeConvexHull_aux(const std::vector<double> &coords,
                           std::vector<size_t> &res,
                           std::string &errMsg) {
  Mpoints multiPoints;
  Polygon hull;
  size_t nbPoint = coords.size() / 2;
  for(size_t i = 0; i < nbPoint; i++) {
    boost::geometry::append(
      multiPoints, Point(coords[2 * i], coords[2 * i + 1]));
  }

  boost::geometry::convex_hull(multiPoints, hull);

  // From the coordinates of the points in the convex hull, we find their
  // indices.
  for(const auto &boostPt : hull.outer()) {
    double coordsCur[2] = {boostPt.get<0>(), boostPt.get<1>()};
    for(size_t j = 0; j < coords.size() / 2; j++) {
      if(fabs(coords[2 * j] - coordsCur[0])
           + fabs(coords[2 * j + 1] - coordsCur[1])
         < Epsilon) {
        res.push_back(j);
        break;
      }
    }
  }
  if(res.size() != hull.outer().size()) {
    errMsg = "Error : could not retrieve all vertices in the convex hull.";
    return false;
  }

  // Boost closes the polygon, hence the first and the last vertices are the
  // identical. We do not want a duplicate vertex.
  res.pop_back();
  return true;
}
#endif
