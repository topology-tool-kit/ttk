/// \ingroup base
/// \class ttk::BoundingVolumeHierarchy
/// \author Rosty Hnatyshyn <rostyslav.hnatyshyn@gmail.com>
/// \date 10.11.2020
///
/// \brief Acceleration structure for native ray tracer.
/// Based on implementation described in Physically Based Rendering:
/// From Theory to Implementation by Matt Pharr, Wenzel Jakob and
/// Greg Humphreys.

#pragma once

#include "Ray.h"
#include <Geometry.h>
#include <algorithm>
#include <limits>
#include <stack>
#include <vector>

namespace ttk {
  template <typename IT>
  class BoundingVolumeHierarchy {
  protected:
    struct Node {
      Node() {
      }
      // leaf constructor
      Node(const std::vector<int> &triangleIndices,
           const size_t nTriangles,
           const float *pMin,
           const float *pMax) {
        m_minX = pMin[0];
        m_minY = pMin[1];
        m_minZ = pMin[2];
        m_maxX = pMax[0];
        m_maxY = pMax[1];
        m_maxZ = pMax[2];
        indices = triangleIndices;
        numTriangles = nTriangles;
        m_left = m_right = nullptr;
      }
      // interior constructor
      Node(const int axis, Node *left, Node *right) {
        numTriangles = 0;
        m_left = left;
        m_right = right;
        m_splitAxis = axis;
        m_minX = std::min(left->m_minX, right->m_minX);
        m_minY = std::min(left->m_minY, right->m_minY);
        m_minZ = std::min(left->m_minZ, right->m_minZ);
        m_maxX = std::max(left->m_maxX, right->m_maxX);
        m_maxY = std::max(left->m_maxY, right->m_maxY);
        m_maxZ = std::max(left->m_maxZ, right->m_maxZ);
      }

      ~Node() {
        if(m_left)
          delete m_left;
        if(m_right)
          delete m_right;
      }

      std::vector<int> indices;
      int numTriangles;
      float m_minX, m_minY, m_minZ;
      float m_maxX, m_maxY, m_maxZ;
      Node *m_left;
      Node *m_right;
      int m_splitAxis;
    };

    struct Triangle {
      int m_index;
      float m_centroid_x, m_centroid_y, m_centroid_z;
      float m_minX, m_minY, m_minZ;
      float m_maxX, m_maxY, m_maxZ;
      Triangle() {
      }

      void init(const int &index,
                const float &centroid_x,
                const float &centroid_y,
                const float &centroid_z,
                const float *pMin,
                const float *pMax) {
        m_index = index;
        m_centroid_x = centroid_x;
        m_centroid_y = centroid_y;
        m_centroid_z = centroid_z;
        m_minX = pMin[0];
        m_minY = pMin[1];
        m_minZ = pMin[2];

        m_maxX = pMax[0];
        m_maxY = pMax[1];
        m_maxZ = pMax[2];
      }
    };

  public:
    BoundingVolumeHierarchy(const float *coords,
                            const IT *connectivityList,
                            const size_t &nTriangles) {
      std::vector<Triangle> triangles;
      buildTriangleList(triangles, coords, connectivityList, nTriangles);

      this->nodes = buildTree(triangles, 0, nTriangles);
    }

    ~BoundingVolumeHierarchy() {
      if(this->nodes)
        delete this->nodes;
    }

    Node *
      buildTree(std::vector<Triangle> &triangles, size_t start, size_t end) {

      float minX, minY, minZ;
      float maxX, maxY, maxZ;
      minX = minY = minZ = std::numeric_limits<float>::max();
      maxX = maxY = maxZ = std::numeric_limits<float>::min();
      for(size_t i = start; i < end; i++) {
        const Triangle &t = triangles[i];
        minX = std::min(t.m_minX, minX);
        minY = std::min(t.m_minY, minY);
        minZ = std::min(t.m_minZ, minZ);

        maxX = std::max(t.m_maxX, maxX);
        maxY = std::max(t.m_maxY, maxY);
        maxZ = std::max(t.m_maxZ, maxZ);
      }
      int numberTriangles = end - start;
      if(numberTriangles == 1) {
        const Triangle &t = triangles[start];
        std::vector<int> indices = {t.m_index};
        float pMin[3] = {t.m_minX, t.m_minY, t.m_minZ};
        float pMax[3] = {t.m_maxX, t.m_maxY, t.m_maxZ};
        return new Node(indices, 1, pMin, pMax);
      } else {
        // find the bounds of the centroids, figure out what dimension to split
        // on
        float cminX, cminY, cminZ;
        float cmaxX, cmaxY, cmaxZ;
        cminX = cminY = cminZ = std::numeric_limits<float>::max();
        cmaxX = cmaxY = cmaxZ = std::numeric_limits<float>::min();
        for(size_t i = start; i < end; i++) {
          const Triangle &t = triangles[i];
          cminX = std::min(t.m_centroid_x, cminX);
          cminY = std::min(t.m_centroid_y, cminY);
          cminZ = std::min(t.m_centroid_z, cminZ);

          cmaxX = std::max(t.m_centroid_x, cmaxX);
          cmaxY = std::max(t.m_centroid_y, cmaxY);
          cmaxZ = std::max(t.m_centroid_z, cmaxZ);
        }
        // figure out the biggest extent, use that dimension
        float diffX = std::abs(cmaxX - cminX);
        float diffY = std::abs(cmaxY - cminY);
        float diffZ = std::abs(cmaxZ - cminZ);
        float maximumExtent = std::max({diffX, diffY, diffZ});
        int axis;
        float minToCheck, maxToCheck;
        if(maximumExtent == diffX) {
          axis = 0;
          minToCheck = cminX;
          maxToCheck = cmaxX;
        } else if(maximumExtent == diffY) {
          axis = 1;
          minToCheck = cminY;
          maxToCheck = cmaxY;
        } else {
          axis = 2;
          minToCheck = cminZ;
          maxToCheck = cmaxZ;
        }
        size_t half = (start + end) / 2;
        // partition triangles into two sets and build children
        if(minToCheck == maxToCheck) {

          std::vector<int> triangleIndices;
          for(size_t i = start; i < end; i++) {
            triangleIndices.push_back(triangles[i].m_index);
          }
          float pMin[3] = {minX, minY, minZ};
          float pMax[3] = {maxX, maxY, maxZ};
          return new Node(triangleIndices, triangleIndices.size(), pMin, pMax);
        } else {
          // partition triangles into equally sized subsets
          std::nth_element(&triangles[start], &triangles[half],
                           &triangles[end - 1] + 1,
                           [axis](const Triangle &t1, const Triangle &t2) {
                             switch(axis) {
                               case 0:
                               default:
                                 return t1.m_centroid_x < t2.m_centroid_x;
                                 break;
                               case 1:
                                 return t1.m_centroid_y < t2.m_centroid_y;
                                 break;
                               case 2:
                                 return t1.m_centroid_z < t2.m_centroid_z;
                                 break;
                             }
                           });

          return new Node(axis, buildTree(triangles, start, half),
                          buildTree(triangles, half, end));
        }
      }

      return nullptr;
    }
    int buildTriangleList(std::vector<Triangle> &triangles,
                          const float *coords,
                          const IT *connectivityList,
                          const size_t &nTriangles) {
      triangles.resize(nTriangles);
      for(size_t ti = 0; ti < nTriangles; ti++) {

        const IT v1 = connectivityList[ti * 3 + 0] * 3;
        const IT v2 = connectivityList[ti * 3 + 1] * 3;
        const IT v3 = connectivityList[ti * 3 + 2] * 3;

        const float &x1 = coords[v1 + 0];
        const float &x2 = coords[v2 + 0];
        const float &x3 = coords[v3 + 0];

        const float &y1 = coords[v1 + 1];
        const float &y2 = coords[v2 + 1];
        const float &y3 = coords[v3 + 1];

        const float &z1 = coords[v1 + 2];
        const float &z2 = coords[v2 + 2];
        const float &z3 = coords[v3 + 2];

        const float pMin[3] = {std::min({x1, x2, x3}), std::min({y1, y2, y3}),
                               std::min({z1, z2, z3})};
        const float pMax[3] = {std::max({x1, x2, x3}), std::max({y1, y2, y3}),
                               std::max({z1, z2, z3})};

        triangles[ti].init(ti, findCentroid(x1, x2, x3),
                           findCentroid(y1, y2, y3), findCentroid(z1, z2, z3),
                           pMin, pMax);
      }

      return 1;
    }

    bool MollerTrumbore(Ray &ray,
                        const IT v0,
                        const IT v1,
                        const IT v2,
                        const float *vertexCoords) const {
      constexpr float kEpsilon = 1e-8;

      float v0v1[3], v0v2[3], pvec[3], tvec[3], qvec[3];
      ttk::Geometry::subtractVectors(
        &vertexCoords[v0], &vertexCoords[v1], v0v1);
      ttk::Geometry::subtractVectors(
        &vertexCoords[v0], &vertexCoords[v2], v0v2);
      ttk::Geometry::crossProduct(ray.m_direction, v0v2, pvec);
      float det = ttk::Geometry::dotProduct(v0v1, pvec);
      if(det > -kEpsilon && det < kEpsilon)
        return false;

      float invDet = 1.0f / det;

      ttk::Geometry::subtractVectors(&vertexCoords[v0], ray.m_origin, tvec);
      float u = ttk::Geometry::dotProduct(tvec, pvec) * invDet;
      if(u < 0.0 || u > 1.0)
        return false;

      ttk::Geometry::crossProduct(tvec, v0v1, qvec);
      float v = ttk::Geometry::dotProduct(ray.m_direction, qvec) * invDet;
      if(v < 0.0 || u + v > 1.0)
        return false;

      float t = ttk::Geometry::dotProduct(v0v2, qvec) * invDet;
      ray.distance = t;
      ray.u = u;
      ray.v = v;

      return true;
    }

    bool intersect(Ray &r,
                   const IT *connectivityList,
                   const float *vertexCoords,
                   int *triangleIndex,
                   float *distance) const {
      bool wasHit = false;
      float nearestTriangle = std::numeric_limits<float>::max();
      std::stack<Node *> stack;
      Node *node = &nodes[0];
      if(!wasNodeHit(r, node)) {
        return false;
      }
      stack.push(node);
      while(stack.size() != 0) {
        node = stack.top();
        stack.pop();
        if(wasNodeHit(r, node)) {

          if(node->numTriangles > 0) {

            for(int i = 0; i < node->numTriangles; i++) {
              bool hasHit = false;
              int triIdx = node->indices[i];

              IT v0 = connectivityList[triIdx * 3 + 0];
              IT v1 = connectivityList[triIdx * 3 + 1];
              IT v2 = connectivityList[triIdx * 3 + 2];
              v0 *= 3;
              v1 *= 3;
              v2 *= 3;
              hasHit = MollerTrumbore(r, v0, v1, v2, vertexCoords);
              if(hasHit && r.distance < nearestTriangle) {
                *triangleIndex = triIdx;
                nearestTriangle = r.distance;
                wasHit = true;
                *distance = r.distance;
              }
            }
          } else {

            if(node->m_right != nullptr) {
              stack.push(node->m_right);
            }
            if(node->m_left != nullptr) {
              stack.push(node->m_left);
            }
          }
        }
      }
      return wasHit;
    }

    bool wasNodeHit(const Ray &r, Node *n) const {
      float tmin = (n->m_minX - r.m_origin[0]) / r.m_direction[0];
      float tmax = (n->m_maxX - r.m_origin[0]) / r.m_direction[0];

      if(tmin > tmax)
        std::swap(tmin, tmax);

      float tymin = (n->m_minY - r.m_origin[1]) / r.m_direction[1];
      float tymax = (n->m_maxY - r.m_origin[1]) / r.m_direction[1];

      if(tymin > tymax)
        std::swap(tymin, tymax);

      if((tmin > tymax) || (tymin > tmax))
        return false;

      if(tymin > tmin)
        tmin = tymin;

      if(tymax < tmax)
        tmax = tymax;

      float tzmin = (n->m_minZ - r.m_origin[2]) / r.m_direction[2];
      float tzmax = (n->m_maxZ - r.m_origin[2]) / r.m_direction[2];

      if(tzmin > tzmax)
        std::swap(tzmin, tzmax);

      if((tmin > tzmax) || (tzmin > tmax))
        return false;

      if(tzmin > tmin)
        tmin = tzmin;

      if(tzmax < tmax)
        tmax = tzmax;

      return true;
    }

  private:
    Node *nodes;
    float findCentroid(const float &v1, const float &v2, const float &v3) {
      return (v1 + v2 + v3) / 3;
    }
  };
} // namespace ttk
