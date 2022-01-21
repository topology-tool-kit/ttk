#include <Triangulation.h>

using namespace std;
using namespace ttk;

Triangulation::Triangulation()
  : AbstractTriangulation{}, abstractTriangulation_{nullptr} {
  debugLevel_ = 0; // overrides the global debug level.
  gridDimensions_ = {-1, -1, -1};
  hasPeriodicBoundaries_ = false;
}

Triangulation::Triangulation(const Triangulation &rhs)
  : AbstractTriangulation(rhs), abstractTriangulation_{nullptr},
    explicitTriangulation_{rhs.explicitTriangulation_},
    implicitTriangulation_{rhs.implicitTriangulation_},
    periodicImplicitTriangulation_{rhs.periodicImplicitTriangulation_},
    compactTriangulation_{rhs.compactTriangulation_} {

  gridDimensions_ = rhs.gridDimensions_;
  hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

  if(rhs.abstractTriangulation_ == &rhs.explicitTriangulation_) {
    abstractTriangulation_ = &explicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.implicitTriangulation_) {
    abstractTriangulation_ = &implicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation_) {
    abstractTriangulation_ = &periodicImplicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.compactTriangulation_) {
    abstractTriangulation_ = &compactTriangulation_;
  }
}

Triangulation::Triangulation(Triangulation &&rhs) noexcept
  : AbstractTriangulation(std::move(rhs)), abstractTriangulation_{nullptr},
    explicitTriangulation_{std::move(rhs.explicitTriangulation_)},
    implicitTriangulation_{std::move(rhs.implicitTriangulation_)},
    periodicImplicitTriangulation_{
      std::move(rhs.periodicImplicitTriangulation_)},
    compactTriangulation_{std::move(rhs.compactTriangulation_)} {

  gridDimensions_ = rhs.gridDimensions_;
  hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

  if(rhs.abstractTriangulation_ == &rhs.explicitTriangulation_) {
    abstractTriangulation_ = &explicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.implicitTriangulation_) {
    abstractTriangulation_ = &implicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation_) {
    abstractTriangulation_ = &periodicImplicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.compactTriangulation_) {
    abstractTriangulation_ = &compactTriangulation_;
  }
}

Triangulation &Triangulation::operator=(const Triangulation &rhs) {
  if(this != &rhs) {
    AbstractTriangulation::operator=(rhs);
    gridDimensions_ = rhs.gridDimensions_;
    abstractTriangulation_ = nullptr;
    explicitTriangulation_ = rhs.explicitTriangulation_;
    implicitTriangulation_ = rhs.implicitTriangulation_;
    periodicImplicitTriangulation_ = rhs.periodicImplicitTriangulation_;
    compactTriangulation_ = rhs.compactTriangulation_;
    hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

    if(rhs.abstractTriangulation_ == &rhs.explicitTriangulation_) {
      abstractTriangulation_ = &explicitTriangulation_;
    } else if(rhs.abstractTriangulation_ == &rhs.implicitTriangulation_) {
      abstractTriangulation_ = &implicitTriangulation_;
    } else if(rhs.abstractTriangulation_
              == &rhs.periodicImplicitTriangulation_) {
      abstractTriangulation_ = &periodicImplicitTriangulation_;
    } else if(rhs.abstractTriangulation_ == &rhs.compactTriangulation_) {
      abstractTriangulation_ = &compactTriangulation_;
    }
  }
  return *this;
}

Triangulation &Triangulation::operator=(Triangulation &&rhs) noexcept {
  if(this != &rhs) {
    gridDimensions_ = rhs.gridDimensions_;
    abstractTriangulation_ = nullptr;
    explicitTriangulation_ = std::move(rhs.explicitTriangulation_);
    implicitTriangulation_ = std::move(rhs.implicitTriangulation_);
    periodicImplicitTriangulation_
      = std::move(rhs.periodicImplicitTriangulation_);
    compactTriangulation_ = std::move(rhs.compactTriangulation_);
    hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

    if(rhs.abstractTriangulation_ == &rhs.explicitTriangulation_) {
      abstractTriangulation_ = &explicitTriangulation_;
    } else if(rhs.abstractTriangulation_ == &rhs.implicitTriangulation_) {
      abstractTriangulation_ = &implicitTriangulation_;
    } else if(rhs.abstractTriangulation_
              == &rhs.periodicImplicitTriangulation_) {
      abstractTriangulation_ = &periodicImplicitTriangulation_;
    } else if(rhs.abstractTriangulation_ == &rhs.compactTriangulation_) {
      abstractTriangulation_ = &compactTriangulation_;
    }
  }
  AbstractTriangulation::operator=(std::move(rhs));

  return *this;
}

Triangulation::~Triangulation() = default;
