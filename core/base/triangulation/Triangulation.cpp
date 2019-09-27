#include <Triangulation.h>

using namespace std;
using namespace ttk;

Triangulation::Triangulation()
  : AbstractTriangulation{}, gridDimensions_{-1, -1, -1},
    abstractTriangulation_{nullptr}, usePeriodicBoundaries_{false} {
  debugLevel_ = 0; // overrides the global debug level.
}

Triangulation::Triangulation(const Triangulation &rhs)
  : AbstractTriangulation(rhs), gridDimensions_{rhs.gridDimensions_},
    abstractTriangulation_{nullptr},
    explicitTriangulation_{rhs.explicitTriangulation_},
    implicitTriangulation_{rhs.implicitTriangulation_},
    periodicImplicitTriangulation_{rhs.periodicImplicitTriangulation_},
    usePeriodicBoundaries_{rhs.usePeriodicBoundaries_} {

  if(rhs.abstractTriangulation_ == &rhs.explicitTriangulation_) {
    abstractTriangulation_ = &explicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.implicitTriangulation_) {
    abstractTriangulation_ = &implicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation_) {
    abstractTriangulation_ = &periodicImplicitTriangulation_;
  }
}

Triangulation::Triangulation(Triangulation &&rhs)
  : AbstractTriangulation(std::move(rhs)), gridDimensions_{std::move(
                                             rhs.gridDimensions_)},
    abstractTriangulation_{nullptr}, explicitTriangulation_{std::move(
                                       rhs.explicitTriangulation_)},
    implicitTriangulation_{std::move(rhs.implicitTriangulation_)},
    periodicImplicitTriangulation_{
      std::move(rhs.periodicImplicitTriangulation_)},
    usePeriodicBoundaries_{std::move(rhs.usePeriodicBoundaries_)} {

  if(rhs.abstractTriangulation_ == &rhs.explicitTriangulation_) {
    abstractTriangulation_ = &explicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.implicitTriangulation_) {
    abstractTriangulation_ = &implicitTriangulation_;
  } else if(rhs.abstractTriangulation_ == &rhs.periodicImplicitTriangulation_) {
    abstractTriangulation_ = &periodicImplicitTriangulation_;
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
    usePeriodicBoundaries_ = rhs.usePeriodicBoundaries_;

    if(rhs.abstractTriangulation_ == &rhs.explicitTriangulation_) {
      abstractTriangulation_ = &explicitTriangulation_;
    } else if(rhs.abstractTriangulation_ == &rhs.implicitTriangulation_) {
      abstractTriangulation_ = &implicitTriangulation_;
    } else if(rhs.abstractTriangulation_
              == &rhs.periodicImplicitTriangulation_) {
      abstractTriangulation_ = &periodicImplicitTriangulation_;
    }
  }
  return *this;
}

Triangulation &Triangulation::operator=(Triangulation &&rhs) {
  if(this != &rhs) {
    AbstractTriangulation::operator=(std::move(rhs));
    gridDimensions_ = std::move(rhs.gridDimensions_);
    abstractTriangulation_ = nullptr;
    explicitTriangulation_ = std::move(rhs.explicitTriangulation_);
    implicitTriangulation_ = std::move(rhs.implicitTriangulation_);
    periodicImplicitTriangulation_
      = std::move(rhs.periodicImplicitTriangulation_);
    usePeriodicBoundaries_ = std::move(rhs.usePeriodicBoundaries_);

    if(rhs.abstractTriangulation_ == &rhs.explicitTriangulation_) {
      abstractTriangulation_ = &explicitTriangulation_;
    } else if(rhs.abstractTriangulation_ == &rhs.implicitTriangulation_) {
      abstractTriangulation_ = &implicitTriangulation_;
    } else if(rhs.abstractTriangulation_
              == &rhs.periodicImplicitTriangulation_) {
      abstractTriangulation_ = &periodicImplicitTriangulation_;
    }
  }
  return *this;
}

Triangulation::~Triangulation() = default;
