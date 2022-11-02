#include <Triangulation.h>

using namespace std;
using namespace ttk;

Triangulation::Triangulation() : abstractTriangulation_{nullptr} {
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

  switch(rhs.getType()) {
    case Type::EXPLICIT:
      this->abstractTriangulation_ = &this->explicitTriangulation_;
      break;
    case Type::COMPACT:
      this->abstractTriangulation_ = &this->compactTriangulation_;
      break;
    case Type::IMPLICIT:
      this->abstractTriangulation_ = &this->implicitTriangulation_;
      break;
    case Type::HYBRID_IMPLICIT:
      this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation_;
      break;
    case Type::PERIODIC:
      this->abstractTriangulation_ = &this->periodicImplicitTriangulation_;
      break;
    case Type::HYBRID_PERIODIC:
      this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation_;
      break;
  }
}

Triangulation::Triangulation(Triangulation &&rhs) noexcept
  : AbstractTriangulation(
    std::move(*static_cast<AbstractTriangulation *>(&rhs))),
    abstractTriangulation_{nullptr}, explicitTriangulation_{std::move(
                                       rhs.explicitTriangulation_)},
    implicitTriangulation_{std::move(rhs.implicitTriangulation_)},
    periodicImplicitTriangulation_{
      std::move(rhs.periodicImplicitTriangulation_)},
    compactTriangulation_{std::move(rhs.compactTriangulation_)} {

  gridDimensions_ = rhs.gridDimensions_;
  hasPeriodicBoundaries_ = rhs.hasPeriodicBoundaries_;

  switch(rhs.getType()) {
    case Type::EXPLICIT:
      this->abstractTriangulation_ = &this->explicitTriangulation_;
      break;
    case Type::COMPACT:
      this->abstractTriangulation_ = &this->compactTriangulation_;
      break;
    case Type::IMPLICIT:
      this->abstractTriangulation_ = &this->implicitTriangulation_;
      break;
    case Type::HYBRID_IMPLICIT:
      this->abstractTriangulation_ = &this->implicitPreconditionsTriangulation_;
      break;
    case Type::PERIODIC:
      this->abstractTriangulation_ = &this->periodicImplicitTriangulation_;
      break;
    case Type::HYBRID_PERIODIC:
      this->abstractTriangulation_ = &this->periodicPreconditionsTriangulation_;
      break;
  }
  rhs.abstractTriangulation_ = nullptr;
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

    switch(rhs.getType()) {
      case Type::EXPLICIT:
        this->abstractTriangulation_ = &this->explicitTriangulation_;
        break;
      case Type::COMPACT:
        this->abstractTriangulation_ = &this->compactTriangulation_;
        break;
      case Type::IMPLICIT:
        this->abstractTriangulation_ = &this->implicitTriangulation_;
        break;
      case Type::HYBRID_IMPLICIT:
        this->abstractTriangulation_
          = &this->implicitPreconditionsTriangulation_;
        break;
      case Type::PERIODIC:
        this->abstractTriangulation_ = &this->periodicImplicitTriangulation_;
        break;
      case Type::HYBRID_PERIODIC:
        this->abstractTriangulation_
          = &this->periodicPreconditionsTriangulation_;
        break;
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

    switch(rhs.getType()) {
      case Type::EXPLICIT:
        this->abstractTriangulation_ = &this->explicitTriangulation_;
        break;
      case Type::COMPACT:
        this->abstractTriangulation_ = &this->compactTriangulation_;
        break;
      case Type::IMPLICIT:
        this->abstractTriangulation_ = &this->implicitTriangulation_;
        break;
      case Type::HYBRID_IMPLICIT:
        this->abstractTriangulation_
          = &this->implicitPreconditionsTriangulation_;
        break;
      case Type::PERIODIC:
        this->abstractTriangulation_ = &this->periodicImplicitTriangulation_;
        break;
      case Type::HYBRID_PERIODIC:
        this->abstractTriangulation_
          = &this->periodicPreconditionsTriangulation_;
        break;
    }
  }
  AbstractTriangulation::operator=(std::move(rhs));

  return *this;
}

Triangulation::~Triangulation() = default;

void Triangulation::switchGrid(const bool usePeriodic,
                               const bool usePreconditions) {
  if(abstractTriangulation_ != nullptr
     && abstractTriangulation_ != &implicitTriangulation_
     && abstractTriangulation_ != &implicitPreconditionsTriangulation_
     && abstractTriangulation_ != &periodicImplicitTriangulation_
     && abstractTriangulation_ != &periodicPreconditionsTriangulation_) {
    return;
  }

  if(abstractTriangulation_ != nullptr) {
    // clear preconditions when switching triangulation type
    if(abstractTriangulation_ == &implicitPreconditionsTriangulation_
       && (usePeriodic || !usePreconditions)) {
      implicitPreconditionsTriangulation_.clear();
    } else if(abstractTriangulation_ == &periodicPreconditionsTriangulation_
              && (!usePeriodic || !usePreconditions)) {
      periodicPreconditionsTriangulation_.clear();
    }
  }

  if(!usePeriodic && !usePreconditions) {
    abstractTriangulation_ = &implicitTriangulation_;
    implicitTriangulation_.preconditionVerticesAndCells();
  } else if(!usePeriodic && usePreconditions) {
    abstractTriangulation_ = &implicitPreconditionsTriangulation_;
    implicitPreconditionsTriangulation_.preconditionVerticesAndCells();
  } else if(usePeriodic && !usePreconditions) {
    abstractTriangulation_ = &periodicImplicitTriangulation_;
    periodicImplicitTriangulation_.preconditionVerticesAndCells();
  } else if(usePeriodic && usePreconditions) {
    abstractTriangulation_ = &periodicPreconditionsTriangulation_;
    periodicPreconditionsTriangulation_.preconditionVerticesAndCells();
  }
}

bool Triangulation::processImplicitStrategy(const STRATEGY strategy) const {

  if(strategy == STRATEGY::DEFAULT) {

#ifndef TTK_IMPLICIT_PRECONDITIONS_THRESHOLD
#define TTK_IMPLICIT_PRECONDITIONS_THRESHOLD 256 * 256 * 256
#endif // TTK_IMPLICIT_PRECONDITIONS_THRESHOLD

    const uint64_t threshold{TTK_IMPLICIT_PRECONDITIONS_THRESHOLD};
    const uint64_t nVerts
      = gridDimensions_[0] * gridDimensions_[1] * gridDimensions_[2];

    // disable preconditioning above TTK_IMPLICIT_PRECONDITIONS_THRESHOLD
    const auto doPreconditioning = nVerts <= threshold;

    if(!doPreconditioning) {

      // ensure that the warning message is printed at creation time
      int prevDebugLevel{-1};
      if(this->debugLevel_ < 2) {
        prevDebugLevel = this->debugLevel_;
        this->debugLevel_ = 2;
      }
      this->printWrn("Large grid detected (> " + std::to_string(threshold)
                     + " vertices)");
      this->printWrn("Defaulting to the fully implicit triangulation");
      if(prevDebugLevel != -1) {
        this->debugLevel_ = prevDebugLevel;
      }
    }
    return doPreconditioning;
  } else if(strategy == STRATEGY::WITH_PRECONDITIONS) {
    return true;
  }
  return false;
}
