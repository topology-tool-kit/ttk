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
