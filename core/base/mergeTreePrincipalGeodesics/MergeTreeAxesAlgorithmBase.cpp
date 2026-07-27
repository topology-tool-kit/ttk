#include <MergeTreeAxesAlgorithmBase.h>

//----------------------------------------------------------------------------
// Setter
//----------------------------------------------------------------------------
void ttk::MergeTreeAxesAlgorithmBase::setDeterministic(
  const bool deterministic) {
  deterministic_ = deterministic;
}

void ttk::MergeTreeAxesAlgorithmBase::setNumberOfProjectionSteps(
  const unsigned int k) {
  k_ = k;
}

void ttk::MergeTreeAxesAlgorithmBase::setBarycenterSizeLimitPercent(
  const double barycenterSizeLimitPercent) {
  barycenterSizeLimitPercent_ = barycenterSizeLimitPercent;
}

void ttk::MergeTreeAxesAlgorithmBase::setProbabilisticVectorsInit(
  const bool probabilisticVectorsInit) {
  probabilisticVectorsInit_ = probabilisticVectorsInit;
}
