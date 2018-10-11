#ifndef _BOTTLENECKDISTANCEIMPLMAIN_H
#define _BOTTLENECKDISTANCEIMPLMAIN_H

//  vector <   -- diagram
//    tuple <    -- pair of critical points
//      int , NodeType
//      int, NodeType
//      dataType  -- persistance of pair
//      int  -- type (0/min, 1/saddle, 2/max)
//      dataType            -- scalar value at vertex 1
//      float, float, float -- vertex 1 coordinates
//      dataType            -- scalar value at vertex 2
//      float, float, float -- vertex 2 coordinates
template <typename dataType>
int BottleneckDistance::computeBottleneck(
    const std::vector<diagramTuple> &d1,
    const std::vector<diagramTuple> &d2,
    std::vector<matchingTuple> &matchings,
    const bool usePersistenceMetric)
{
  auto d1Size = (int) d1.size();
  auto d2Size = (int) d2.size();

  bool transposeOriginal = d1Size > d2Size;
  const std::vector<diagramTuple> &CTDiagram1 = transposeOriginal ? d2 : d1;
  const std::vector<diagramTuple> &CTDiagram2 = transposeOriginal ? d1 : d2;
  if (transposeOriginal) {
    int temp = d1Size;
    d1Size = d2Size;
    d2Size = temp;
  }

  if (transposeOriginal) {
    std::stringstream msg;
    msg << "[BottleneckDistance] The first persistence diagram is larger than the second." << std::endl;
    msg << "[BottleneckDistance] Solving the transposed problem." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  // Check user parameters.
  const int wasserstein = (wasserstein_ == "inf") ? -1 : stoi(wasserstein_);
  const std::string method = method_;
  if (wasserstein < 0 && wasserstein != -1) return -4;

  // Needed for distance computation.
  const double maxDistance =
      this->computeGeometricalRange<dataType>(CTDiagram1, CTDiagram2, d1Size, d2Size);

  // Needed to limit computation time.
  const dataType zeroThresh =
      this->computeMinimumRelevantPersistence<dataType>(CTDiagram1, CTDiagram2, d1Size, d2Size);

  // Initialize solvers.
  std::vector<matchingTuple> minMatchings;
  std::vector<matchingTuple> maxMatchings;
  std::vector<matchingTuple> sadMatchings;

  // Initialize cost matrices.
  int nbRowMin = 0, nbColMin = 0;
  int maxRowColMin = 0, minRowColMin = 0;
  int nbRowMax = 0, nbColMax = 0;
  int maxRowColMax = 0, minRowColMax = 0;
  int nbRowSad = 0, nbColSad = 0;
  int maxRowColSad = 0, minRowColSad = 0;

  // Remap for matchings.
  std::vector<int> minMap1; std::vector<int> minMap2;
  std::vector<int> maxMap1; std::vector<int> maxMap2;
  std::vector<int> sadMap1; std::vector<int> sadMap2;

  this->computeMinMaxSaddleNumberAndMapping(
      CTDiagram1, d1Size, nbRowMin, nbRowMax, nbRowSad,
      minMap1, maxMap1, sadMap1, zeroThresh);
  this->computeMinMaxSaddleNumberAndMapping(
      CTDiagram2, d2Size, nbColMin, nbColMax, nbColSad,
      minMap2, maxMap2, sadMap2, zeroThresh);

  // Automatically transpose if nb rows > nb cols
  maxRowColMin = std::max(nbRowMin+1, nbColMin+1);
  maxRowColMax = std::max(nbRowMax+1, nbColMax+1);
  maxRowColSad = std::max(nbRowSad+1, nbColSad+1);

  minRowColMin = std::min(nbRowMin+1, nbColMin+1);
  minRowColMax = std::min(nbRowMax+1, nbColMax+1);
  minRowColSad = std::min(nbRowSad+1, nbColSad+1);

  std::vector<std::vector<dataType>> minMatrix((unsigned long) minRowColMin, std::vector<dataType>(maxRowColMin));
  std::vector<std::vector<dataType>> maxMatrix((unsigned long) minRowColMax, std::vector<dataType>(maxRowColMax));
  std::vector<std::vector<dataType>> sadMatrix((unsigned long) minRowColSad, std::vector<dataType>(maxRowColSad));

  double px = px_;
  double py = py_;
  double pz = pz_;
  double pe = pe_;
  double ps = ps_;

  std::function<dataType (const diagramTuple, const diagramTuple)>
      distanceFunction =
      [maxDistance, wasserstein, px, py, pz, pe, ps]
          (const diagramTuple a, const diagramTuple b) -> dataType
      {
        BNodeType ta1 = std::get<1>(a);
        BNodeType ta2 = std::get<3>(a);
        double w = 1; // wasserstein > 1 ? wasserstein : 1;

        // We don't match critical points of different index.
        // This must be ensured before calling the distance function.
        // BNodeType tb1 = get<1>(b);
        // BNodeType tb2 = get<3>(b);
        bool isMin1 = ta1 == BLocalMin;
        bool isMax1 = ta2 == BLocalMax;
        // bool isBoth = isMin1 && isMax1;

        dataType rX = std::get<6>(a);
        dataType rY = std::get<10>(a);
        dataType cX = std::get<6>(b);
        dataType cY = std::get<10>(b);
        dataType x = (isMin1 ? pe : ps) * abs_diff<dataType>(rX, cX);
        dataType y = (isMax1 ? pe : ps) * abs_diff<dataType>(rY, cY);
        double geoDistance =
          isMax1 ?
          (
            px * pow(abs(std::get<11>(a) - std::get<11>(b)), w) +
            py * pow(abs(std::get<12>(a) - std::get<12>(b)), w) +
            pz * pow(abs(std::get<13>(a) - std::get<13>(b)), w)
          ) : isMin1 ?
          (
            px * pow(abs(std::get<7>(a) - std::get<7>(b)), w) +
            py * pow(abs(std::get<8>(a) - std::get<8>(b)), w) +
            pz * pow(abs(std::get<9>(a) - std::get<9>(b)), w)
          ) :
          (
            px * pow(abs(std::get<7>(a)+std::get<11>(a))/2 - abs(std::get<7>(b)+std::get<11>(b))/2, w) +
            py * pow(abs(std::get<8>(a)+std::get<12>(a))/2 - abs(std::get<8>(b)+std::get<12>(b))/2, w) +
            pz * pow(abs(std::get<9>(a)+std::get<13>(a))/2 - abs(std::get<9>(b)+std::get<13>(b))/2, w)
          );


        // dist /= maxDistance;
        double infDistance = std::max(x, y);
        double val = infDistance + geoDistance;
        val = pow(val, 1/w);
        return val;

      };

  std::function<dataType (const diagramTuple)>
      diagonalDistanceFunction =
      [wasserstein, px, py, pz, ps, pe]
          (const diagramTuple a) -> dataType
      {
        BNodeType ta1 = std::get<1>(a);
        BNodeType ta2 = std::get<3>(a);
        double w = 1; // wasserstein > 1 ? wasserstein : 1;
        bool isMin1 = ta1 == BLocalMin;
        bool isMax1 = ta2 == BLocalMax;

        dataType rX = std::get<6>(a);
        dataType rY = std::get<10>(a);

        double x1 = std::get<7>(a);
        double y1 = std::get<8>(a);
        double z1 = std::get<9>(a);
        double x2 = std::get<11>(a);
        double y2 = std::get<12>(a);
        double z2 = std::get<13>(a);

        double infDistance = (isMin1 || isMax1 ? pe : ps) * abs_diff<dataType>(rX, rY);
        double geoDistance =
            (px * pow(abs(x2 - x1), w) + py * pow(abs(y2 - y1), w) + pz * pow(abs(z2 - z1), w));
        double val = infDistance + geoDistance;
        return pow(val, 1/w);
      };

  const bool transposeMin = nbRowMin > nbColMin;
  const bool transposeMax = nbRowMax > nbColMax;
  const bool transposeSad = nbRowSad > nbColSad;

  Timer t;

  this->buildCostMatrices(
      CTDiagram1, CTDiagram2, d1Size, d2Size,
      distanceFunction, diagonalDistanceFunction,
      zeroThresh, minMatrix, maxMatrix, sadMatrix,
      transposeMin,
      transposeMax,
      transposeSad,
      wasserstein
  );


    if (nbRowMin > 0 && nbColMin > 0) {
      Munkres solverMin;
      dMsg(std::cout, "[BottleneckDistance] Affecting minima...\n", timeMsg);
      this->solvePWasserstein(
        minRowColMin, maxRowColMin,
        minMatrix, minMatchings, solverMin);
    }

    if (nbRowMax > 0 && nbColMax > 0) {
      Munkres solverMax;
      dMsg(std::cout, "[BottleneckDistance] Affecting maxima...\n", timeMsg);
      this->solvePWasserstein(
        minRowColMax, maxRowColMax,
        maxMatrix, maxMatchings, solverMax);
    }

    if (nbRowSad > 0 && nbColSad > 0) {
      Munkres solverSad;
      dMsg(std::cout, "[BottleneckDistance] Affecting saddles...\n", timeMsg);
      this->solvePWasserstein(
        minRowColSad, maxRowColSad,
        sadMatrix, sadMatchings, solverSad);
    }

  } else {

    // Launch solving for minima.
    if (nbRowMin > 0 && nbColMin > 0) {
      GabowTarjan solverMin;
      dMsg(std::cout, "[BottleneckDistance] Affecting minima...\n", timeMsg);
      this->solveInfinityWasserstein(
        minRowColMin, maxRowColMin, nbRowMin, nbColMin,
        minMatrix, minMatchings, solverMin);
    }

    // Launch solving for maxima.
    if (nbRowMax > 0 && nbColMax > 0) {
      GabowTarjan solverMax;
      dMsg(std::cout, "[BottleneckDistance] Affecting maxima...\n", timeMsg);
      this->solveInfinityWasserstein(
        minRowColMax, maxRowColMax, nbRowMax, nbColMax,
        maxMatrix, maxMatchings, solverMax);
    }

    // Launch solving for saddles.
    if (nbRowSad > 0 && nbColSad > 0) {
      GabowTarjan solverSad;
      dMsg(std::cout, "[BottleneckDistance] Affecting saddles...\n", timeMsg);
      this->solveInfinityWasserstein(
        minRowColSad, maxRowColSad, nbRowSad, nbColSad,
        sadMatrix, sadMatchings, solverSad);
    }
  }

  {
    std::stringstream msg;
    msg << "[BottleneckDistance] TTK CORE DONE IN "
        << t.getElapsedTime() << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  // Rebuild mappings.
  // Begin cost computation for unpaired vertices.

  // std::cout << "Min" << std::endl;
  dataType addedMinPersistence = this->buildMappings<dataType>(
    minMatchings, transposeOriginal, transposeMin,
    matchings, minMap1, minMap2, wasserstein);

  // std::cout << "Max" << std::endl;
  dataType addedMaxPersistence = this->buildMappings<dataType>(
    maxMatchings, transposeOriginal, transposeMax,
    matchings, maxMap1, maxMap2, wasserstein);

  // std::cout << "Sad" << std::endl;
  dataType addedSadPersistence = this->buildMappings<dataType>(
    sadMatchings, transposeOriginal, transposeSad,
    matchings, sadMap1, sadMap2, wasserstein);

  // TODO [HIGH] do that for embeddings
  // Recompute matching weights for user-friendly distance.

  dataType d = 0;
  std::vector<bool> paired1(d1Size);
  std::vector<bool> paired2(d2Size);
  for (int b = 0; b < d1Size; ++b) paired1[b] = false;
  for (int b = 0; b < d2Size; ++b) paired2[b] = false;

  int numberOfMismatches = 0;
  for (int m = 0, ms = (int) matchings.size(); m < ms; ++m)
  {
    matchingTuple t = matchings[m];
    int i = transposeOriginal ? std::get<1>(t) : std::get<0>(t);
    int j = transposeOriginal ? std::get<0>(t) : std::get<1>(t);
    //dataType val = std::get<2>(t);

    diagramTuple t1 = CTDiagram1[i];
    diagramTuple t2 = CTDiagram2[j];
    //dataType rX = std::get<6>(t1); dataType rY = std::get<10>(t1);
    //dataType cX = std::get<6>(t2); dataType cY = std::get<10>(t2);
    //dataType x = rX - cX; dataType y = rY - cY;
    paired1[i] = true;
    paired2[j] = true;

    //dataType lInf = std::max(abs<dataType>(x), abs<dataType>(y));

    //if (((wasserstein < 0 && lInf != val) || (wasserstein > 0 && pow(lInf, wasserstein) != val)))
      //++numberOfMismatches;

    dataType partialDistance = distanceFunction(t1, t2);
    //wasserstein > 0 ? pow(lInf, wasserstein) : std::max(d, lInf);

    if (wasserstein > 0) d += partialDistance;
    else d = partialDistance;
  }

  if (numberOfMismatches > 0) {
    std::stringstream msg;
    msg << "[BottleneckDistance] Distance mismatch when rebuilding " << numberOfMismatches << " matchings." << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  dataType affectationD = d;
  d = wasserstein > 0 ?
      pow(d + addedMaxPersistence + addedMinPersistence + addedSadPersistence, (1.0 / (double) wasserstein)) :
      std::max(d, std::max(addedMaxPersistence, std::max(addedMinPersistence, addedSadPersistence)));

  {
    std::stringstream msg;
    msg << "[BottleneckDistance] Computed distance diagMax("
        << addedMaxPersistence << "), diagMin("
        << addedMinPersistence << "), diagSad("
        << addedSadPersistence << "), affAll(" << affectationD << "), res(" << d << ")" << std::endl;
    dMsg(std::cout, msg.str(), timeMsg);
  }

  distance_ = (double) d;
  return 0;
}




template <typename dataType>
int ttk::BottleneckDistance::computeAuction(
  const std::vector<diagramTuple> *CTDiagram1,
  const std::vector<diagramTuple> *CTDiagram2,
  std::vector<matchingTuple> *matchings,
  const double alpha,
  const double delta_lim
)
{
	const double geometricalFactor = (alpha < 0.0 || alpha > 1.0) ? 1.0 : alpha;
	const int wasserstein = (wasserstein_ == "inf") ? -1 : stoi(wasserstein_);
	auto* distance = new dataType;

	const int d1Size = (int) CTDiagram1->size();
	const int d2Size = (int) CTDiagram2->size();
	const dataType zeroThresh = this->computeMinimumRelevantPersistence(CTDiagram1, CTDiagram2, d1Size, d2Size);

	std::vector<diagramTuple> D1Min, D1Sad, D1Max, D2Min, D2Sad, D2Max;
	std::vector<int> D1Min_idx, D1Sad_idx, D1Max_idx, D2Min_idx, D2Sad_idx, D2Max_idx;
	// Creating the subdiagrams for the three types of pairs
	// Creating a table to identify for each of these subdiagram point, its corresponding point in the original diagram
	for(int i=0; i < (int) CTDiagram1->size(); ++i){
		diagramTuple t = CTDiagram1->at(i);
		BNodeType nt1 = std::get<1>(t);
		BNodeType nt2 = std::get<3>(t);
    if( std::get<6>(t) > std::get<10>(t) ){
      std::cout << "\n###########\n ERROR : wrong match \n###########\n" << '\n';
      return -1;
    }
		dataType dt = std::get<4>(t);
		if (abs<dataType>(dt) < zeroThresh) continue;

		if (nt1 == BLocalMin && nt2 == BLocalMax) {
			D1Max.push_back(t);
			D1Max_idx.push_back(i);
		}
		else {
			if (nt1 == BLocalMax || nt2 == BLocalMax) {
				D1Max.push_back(t);
				D1Max_idx.push_back(i);
			}
			if (nt1 == BLocalMin || nt2 == BLocalMin) {
				D1Min.push_back(t);
				D1Min_idx.push_back(i);
			}
			if ((nt1 == BSaddle1 && nt2 == BSaddle2)
				|| (nt1 == BSaddle2 && nt2 == BSaddle1)) {
				D1Sad.push_back(t);
				D1Sad_idx.push_back(i);
			}
		}
	}

	for(int i=0; i < (int) CTDiagram2->size(); ++i){
		diagramTuple t = CTDiagram2->at(i);
		BNodeType nt1 = std::get<1>(t);
		BNodeType nt2 = std::get<3>(t);
		dataType dt = std::get<4>(t);
		if (abs<dataType>(dt) < zeroThresh) continue;

		if (nt1 == BLocalMin && nt2 == BLocalMax) {
			D2Max.push_back(t);
			D2Max_idx.push_back(i);
		}
		else {
			if (nt1 == BLocalMax || nt2 == BLocalMax) {
				D2Max.push_back(t);
				D2Max_idx.push_back(i);
			}
			if (nt1 == BLocalMin || nt2 == BLocalMin) {
				D2Min.push_back(t);
				D2Min_idx.push_back(i);
			}
			if ((nt1 == BSaddle1 && nt2 == BSaddle2)
				|| (nt1 == BSaddle2 && nt2 == BSaddle1)) {
				D2Sad.push_back(t);
				D2Sad_idx.push_back(i);
			}
		}
	}


	std::vector<matchingTuple> minMatchings;
	std::vector<matchingTuple> maxMatchings;
	std::vector<matchingTuple> sadMatchings;
	dataType d = 0;

	#ifdef TTK_ENABLE_OPENMP
	#pragma omp parallel sections
	#endif
	{
		#ifdef TTK_ENABLE_OPENMP
		#pragma omp section
		#endif
		{
			if(D1Min.size()+D2Min.size()>0){
				Auction<dataType> auctionMin(wasserstein, geometricalFactor, delta_lim, use_kdtree_);
				dMsg(std::cout, "[BottleneckDistance] Affecting minima...\n", timeMsg);
				auctionMin.BuildAuctionDiagrams(D1Min, D2Min);
				dataType cost = auctionMin.run(&minMatchings);
				std::stringstream msg;
				msg << "[Auction] Total cost for minima = " << cost << std::endl;
				dMsg(std::cout, msg.str(), timeMsg);
				d += cost;
			}
		}

		#ifdef TTK_ENABLE_OPENMP
		#pragma omp section
		#endif
		{
			if(D1Sad.size()+D2Sad.size()>0){
				Auction<dataType> auctionSad(wasserstein, geometricalFactor, delta_lim, use_kdtree_);
				dMsg(std::cout, "[BottleneckDistance] Affecting saddles...\n", timeMsg);
				auctionSad.BuildAuctionDiagrams(D1Sad, D2Sad);
				dataType cost = auctionSad.run(&sadMatchings);
				std::stringstream msg;
				msg << "[Auction] Total cost for saddles = " << cost << std::endl;
				dMsg(std::cout, msg.str(), timeMsg);
				d += cost;
			}
		}

		#ifdef TTK_ENABLE_OPENMP
		#pragma omp section
		#endif
		{
			if(D1Max.size()+D2Max.size()>0){
				Auction<dataType> auctionMax(wasserstein, geometricalFactor, delta_lim, use_kdtree_);
				dMsg(std::cout, "[BottleneckDistance] Affecting maxima...\n", timeMsg);
				auctionMax.BuildAuctionDiagrams(D1Max, D2Max);
				dataType cost = auctionMax.run(&maxMatchings);
				std::stringstream msg;
				msg << "[Auction] Total cost for maxima = " << cost << std::endl;
				dMsg(std::cout, msg.str(), timeMsg);
				d += cost;
			}
		}
	}
	// Construct the matching vectors as the concatenation of all three matchings
	for(int k=0; k < (int) minMatchings.size(); k++){
		matchingTuple t = minMatchings[k];
		int i = std::get<0>(t);
		int j = std::get<1>(t);
		dataType cost = std::get<2>(t);
		matchings->push_back(std::make_tuple(D1Min_idx[i], D2Min_idx[j], cost));
	}

	for(int k=0; k < (int) maxMatchings.size(); k++){
		matchingTuple t = maxMatchings[k];
		int i = std::get<0>(t);
		int j = std::get<1>(t);
		dataType cost = std::get<2>(t);
		matchings->push_back(std::make_tuple(D1Max_idx[i], D2Max_idx[j], cost));
	}

	for(int k=0; k < (int) sadMatchings.size(); k++){
		matchingTuple t = sadMatchings[k];
		int i = std::get<0>(t);
		int j = std::get<1>(t);
		dataType cost = std::get<2>(t);
		matchings->push_back(std::make_tuple(D1Sad_idx[i], D2Sad_idx[j], cost));

	}


	*distance = pow(d, (1.0 / (double) wasserstein));
	distance_ = (void*)(distance);
	return 0;
  }


#endif
