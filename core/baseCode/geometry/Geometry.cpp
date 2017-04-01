#include                <Geometry.h>

// TODO:
// for faster computations, remove all pow10

Geometry::Geometry(){ 
  
}

Geometry::~Geometry(){
}

double Geometry::angle(
  const double *vA0, const double *vA1, 
  const double *vB0, const double *vB1){

  return M_PI - acos(
    dotProduct(vA0, vA1, vB0, vB1)/
    (magnitude(vA0, vA1)*magnitude(vB0, vB1)));
}

bool Geometry::areVectorsColinear(const double *vA0, const double *vA1, 
  const double *vB0, const double *vB1, 
  vector<double> *coefficients,
  const double *tolerance){

  int aNullComponents = 0, bNullComponents = 0;
  vector<double> a(3), b(3);
  for(int i = 0; i < 3; i++){
    a[i] = vA1[i] - vA0[i];
    if(fabs(a[i]) < pow10(-FLT_DIG))
      aNullComponents ++;
    b[i] = vB1[i] - vB0[i];
    if(fabs(b[i]) < pow10(-FLT_DIG))
      bNullComponents ++;
  }
  
  if((aNullComponents == 3)||(bNullComponents == 3))
    return true;
  
  // check for axis aligned vectors
  if((aNullComponents > 1)||(bNullComponents > 1)){
    if(aNullComponents == bNullComponents){
      // only one non-null component for both vectors
      return true;
    }
  }
 
  bool useDenominatorA = false;
  double sumA = 0, sumB = 0;
  for(int i = 0; i < 3; i++){
    sumA += fabs(a[i]);
    sumB += fabs(b[i]);
  }
  if(sumA > sumB){
    useDenominatorA = true;
  }
  
  vector<double> k(3, 0);
 
  double maxDenominator = 0;
  int isNan = -1, maximizer = 0;
  for(int i = 0; i < 3; i++){
    if(useDenominatorA){
      if(fabs(a[i]) > pow10(-FLT_DIG)){
        k[i] = b[i]/a[i];
      }
      else{
        isNan = i;
      } 
    }
    else{
      if(fabs(b[i]) > pow10(-FLT_DIG)){
        k[i] = a[i]/b[i];
      }
      else{
        isNan = i;
      }
    }
    
    if(!i){
      maxDenominator = fabs(k[i]);
      maximizer = i;
    }
    else{
      if(fabs(k[i]) > maxDenominator){
        maxDenominator = fabs(k[i]);
        maximizer = i;
      }
    }
  }
  
  double colinearityThreshold;
  
  colinearityThreshold = pow10(-FLT_DIG);
  if(tolerance)
    colinearityThreshold = *tolerance;
  
  if(coefficients){
    (*coefficients) = k;
  }
  
  if(isNan == -1){
   
    if((fabs(1 - fabs(k[(maximizer+1)%3]/k[maximizer])) < colinearityThreshold)
      &&(fabs(1 - fabs(k[(maximizer+2)%3]/k[maximizer])) 
        < colinearityThreshold)){
      return true;
    }
  }
  else{
    if(fabs(1 - fabs(k[(isNan+1)%3]/k[(isNan+2)%3])) < colinearityThreshold)
      return true;
  }

  k[0] = k[1] = k[2] = 0;
  
  return false;
}

int Geometry::computeBarycentricCoordinates(
  const double *p0, const double *p1, const double *p, 
  vector<double> &baryCentrics, const int &dimension){

  baryCentrics.resize(2);
  
  int bestI = 0;
  double maxDenominator = 0;
  
  for(int i = 0; i < dimension; i++){
    
    double denominator = fabs(p0[i] - p1[i]);
    if(!i){
      maxDenominator = denominator;
      bestI = i;
    }
    else{
      if(denominator > maxDenominator){
        maxDenominator = denominator;
        bestI = i;
      }
    }
  }
  
  baryCentrics[0] = p0[bestI] - p1[bestI];
  baryCentrics[0] = (p[bestI] - p1[bestI])/baryCentrics[0];
  
  baryCentrics[1] = 1 - baryCentrics[0];
  
  
  // check if the point lies in the edge
  vector<float> test(dimension);
  for(int i = 0; i < dimension; i++)
    test[i] = 
      baryCentrics[0]*p0[i] + baryCentrics[1]*p1[i];
     
  if((!((fabs(test[0] - p[0]) < pow(10, -FLT_DIG+1))
    &&(fabs(test[1] - p[1]) < pow(10, -FLT_DIG+1))))
    ){
    for(int i = 0; i < 2; i++) baryCentrics[i] = -baryCentrics[i];
  }
  
  return 0;
}


int Geometry::computeBarycentricCoordinates(
  const double *p0, const double *p1, const double *p2, const double *p, 
  vector<double> &baryCentrics){

  baryCentrics.resize(3);
 
  // find the pair of coordinates that maximize the sum of the denominators
  // (more stable computations)
  int bestI = 0, bestJ = 1;
  double maxDenominator = 0;
  
  for(int i = 0; i < 2; i++){
    for(int j = i + 1; j < 3; j++){
      
      baryCentrics[0] = (p1[j] - p2[j])*(p0[i] - p2[i])
        + (p2[i] - p1[i])*(p0[j] - p2[j]);
      baryCentrics[1] = (p1[j] - p2[j])*(p0[i] - p2[i])
        + (p2[i] - p1[i])*(p0[j] - p2[j]);
        
      double denominator = fabs(baryCentrics[0]);
      
      if(fabs(baryCentrics[1]) < denominator)
        denominator = fabs(baryCentrics[1]);
        
      if((i == 0)&&(j == 1)){
        maxDenominator = denominator;
      }
      else{
        if(denominator > maxDenominator){
          maxDenominator = denominator;
          bestI = i;
          bestJ = j;
        }
      }
    }
  }
  
  baryCentrics[0] = (p1[bestJ] - p2[bestJ])*(p0[bestI] - p2[bestI])
    + (p2[bestI] - p1[bestI])*(p0[bestJ] - p2[bestJ]);
  // (y1 - y2)*(x - x2) + (x2 - x1)*(y - y2)
  baryCentrics[0] = ((p1[bestJ] - p2[bestJ])*(p[bestI] - p2[bestI]) 
    + (p2[bestI] - p1[bestI])*(p[bestJ] - p2[bestJ]))/baryCentrics[0];
  
  // (y1 - y2)*(x0 - x2) + (x2 - x1)*(y0 - y2)
  baryCentrics[1] = (p1[bestJ] - p2[bestJ])*(p0[bestI] - p2[bestI])
    + (p2[bestI] - p1[bestI])*(p0[bestJ] - p2[bestJ]);
  // (y2 - y0)*(x - x2) + (x0 - x2)*(y - y2)
  baryCentrics[1] = ((p2[bestJ] - p0[bestJ])*(p[bestI] - p2[bestI])
    + (p0[bestI] - p2[bestI])*(p[bestJ] - p2[bestJ]))/baryCentrics[1];
    
  baryCentrics[2] = 1 - baryCentrics[0] - baryCentrics[1];
  
  // check if the point lies in the triangle
  vector<float> test(3);
  for(int i = 0; i < 3; i++)
    test[i] = 
      baryCentrics[0]*p0[i] + baryCentrics[1]*p1[i] + baryCentrics[2]*p2[i];
     
  if(!((fabs(test[0] - p[0]) < pow(10, -FLT_DIG))
    &&(fabs(test[1] - p[1]) < pow(10  , -FLT_DIG))
    &&(fabs(test[2] - p[2]) < pow(10, -FLT_DIG)))){
    for(int i = 0; i < 3; i++) baryCentrics[i] = -1 - baryCentrics[i];
  }
  
  return 0;
}

int Geometry::computeBarycentricCoordinates(
  const float *p0, const float *p1, const float *p2, const float *p, 
  vector<double> &baryCentrics){

  vector<double> P0(3), P1(3), P2(3), P(3);
  
  for(int i = 0; i < 3; i++){
    P0[i] = p0[i];
    P1[i] = p1[i];
    P2[i] = p2[i];
    P[i] = p[i];
  }
  
  return computeBarycentricCoordinates(
    P0.data(), P1.data(), P2.data(), P.data(), baryCentrics);
}

bool Geometry::computeSegmentIntersection(
  const double &xA, const double &yA, const double &xB, const double &yB,
  const double &xC, const double &yC, const double &xD, const double &yD,
  double &x, double &y){

  double d = (xA - xB)
    *(yC - yD) 
    - 
    (yA - yB)
    *(xC- xD);
 
  if(fabs(d) < pow(10, -DBL_DIG))
    return false;
  
  x = 
    ((xC - xD)
      *(xA*yB 
        - yA*xB)
      - (xA - xB)
        *(xC*yD 
          - yC*xD))/d;
          
  y = 
    ((yC - yD)
      *(xA*yB
        - yA*xB)
      - (yA - yB)
        *(xC*yD
          - yC*xD))/d;
  
  if((x < std::min(xA, xB) - pow10(-FLT_DIG))
    ||(x > std::max(xA, xB) + pow10(-FLT_DIG)))
    return false;
          
  if((x < std::min(xC, xD) - pow10(-FLT_DIG))
    ||(x > std::max(xC, xD) + pow10(-FLT_DIG)))
    return false;
    
  return true;
}

int Geometry::computeTriangleArea(
  const double *p0, const double *p1, const double *p2, 
  double &area){
  
  vector<double> cross;
  
  crossProduct(p0, p1, p1, p2, cross);
  
  area = 0.5*magnitude(cross.data());
  
  return 0;
}

int Geometry::computeTriangleAngles(
  const double *p0, const double *p1, const double *p2, 
  vector<double> &angles){

  angles.resize(3);
 
  angles[0] = angle(p0, p1, p1, p2);
  angles[1] = angle(p1, p2, p2, p0);
  angles[2] = angle(p2, p0, p0, p1);
  
  return 0;
}

int Geometry::crossProduct(const double *vA0, const double *vA1,
  const double *vB0, const double *vB1,
  vector<double> &crossProduct){
  
  crossProduct.resize(3);
 
  vector<double> a(3), b(3);
  
  for(int i = 0; i < 3; i++){
    a[i] = vA1[i] - vA0[i];
    b[i] = vB1[i] - vB0[i];
  }
  
  for(int i = 0; i < 3; i++){
    crossProduct[i] = 
      a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3];
  }
  
  return 0;
}

int Geometry::crossProduct(const double *vA, const double *vB,
		double *crossProduct){
	crossProduct[0]=vA[1]*vB[2]-vA[2]*vB[1];
	crossProduct[1]=vA[2]*vB[0]-vA[0]*vB[2];
	crossProduct[2]=vA[0]*vB[1]-vA[1]*vB[0];
	return 0;
}

double Geometry::distance(const double *p0, const double *p1, 
  const int &dimension){

  double distance = 0;
  
  for(int i = 0; i < dimension; i++){
    distance += (p0[i] - p1[i])*(p0[i] - p1[i]);
  }
  
  return sqrt(distance);
}

double Geometry::distance(const float *p0, const float *p1, 
  const int &dimension){

  double distance = 0;
  
  for(int i = 0; i < dimension; i++){
    distance += (p0[i] - p1[i])*(p0[i] - p1[i]);
  }
  
  return sqrt(distance);
}

double Geometry::dotProduct(
  const double *vA0, const double *vA1, 
  const double *vB0, const double *vB1){

  double dotProduct = 0;
  for(int i = 0; i < 3; i++)
    dotProduct += (vA1[i] - vA0[i])*(vB1[i] - vB0[i]);
    
  return dotProduct;
}

double Geometry::dotProduct(const double *vA,
		const double *vB){
		return vA[0]*vB[0]+vA[1]*vB[1]+vA[2]*vB[2];
}

int Geometry::getBoundingBox(const vector<vector<float> > &points, 
  vector<pair<double, double>> &bBox){

  if(!points.size())
    return -1;
  
  int dimension = points[0].size();
  
  bBox.resize(dimension);
  
  for(int i = 0; i < (int) points.size(); i++){
      
    if(!i){
      for(int j = 0; j < dimension; j++){
        bBox[j].first = points[i][j];
        bBox[j].second = points[i][j];
      }
    }
    else{
      for(int j = 0; j < dimension; j++){
        if(points[i][j] < bBox[j].first){
          bBox[j].first = points[i][j];
        }
        if(points[i][j] > bBox[j].second){
          bBox[j].second = points[i][j];
        }
      }
    }
  }
  
  return 0;
}

bool Geometry::isPointInTriangle(const double *p0, 
  const double *p1, const double *p2, const double *p){

  vector<double> barycentrics;
  
  Geometry::computeBarycentricCoordinates(p0, p1, p2, p, barycentrics);
  
  for(int i = 0; i < (int) barycentrics.size(); i++){
    if(barycentrics[i] < -pow10(-DBL_DIG))
      return false;
    if(barycentrics[i] > 1 + pow10(-DBL_DIG))
      return false;
  }
    
  return true;
}

bool Geometry::isPointInTriangle(const float *p0, 
  const float *p1, const float *p2, const float *p){

  vector<double> barycentrics;
  
  Geometry::computeBarycentricCoordinates(p0, p1, p2, p, barycentrics);
  
  for(int i = 0; i < (int) barycentrics.size(); i++){
    if(barycentrics[i] < -pow10(-FLT_DIG))
      return false;
    if(barycentrics[i] > 1 + pow10(-FLT_DIG))
      return false;
  }
    
  return true;
}

bool Geometry::isPointOnSegment(const double &x, const double &y,
  const double &xA, const double &yA, const double &xB, const double &yB){
    
  
  vector<double> pA(2), pB(2), p(2);
  
  pA[0] = xA;
  pA[1] = yA;
  
  pB[0] = xB;
  pB[1] = yB;
  
  p[0] = x;
  p[1] = y;
  
  
  return Geometry::isPointOnSegment(p.data(), pA.data(), pB.data(), 2);
}

bool Geometry::isPointOnSegment(const double *p, 
  const double *pA, const double *pB, const int &dimension){

  vector<double> baryCentrics;
  
  Geometry::computeBarycentricCoordinates(pA, pB, p,
    baryCentrics, dimension);
  
  return (((baryCentrics[0] > -pow10(-DBL_DIG))
    &&(baryCentrics[0] < 1 + pow10(-DBL_DIG)))
    &&
    ((baryCentrics[1] > -pow10(-DBL_DIG))
    &&(baryCentrics[1] < 1 + pow10(-DBL_DIG))));
}


bool Geometry::isTriangleColinear(
  const double *p0, const double *p1, const double *p2, 
  const double *tolerance){
 
  bool maxDecision = false;
  double maxCoefficient = 0;
  vector<double> coefficients(3);
  
  bool decision = areVectorsColinear(p0, p1, p1, p2, &coefficients, tolerance);
  maxDecision = decision;
  for(int i = 0; i < 3; i++){
    if(!i){
      maxCoefficient = fabs(coefficients[i]);
      maxDecision = decision;
    }
    else{
      if(fabs(coefficients[i]) > maxCoefficient){
        maxCoefficient = fabs(coefficients[i]);
        maxDecision = decision;
      }
    }
  }
  
  decision = areVectorsColinear(p0, p2, p2, p1, &coefficients, tolerance);
  for(int i = 0; i < 3; i++){
    if(fabs(coefficients[i]) > maxCoefficient){
      maxCoefficient = fabs(coefficients[i]);
      maxDecision = decision;
    }
  }
  
  decision = areVectorsColinear(p1, p0, p0, p2, &coefficients, tolerance);
  for(int i = 0; i < 3; i++){
    if(fabs(coefficients[i]) > maxCoefficient){
      maxCoefficient = fabs(coefficients[i]);
      maxDecision = decision;
    }
  }
  
  return maxDecision;
}

double Geometry::magnitude(const double *v){
  
  double mag = 0;
  
  for(int i = 0; i < 3; i++){
    mag += v[i]*v[i];
  }
  
  return sqrt(mag);
}

double Geometry::magnitude(const double *o, const double *d){

  double mag = 0;
  
  for(int i = 0; i < 3; i++){
    mag += (o[i] - d[i])*(o[i] - d[i]);
  }
  
  return sqrt(mag);
}
