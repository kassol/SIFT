#pragma once
#include "parameteresitmator.h"
#include "sift.h"

class MatParamEstimator : public ParameterEsitmator<SamePoint, double>
{
public:
	MatParamEstimator(double delta);

	virtual void estimate(std::vector<SamePoint> &data, std::vector<double> &parameters);

	virtual void leastSquaresEstimate(std::vector<SamePoint> &data, std::vector<double> &parameters);

	virtual bool agree(std::vector<double> &parameters, SamePoint &data);
private:
	double m_deltaSquared;
};

