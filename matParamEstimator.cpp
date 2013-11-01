#include "matParamEstimator.h"

MatParamEstimator::MatParamEstimator(double delta):m_deltaSquared(delta*delta)
{

}

void MatParamEstimator::estimate(std::vector<SamePoint> &data, std::vector<double> &parameters)
{
	parameters.clear();
	if (data.size() < 3)
	{
		return;
	}
	double a, b, c, d, e, f;
	a =((data[0].rx-data[2].rx)*(data[1].ly-data[2].ly)/(data[0].ly-data[2].ly)-(data[1].rx-data[2].rx))/
		((data[0].lx-data[2].lx)*(data[1].ly-data[2].ly)/(data[0].ly-data[2].ly)-(data[1].lx-data[2].lx));
	b = ((data[0].rx-data[2].rx)-a*(data[0].lx-data[2].lx))/(data[0].ly-data[2].ly);
	c = data[0].rx-a*data[0].lx-b*data[0].ly;
	d =((data[0].ry-data[2].ry)*(data[1].ly-data[2].ly)/(data[0].ly-data[2].ly)-(data[1].ry-data[2].ry))/
		((data[0].lx-data[2].lx)*(data[1].ly-data[2].ly)/(data[0].ly-data[2].ly)-(data[1].lx-data[2].lx));
	e = ((data[0].ry-data[2].ry)-d*(data[0].lx-data[2].lx))/(data[0].ly-data[2].ly);
	f = data[0].ry-d*data[0].lx-e*data[0].ly;
	parameters.push_back(a);
	parameters.push_back(b);
	parameters.push_back(c);
	parameters.push_back(d);
	parameters.push_back(e);
	parameters.push_back(f);
}

void MatParamEstimator::leastSquaresEstimate(std::vector<SamePoint> &data, std::vector<double> &parameters)
{
	parameters.clear();
	double ratio[3][3] = {0};
	double resultx[3] = {0};
	double resulty[3] = {0};
	double a0, a1, a2, b0, b1, b2;
	std::vector<SamePoint>::iterator temIte = data.begin();
	int i = 0;
	while(temIte != data.end())
	{
		ratio[0][1] += temIte->lx;
		ratio[0][2] += temIte->ly;
		ratio[1][0] += temIte->lx;
		ratio[1][1] += temIte->lx*temIte->lx;
		ratio[1][2] += temIte->lx*temIte->ly;
		ratio[2][0] += temIte->ly;
		ratio[2][1] += temIte->lx*temIte->ly;
		ratio[2][2] += temIte->ly*temIte->ly;
		resultx[0] += temIte->rx;
		resultx[1] += temIte->lx*temIte->rx;
		resultx[2] += temIte->ly*temIte->rx;
		resulty[0] += temIte->ry;
		resulty[1] += temIte->lx*temIte->ry;
		resulty[2] += temIte->ly*temIte->ry;
		++i;
		++temIte;
	}
	ratio[0][0] = i;

	for (int j = 0; j < 2; ++j)
	{
		resultx[j] = resultx[j]*(-ratio[2][0]/ratio[j][0])+resultx[2];
		resulty[j] = resulty[j]*(-ratio[2][0]/ratio[j][0])+resulty[2];
		for (int i = 2; i >= 0; --i)
		{
			ratio[j][i] = ratio[j][i]*(-ratio[2][0]/ratio[j][0]) + ratio[2][i];
		}
	}
	a2 = (resultx[0]*(-ratio[1][1]/ratio[0][1])+resultx[1])/(ratio[0][2]*(-ratio[1][1]/ratio[0][1])+ratio[1][2]);
	a1 = (resultx[1]-a2*ratio[1][2])/ratio[1][1];
	a0 = (resultx[2]-a2*ratio[2][2]-a1*ratio[2][1])/ratio[2][0];


	b2 = (resulty[0]*(-ratio[1][1]/ratio[0][1])+resulty[1])/(ratio[0][2]*(-ratio[1][1]/ratio[0][1])+ratio[1][2]);
	b1 = (resulty[1]-b2*ratio[1][2])/ratio[1][1];
	b0 = (resulty[2]-b2*ratio[2][2]-b1*ratio[2][1])/ratio[2][0];

	parameters.push_back(a1);
	parameters.push_back(a2);
	parameters.push_back(a0);
	parameters.push_back(b1);
	parameters.push_back(b2);
	parameters.push_back(b0);
}

bool MatParamEstimator::agree(std::vector<double> &parameters, SamePoint &data)
{
	double dx = data.lx*parameters[0]+data.ly*parameters[1]+parameters[2] - data.rx;
	double dy = data.lx*parameters[3]+data.ly*parameters[4]+parameters[5]-data.ry;
	double distance = dx*dx+dy*dy;
	return (distance<m_deltaSquared);
}