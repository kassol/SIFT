#pragma once

#include <set>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <time.h>



#include "parameterEsitmator.h"

template<class T, class S>
class Ransac
{
public:
	static double compute(std::vector<S> &parameters, ParameterEsitmator<T, S>* paramEstimator,
		std::vector<T> &data, int numForEstimate, double desiredProbabilityForNoOutliers,
		double maximalOutlierPercentage, std::vector<T>& outData);

	static double compute(std::vector<S> &parameters, ParameterEsitmator<T, S>* paramEstimator,
		std::vector<T> &data, int numForEstimate, std::vector<T>&outdata);

private:
	static void computeAllChoices(ParameterEsitmator<T, S>* paramEstimator, std::vector<T> &data, 
		int numForEstimate, short* bestVotes, short* curVotes, int &numVotesForBest, int startIndex,
		int n, int k, int arrIndex, int* arr);

	static void estimate(ParameterEsitmator<T, S>* paramEstimator, std::vector<T> &data, 
		int numForEstimate, short* bestVotes, short* curVotes, int &numVotesForBest, int* arr);

	class SubSetIndexComparator
	{
	public:
		SubSetIndexComparator(int arrayLength) : m_length(arrayLength){}
		bool operator()(const int *arr1, const int *arr2) const
		{
			for(int i=0; i<m_length; i++)
			{
				if(arr1[i] < arr2[i])
				{
					return true;
				}
			}
			return false;			
		}
	private:
		int m_length;
	};
};



template<class T, class S>
double Ransac<T, S>::compute(std::vector<S> &parameters, ParameterEsitmator<T, S>* paramEstimator, 
	std::vector<T> &data, int numForEstimate, double desiredProbabilityForNoOutliers,
	double maximalOutlierPercentage, std::vector<T>& outData)
{
	int numDataObjects = data.size();
	if (numDataObjects < numForEstimate || maximalOutlierPercentage >= 1.0)
	{
		return 0;
	}
	std::vector<T> exactEstimateData;
	std::vector<S> exactEstimateParameters;
	int i, j, k, l, numVotesForBest, numVotesForCur, maxIndex, numTries;
	short* bestVotes = new short[numDataObjects];
	short* curVotes = new short[numDataObjects];
	short* notChosen = new short[numDataObjects];
	SubSetIndexComparator subSetIndexComparator(numForEstimate);
	std::set<int*, SubSetIndexComparator> chosenSubSets(subSetIndexComparator);
	int* curSubSetIndexes;
	double outlierPercentage = maximalOutlierPercentage;
	double numerator = log(1.0-desiredProbabilityForNoOutliers);
	double denominator = log(1- pow(1-maximalOutlierPercentage, numForEstimate));

	parameters.clear();

	numVotesForBest = -1;
	srand((unsigned int)time(NULL));
	numTries = (int)(numerator/denominator+0.5);
	for (i = 0; i < numTries; ++i)
	{
		memset(notChosen, '1', numDataObjects*sizeof(short));
		curSubSetIndexes = new int[numForEstimate];
		exactEstimateData.clear();

		maxIndex = numDataObjects-1;
		for (l = 0; l < numForEstimate; ++l)
		{
			int selectedIndex = (int)(((float)rand()/(float)RAND_MAX)*maxIndex + 0.5);
			for (j = -1, k = 0; k < numDataObjects && j < selectedIndex; ++k)
			{
				if (notChosen[k])
				{
					++j;
				}
			}
			--k;
			exactEstimateData.push_back(data[k]);
			notChosen[k] = 0;
			--maxIndex;
		}
		for (l = 0, j = 0; j < numDataObjects; ++j)
		{
			if (!notChosen[j])
			{
				curSubSetIndexes[l] = j+1;
				++l;
			}
		}

		std::pair< std::set<int *, SubSetIndexComparator >::iterator, bool > res = chosenSubSets.insert(curSubSetIndexes);

		if (res.second == true)
		{
			paramEstimator->estimate(exactEstimateData, exactEstimateParameters);

			numVotesForCur = 0;
			memset(curVotes, '\0', numDataObjects*sizeof(short));
			for (j = 0; j < numDataObjects; ++j)
			{
				if (paramEstimator->agree(exactEstimateParameters, data[j]))
				{
					curVotes[j] = 1;
					++numVotesForCur;
				}
			}
			if (numVotesForCur > numVotesForBest)
			{
				numVotesForBest = numVotesForCur;
				memcpy(bestVotes, curVotes, numDataObjects*sizeof(short));
			}

			outlierPercentage = 1-(double)numVotesForCur/(double)numDataObjects;
			if (outlierPercentage < maximalOutlierPercentage)
			{
				maximalOutlierPercentage = outlierPercentage;
				denominator = log(1-pow(1-maximalOutlierPercentage, numForEstimate));
				numTries = (int)(numerator/denominator+0.5);
			}
		}
		else
		{
			delete []curSubSetIndexes;
			--i;
		}
	}
	std::set<int *, SubSetIndexComparator >::iterator it = chosenSubSets.begin();
	while (it != chosenSubSets.end())
	{
		delete [](*it);
		++it;
	}
	chosenSubSets.clear();

	for (j = 0; j < numDataObjects; ++j)
	{
		if (bestVotes[j])
		{
			outData.push_back(data[j]);
		}
	}
	paramEstimator->leastSquaresEstimate(outData, parameters);

	delete []bestVotes;
	delete []curVotes;
	delete []notChosen;

	return (double)numVotesForBest/(double)numDataObjects;
}

template<class T, class S>
double Ransac<T, S>::compute(std::vector<S> &parameters, ParameterEsitmator<T, S>* paramEstimator, 
	std::vector<T> &data, int numForEstimate, std::vector<T>&outdata)
{
	int numDataObjects = data.size();
	int numVotesForBest = -1;
	int* arr = new int[numDataObjects];
	short* curVotes = new short[numDataObjects];
	short* bestVotes = new short[numDataObjects];

	if (numDataObjects < numForEstimate)
	{
		return 0;
	}
	computeAllChoices(paramEstimator, data, numForEstimate, bestVotes, curVotes, numVotesForBest, 0,
		data.size(), numForEstimate, 0, arr);

	for (int j = 0; j < numDataObjects; ++j)
	{
		if (bestVotes[j])
		{
			outdata.push_back(data[j]);
		}
	}
	paramEstimator->leastSquaresEstimate(outdata, parameters);
	delete []arr;
	delete []bestVotes;
	delete []curVotes;

	return (double)outdata.size()/(double)numDataObjects;
}

template<class T, class S>
void Ransac<T, S>::computeAllChoices(ParameterEsitmator<T, S>* paramEstimator, std::vector<T> &data, 
	int numForEstimate, short* bestVotes, short* curVotes, int &numVotesForBest, int startIndex, 
	int n, int k, int arrIndex, int* arr)
{
	if (k == 0)
	{
		estimate(paramEstimator, data, numForEstimate, bestVotes, curVotes, numVotesForBest, arr);
		return;
	}
	int endIndex = n-k;
	for (int i = startIndex; i <= endIndex; ++i)
	{
		arr[arrIndex] = i;
		computeAllChoices(paramEstimator, data, numForEstimate, bestVotes, curVotes, numVotesForBest,
			i+1, n, k-1, arrIndex+1, arr);
	}
}

template<class T, class S>
void Ransac<T, S>::estimate(ParameterEsitmator<T, S>* paramEstimator, std::vector<T> &data, 
	int numForEstimate, short* bestVotes, short* curVotes, int &numVotesForBest, int* arr)
{
	std::vector<T> exactEstimateData;
	std::vector<S> exactEstimateParameters;
	int numDataObjects;
	int numVotesForCur;
	int j;

	numDataObjects = data.size();
	memset(curVotes, '\0', numDataObjects*sizeof(short));
	numVotesForCur = 0;
	for (j = 0; j < numForEstimate; ++j)
	{
		exactEstimateData.push_back(data[arr[j]]);
	}
	paramEstimator->estimate(exactEstimateData, exactEstimateParameters);
	for (j = 0; j < numDataObjects; ++j)
	{
		if (paramEstimator->agree(exactEstimateParameters, data[j]))
		{
			curVotes[j] = 1;
			++numVotesForCur;
		}
	}
	if (numVotesForCur > numVotesForBest)
	{
		numVotesForBest = numVotesForCur;
		memcpy(bestVotes, curVotes, numDataObjects*sizeof(short));
	}
}

