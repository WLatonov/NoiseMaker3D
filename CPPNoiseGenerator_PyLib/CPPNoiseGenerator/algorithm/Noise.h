#ifndef NOISE_H
#define NOISE_H

#include "Mesh.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <random>

std::vector<int> GetRandomBoolVector(size_t length, double percentage);
int GetRandomVertexId(size_t length);
double GetNormalDistributedProbability(double x, double stDiv);
bool GetBoolWithSpecifiedProbabiliy(double probability);
std::vector<bool> GetBoolVectorWithSpecifiedProbabiliy(int length, double probability);

class Noise
{
public:
	enum NoiseType { kGaussian, kImpulsive, kLaplace };
	enum NoiseDirection { kNormal, kRandom };

	explicit Noise() {};
	~Noise() {}

public:
	void initParameters(int noise_type_, double noise_level_);

	double generateRandomGaussian(double mean, double StandardDerivation);

	TriMesh::Normal generateRandomDirection();

	void randomGaussianNumbers(double mean, double StandardDerivation, int number, std::vector<double> &RandomNumbers);
	void randomImpulsiveNumbers(int min, int max, int number,
		double mean, double StandardDerivation,
		std::vector< std::pair<int, double> > &VertexListAndRandomNumbers);
	void randomDirections(int number, std::vector<TriMesh::Normal> &RandomDirections);
	void randomLaplaceNumbers(double mu, double b, int number, std::vector<double>& RandomNumbers);
	void randomExponentialNumbers(double lambda, int number, std::vector<double>& RandomNumbers);
	void randomExtremeValueNumbers(double a, double b, int number, std::vector<double>& RandomNumbers);
	void randomUniformNumbers(double a, double b, int number, std::vector<double>& RandomNumbers);
	void randomWeibullNumbers(double a, double b, int number, std::vector<double>& RandomNumbers);
	void randomGammaNumbers(double alpha, double beta, int number, std::vector<double>& RandomNumbers);
	void randomLogNormalNumbers(double m, double s, int number, std::vector<double>& RandomNumbers);

private:

	double noise_level, impulsive_level;
	int noise_type_index, noise_direction_index;
	double m_noise_level;
	int m_noise_type;
};

#endif // NOISE_H