#ifndef NOISE_H
#define NOISE_H

#include "Mesh.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <random>

std::vector<int> GetRandomBoolVector(unsigned seed, size_t length, double percentage);
int GetRandomVertexId(unsigned seed, size_t length);
double GetNormalDistributedProbability(double x, double stDiv);
bool GetBoolWithSpecifiedProbabiliy(unsigned seed, double probability);
std::vector<bool> GetBoolVectorWithSpecifiedProbabiliy(unsigned seed, int length, double probability);

class Noise
{
public:
	enum NoiseType { kGaussian, kImpulsive, kLaplace };
	enum NoiseDirection { kNormal, kRandom };

	explicit Noise() {};
	~Noise() {}

public:
	double generateRandomGaussian(double mean, double StandardDerivation);
	TriMesh::Normal generateRandomDirection();

	void randomGaussianNumbers(unsigned seed, double mean, double StandardDerivation, int number, std::vector<double> &RandomNumbers);
	void randomImpulsiveNumbers(unsigned seed, int min, int max, int number,
		double mean, double StandardDerivation,
		std::vector< std::pair<int, double> > &VertexListAndRandomNumbers);
	void randomDirections(unsigned seed, int number, std::vector<TriMesh::Normal> &RandomDirections);
	void randomLaplaceNumbers(unsigned seed, double mu, double b, int number, std::vector<double>& RandomNumbers);
	void randomExponentialNumbers(unsigned seed, double lambda, int number, std::vector<double>& RandomNumbers);
	void randomExtremeValueNumbers(unsigned seed, double a, double b, int number, std::vector<double>& RandomNumbers);
	void randomUniformNumbers(unsigned seed, double a, double b, int number, std::vector<double>& RandomNumbers);
	void randomWeibullNumbers(unsigned seed, double a, double b, int number, std::vector<double>& RandomNumbers);
	void randomGammaNumbers(unsigned seed, double alpha, double beta, int number, std::vector<double>& RandomNumbers);
	void randomLogNormalNumbers(unsigned seed, double m, double s, int number, std::vector<double>& RandomNumbers);

private:

	double noise_level, impulsive_level;
	int noise_type_index, noise_direction_index;
	double m_noise_level;
	int m_noise_type;
};

#endif // NOISE_H