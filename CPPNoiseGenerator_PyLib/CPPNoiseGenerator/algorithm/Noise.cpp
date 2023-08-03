#include "Noise.h"

void Noise::initParameters(int noise_type_, double noise_level_)
{
	noise_level = noise_level_;
	impulsive_level = noise_level_;

	noise_type_index = noise_type_;
	noise_direction_index = 0;
}

// Box-Muller method
double Noise::generateRandomGaussian(double mean, double StandardDerivation)
{
	static double v1, v2, s;
	static int phase = 0;
	double x;

	if (phase == 0)
	{
		do
		{
			v1 = -1 + 2 * (double)rand() / (double)RAND_MAX;
			v2 = -1 + 2 * (double)rand() / (double)RAND_MAX;
			s = v1 * v1 + v2 * v2;
		} while (s >= 1 || s == 0);

		x = v1 * sqrt(-2 * log(s) / s);
	}
	else
		x = v2 * sqrt(-2 * log(s) / s);

	phase = 1 - phase;

	return x * StandardDerivation + mean;
}

TriMesh::Normal Noise::generateRandomDirection()
{
	double x, y, z, length;
	do {
		x = -1 + 2 * (double)rand() / (double)RAND_MAX;
		y = -1 + 2 * (double)rand() / (double)RAND_MAX;
		length = x * x + y * y;
	} while (length > 1);

	const double r = 2 * std::sqrt(1 - length);

	x *= r;
	y *= r;
	z = 1 - 2 * length;

	return TriMesh::Normal(x, y, z);
}

void Noise::randomGaussianNumbers(double mean, double StandardDerivation, int number, std::vector<double> &RandomNumbers)
{
	RandomNumbers.resize(number, 0.0);

	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		RandomNumbers[i] = generateRandomGaussian(mean, StandardDerivation);
	}
}

void Noise::randomImpulsiveNumbers(int min, int max, int number, double mean, double StandardDerivation,
	std::vector<std::pair<int, double> > &VertexListAndRandomNumbers)
{
	int range = max - min + 1;
	if (number > range) return;

	VertexListAndRandomNumbers.resize(number, std::make_pair(0, 0.0));

	std::vector<double> randomNumbers;
	randomGaussianNumbers(mean, StandardDerivation, number, randomNumbers);

	srand((unsigned int)time(NULL));
	std::vector<int> rangeVector(range);
	for (int i = 0; i < range; i++)
		rangeVector[i] = min + i;

	srand((unsigned int)time(NULL));
	std::vector<int> vertexIndexList(number);
	for (int i = 0; i < number; i++) {
		int pos = (int)((double)rand() / RAND_MAX * range);
		vertexIndexList[i] = rangeVector[pos];
		range--;
		std::swap(rangeVector[pos], rangeVector[range]);
	}

	for (int i = 0; i < number; i++)
		VertexListAndRandomNumbers[i] = std::make_pair(vertexIndexList[i], randomNumbers[i]);
}

void Noise::randomDirections(int number, std::vector<TriMesh::Normal> &RandomDirections)
{
	RandomDirections.resize(number, TriMesh::Normal(0.0, 0.0, 0.0));

	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		RandomDirections[i] = generateRandomDirection();
	}
}

void Noise::randomLaplaceNumbers(double mu, double b, int number, std::vector<double>& RandomNumbers)
{
	RandomNumbers.resize(number, 0.0);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(-0.5, 0.5);


	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		double value = distribution(generator);
		double sign = (value >= 0) ? 1.0 : -1.0;
		if (value < 0) {
			value *= -1.0;
		}
		RandomNumbers[i] = mu - b * sign * std::log(1 - 2 * value);
	}
}

void Noise::randomExponentialNumbers(double lambda, int number, std::vector<double>& RandomNumbers)
{
	RandomNumbers.resize(number, 0.0);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::exponential_distribution<double> distribution(lambda);

	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		RandomNumbers[i] = distribution(generator);
	}
}

void Noise::randomExtremeValueNumbers(double a, double b, int number, std::vector<double>& RandomNumbers)
{
	RandomNumbers.resize(number, 0.0);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::extreme_value_distribution<double> distribution(a, b);

	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		RandomNumbers[i] = distribution(generator);
	}
}

void Noise::randomUniformNumbers(double a, double b, int number, std::vector<double>& RandomNumbers)
{
	RandomNumbers.resize(number, 0.0);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(a, b);

	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		RandomNumbers[i] = distribution(generator);
	}
}

void Noise::randomWeibullNumbers(double a, double b, int number, std::vector<double>& RandomNumbers)
{
	RandomNumbers.resize(number, 0.0);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::weibull_distribution<double> distribution(a, b);

	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		RandomNumbers[i] = distribution(generator);
	}
}

void Noise::randomGammaNumbers(double alpha, double beta, int number, std::vector<double>& RandomNumbers)
{
	RandomNumbers.resize(number, 0.0);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::gamma_distribution<double> distribution(alpha, beta);

	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		RandomNumbers[i] = distribution(generator);
	}
}

void Noise::randomLogNormalNumbers(double m, double s, int number, std::vector<double>& RandomNumbers) 
{
	RandomNumbers.resize(number, 0.0);

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::lognormal_distribution<double> distribution(m, s);

	srand((unsigned int)time(NULL));
	for (int i = 0; i < number; i++) {
		RandomNumbers[i] = distribution(generator);
	}
}

std::vector<int> GetRandomBoolVector(size_t length, double percentage)
{
	std::vector<int> result;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::bernoulli_distribution distribution(percentage);

	srand((unsigned int)time(NULL));
	for (size_t i = 0; i < length; ++i) {
		auto randNumb = distribution(generator);
		if (randNumb) {
			result.push_back((int)i);
		}		
	}
	return result;
}

int GetRandomVertexId(size_t length)
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_int_distribution<> distribution(0, length-1);
	return distribution(generator);
}

double GetNormalDistributedProbability(double x, double stDiv) {
	return (1 / (stDiv * sqrt(2 * M_PI))) * exp(-(x / stDiv)* (x / stDiv) / 2);
}

bool GetBoolWithSpecifiedProbabiliy(double probability) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::bernoulli_distribution distribution(probability);
	return distribution(generator);
}

std::vector<bool> GetBoolVectorWithSpecifiedProbabiliy(int length, double probability) {
	std::vector<bool> result;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::bernoulli_distribution distribution(probability);
	for (int k = 0; k < length; ++k) {
		bool currentRes = distribution(generator);
		result.push_back(currentRes);
	}
	return result;
}
