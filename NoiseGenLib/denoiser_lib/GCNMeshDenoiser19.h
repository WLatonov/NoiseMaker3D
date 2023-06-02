#pragma once
#include <string>

#ifdef CALCULATIONDLL_EXPORTS  
#define CALCULATION_API __declspec(dllexport)   
#else  
#define CALCULATION_API __declspec(dllimport)   
#endif  
class  CALCULATION_API GCNMeshDenoiser19
{
public:
    void GenerateGaussianNoise(std::string meshDir, double noise_level, int noise_direction);
    void GenerateImpulseNoise(std::string meshDir, double noise_level, double impulsive_level, int noise_direction);
    
	void GenerateExponentialNoise(std::string meshDir, double lambda, int noise_direction);
	void GenerateExtremeValueNoise(std::string meshDir, double a, double b, int noise_direction);
	void GenerateGammaNoise(std::string meshDir, double alpha, double beta, int noise_direction);
	void GenerateLaplaceNoise(std::string meshDir, double mu, double b, int noise_direction);
	void GenerateLogNormalNoise(std::string meshDir, double m, double s, int noise_direction);
	void GenerateUniformNoise(std::string meshDir, double a, double b, int noise_direction);
	void GenerateWeibullNoise(std::string meshDir, double a, double b, int noise_direction);
};
