#pragma once
#include <string>

#ifdef CALCULATIONDLL_EXPORTS  
#define CALCULATION_API __declspec(dllexport)   
#else  
#define CALCULATION_API __declspec(dllimport)   
#endif  
class  CALCULATION_API CPPNoiseGenerator
{
public:
    void GenerateGaussianNoise(std::string meshDir, double noise_level, int noise_direction = 0);
    void GenerateImpulseNoise(std::string meshDir, double noise_level, double impulsive_level, int noise_direction = 0);
    
	void GenerateExponentialNoise(std::string meshDir, double lambda, int noise_direction = 0);
	void GenerateExtremeValueNoise(std::string meshDir, double a, double b, int noise_direction = 0);
	void GenerateGammaNoise(std::string meshDir, double alpha, double beta, int noise_direction = 0);
	void GenerateLaplaceNoise(std::string meshDir, double mu, double b, int noise_direction = 0);
	void GenerateLogNormalNoise(std::string meshDir, double m, double s, int noise_direction = 0);
	void GenerateUniformNoise(std::string meshDir, double a, double b, int noise_direction = 0);
	void GenerateWeibullNoise(std::string meshDir, double a, double b, int noise_direction = 0);

	void GenerateRandomVerticesTopologyNoise(std::string meshDir, std::string meshNoisedDir, double percentage = 5);
	void GenerateRandomFacesTopologyNoise(std::string meshDir, std::string meshNoisedDir, double percentage = 5);
	void GenerateOneRandomClusterTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation = 0.5, int maxDistance = 50, double divider = 40.0);
	void GenerateOneSpecifiedClusterTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation = 0.5, int maxDistance = 50, double divider = 40.0, int clusterCenterVertexId = 0);
	void SetOfRandomClustersTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation = 0.5, int maxDistance = 50, double divider = 40.0, int clustersNumber = 3, int interClusterDistance = 100, int makeClustersFar = 0);
	void SetPatchTopologyNoise(std::string meshDir, std::string meshNoisedDir, double c1 = 0.005, double c2 = 0.01, double c3 = 0.01, double c4 = 0.01, double c5 = 0.05, int type = 3, double percentage = 30, int ringsNumber = 2, int ringsNumberToDelete = 4, int regionRadius = 16);

};
