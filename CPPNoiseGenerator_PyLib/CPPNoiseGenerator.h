#pragma once
#include <string>

#ifdef CALCULATIONDLL_EXPORTS  
#define CALCULATION_API __declspec(dllexport)   
#else  
#define CALCULATION_API __declspec(dllimport)   
#endif

extern "C" CALCULATION_API void GenerateGaussianNoise(std::string meshDir, double noiseLevel, int noiseDirection = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateImpulseNoise(std::string meshDir, double noiseLevel, double impulsiveLevel, int noiseDirection = 0, int seed = 0);

extern "C" CALCULATION_API void GenerateExponentialNoise(std::string meshDir, double lambda, int noiseDirection = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateExtremeValueNoise(std::string meshDir, double a, double b, int noiseDirection = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateGammaNoise(std::string meshDir, double alpha, double beta, int noiseDirection = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateLaplaceNoise(std::string meshDir, double mu, double b, int noiseDirection = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateLogNormalNoise(std::string meshDir, double m, double s, int noiseDirection = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateUniformNoise(std::string meshDir, double a, double b, int noiseDirection = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateWeibullNoise(std::string meshDir, double a, double b, int noiseDirection = 0, int seed = 0);

extern "C" CALCULATION_API void GenerateRandomVerticesTopologyNoise(std::string meshDir, std::string meshNoisedDir, double percentage = 5, int seed = 0);
extern "C" CALCULATION_API void GenerateRandomFacesTopologyNoise(std::string meshDir, std::string meshNoisedDir, double percentage = 5, int seed = 0);
extern "C" CALCULATION_API void GenerateOneRandomClusterTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation = 0.5, int maxDistance = 50, double divider = 40.0, int seed = 0);
extern "C" CALCULATION_API void GenerateOneSpecifiedClusterTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation = 0.5, int maxDistance = 50, double divider = 40.0, int clusterCenterVertexId = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateSetOfRandomClustersTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation = 0.5, int maxDistance = 50, double divider = 40.0, int clustersNumber = 3, int interClusterDistance = 100, int makeClustersFar = 0, int seed = 0);
extern "C" CALCULATION_API void GenerateSetPatchTopologyNoise(std::string meshDir, std::string meshNoisedDir, double c1 = 0.005, double c2 = 0.01, double c3 = 0.01, double c4 = 0.01, double c5 = 0.05, int type = 3, double percentage = 30, int ringsNumber = 2, int ringsNumberToDelete = 4, int regionRadius = 16, int seed = 0);

