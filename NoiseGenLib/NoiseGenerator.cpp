#include <pybind11/pybind11.h>
#include <Python.h>
#include "CPPNoiseGenerator.h"

// Node noise algorithms

void GenerateGaussianNoisePy(std::string meshDir, double noiseLevel, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateGaussianNoise(meshDir, noiseLevel, noiseDirection);
}

void GenerateImpulsiveNoisePy(std::string meshDir, double noiseLevel, double impulsiveLevel, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateImpulseNoise(meshDir, noiseLevel, impulsiveLevel, noiseDirection);
}

void GenerateExponentialNoisePy(std::string meshDir, double lambda, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateExponentialNoise(meshDir, lambda, noiseDirection);
}

void GenerateExtremeValueNoisePy(std::string meshDir, double a, double b, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateExtremeValueNoise(meshDir, a, b, noiseDirection);
}

void GenerateGammaNoisePy(std::string meshDir, double alpha, double beta, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateGammaNoise(meshDir, alpha, beta, noiseDirection);
}

void GenerateLaplaceNoisePy(std::string meshDir, double mu, double b, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateLaplaceNoise(meshDir, mu, b, noiseDirection);
}

void GenerateLogNormalNoisePy(std::string meshDir, double m, double s, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateLogNormalNoise(meshDir, m, s, noiseDirection);
}

void GenerateUniformNoisePy(std::string meshDir, double a, double b, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateUniformNoise(meshDir, a, b, noiseDirection);
}

void GenerateWeibullNoisePy(std::string meshDir, double a, double b, int noiseDirection) {
    CPPNoiseGenerator obj;
    obj.GenerateWeibullNoise(meshDir, a, b, noiseDirection);
}

// Topology noise algorithms

void GenerateRandomVerticesTopologyNoisePy(std::string meshDir, std::string meshNoisedDir, double percentage) {
    CPPNoiseGenerator obj;
    obj.GenerateRandomVerticesTopologyNoise(meshDir, meshNoisedDir, percentage);
}

void GenerateRandomFacesTopologyNoisePy(std::string meshDir, std::string meshNoisedDir, double percentage) {
    CPPNoiseGenerator obj;
    obj.GenerateRandomFacesTopologyNoise(meshDir, meshNoisedDir, percentage);
}

void GenerateOneRandomClusterTopologyNoisePy(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider) {
    CPPNoiseGenerator obj;
    obj.GenerateOneRandomClusterTopologyNoise(meshDir, meshNoisedDir, standadDeviation, maxDistance, divider);
}

void GenerateOneSpecifiedClusterTopologyNoisePy(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int clusterCenterVertexId) {
    CPPNoiseGenerator obj;
    obj.GenerateOneSpecifiedClusterTopologyNoise(meshDir, meshNoisedDir, standadDeviation, maxDistance, divider, clusterCenterVertexId);
}

void GenerateSetOfRandomClustersTopologyNoisePy(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int clustersNumber, int interClusterDistance, int makeClustersFar) {
    CPPNoiseGenerator obj;
    obj.SetOfRandomClustersTopologyNoise(meshDir, meshNoisedDir, standadDeviation, maxDistance, divider, clustersNumber, interClusterDistance, makeClustersFar);
}

void GenerateSetPatchTopologyNoisePy(std::string meshDir, std::string meshNoisedDir, double c1, double c2, double c3, double c4, double c5, int type, double percentage, int ringsNumber, int ringsNumberToDelete, int regionRadius) {
    CPPNoiseGenerator obj;
    obj.SetPatchTopologyNoise(meshDir, meshNoisedDir, c1, c2, c3, c4, c5, type, percentage, ringsNumber, ringsNumberToDelete, regionRadius);
}




PYBIND11_MODULE(NoiseGenerator, m) {
    m.doc() = "pybind11 example plugin";

    m.def("GenerateGaussianNoisePy", &GenerateGaussianNoisePy, "Generates Gaussian noise");
    m.def("GenerateImpulsiveNoisePy", &GenerateImpulsiveNoisePy, "Generates impulsive noise");
    m.def("GenerateExponentialNoisePy", &GenerateExponentialNoisePy, "Generates Exponential distributed noise");
    m.def("GenerateExtremeValueNoisePy", &GenerateExtremeValueNoisePy, "Generates ExtremeValue distributed noise");
    m.def("GenerateGammaNoisePy", &GenerateGammaNoisePy, "Generates Gamma distributed noise");
    m.def("GenerateLaplaceNoisePy", &GenerateLaplaceNoisePy, "Generates Laplace distributed noise");
    m.def("GenerateLogNormalNoisePy", &GenerateLogNormalNoisePy, "Generates LogNormal distributed noise");
    m.def("GenerateUniformNoisePy", &GenerateUniformNoisePy, "Generates Unifromly distributed noise");
    m.def("GenerateWeibullNoisePy", &GenerateWeibullNoisePy, "Generates Weibull distributed noise");

    m.def("GenerateRandomVerticesTopologyNoisePy", &GenerateRandomVerticesTopologyNoisePy, "Removes random verteces with adjacent faces");
    m.def("GenerateRandomFacesTopologyNoisePy", &GenerateRandomFacesTopologyNoisePy, "Removes random faces");
    m.def("GenerateOneRandomClusterTopologyNoisePy", &GenerateOneRandomClusterTopologyNoisePy, "Removes one normally distributed cluster with randomly selected center");
    m.def("GenerateOneSpecifiedClusterTopologyNoisePy", &GenerateOneSpecifiedClusterTopologyNoisePy, "Removes one normally distributed cluster with specified center");
    m.def("GenerateSetOfRandomClustersTopologyNoisePy", &GenerateSetOfRandomClustersTopologyNoisePy, "Removes set of normally distributed clusters");
    m.def("GenerateSetPatchTopologyNoisePy", &GenerateSetPatchTopologyNoisePy, "Removes specified patches selected by restricting constants");
}






//Работает!

//int main()
//{
//    GCNMeshDenoiser19 obj;
//    //obj.GenerateGaussianNoise("C:\\Users\\20962270\\source\\repos\\GCNMeshDenoiser19\\GCNMeshDenoiser19\\examples\\example.obj", 0.9, 0.9);
//    //obj.GenerateLaplaceNoise("C:\\Users\\20962270\\source\\repos\\GCNMeshDenoiser19\\GCNMeshDenoiser19\\examples\\example.obj", 1, 0.1, 1);
//    //obj.GenerateImpulseNoise("C:\\Users\\20962270\\source\\repos\\GCNMeshDenoiser19\\GCNMeshDenoiser19\\examples\\example.obj", 0.1, 0.2, 0);
//    obj.GenerateGaussianNoise("C:\\Users\\20962270\\source\\repos\\GCNMeshDenoiser19\\GCNMeshDenoiser19\\examples\\example.obj", 0.3, 0);
//}