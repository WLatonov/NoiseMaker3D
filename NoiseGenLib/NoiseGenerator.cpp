#include <pybind11/pybind11.h>
#include <Python.h>
#include "GCNMeshDenoiser19.h"

void GenerateGaussianNoisePy(std::string meshDir, double noiseLevel, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateGaussianNoise(meshDir, noiseLevel, noiseDirection);
}

void GenerateImpulsiveNoisePy(std::string meshDir, double noiseLevel, double impulsiveLevel, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateImpulseNoise(meshDir, noiseLevel, impulsiveLevel, noiseDirection);
}

void GenerateExponentialNoisePy(std::string meshDir, double lambda, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateExponentialNoise(meshDir, lambda, noiseDirection);
}

void GenerateExtremeValueNoisePy(std::string meshDir, double a, double b, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateExtremeValueNoise(meshDir, a, b, noiseDirection);
}

void GenerateGammaNoisePy(std::string meshDir, double alpha, double beta, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateGammaNoise(meshDir, alpha, beta, noiseDirection);
}

void GenerateLaplaceNoisePy(std::string meshDir, double mu, double b, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateLaplaceNoise(meshDir, mu, b, noiseDirection);
}

void GenerateLogNormalNoisePy(std::string meshDir, double m, double s, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateLogNormalNoise(meshDir, m, s, noiseDirection);
}

void GenerateUniformNoisePy(std::string meshDir, double a, double b, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateUniformNoise(meshDir, a, b, noiseDirection);
}

void GenerateWeibullNoisePy(std::string meshDir, double a, double b, int noiseDirection) {
    GCNMeshDenoiser19 obj;
    obj.GenerateWeibullNoise(meshDir, a, b, noiseDirection);
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