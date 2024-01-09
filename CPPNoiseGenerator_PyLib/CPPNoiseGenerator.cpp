#include "CPPNoiseGenerator.h"
#include "MeshProcessor.h"
#include <Python.h>
#include <pybind11/pybind11.h>

namespace pb = pybind11;

void GenerateGaussianNoise(std::string meshDir, double noiseLevel, int noiseDirection, int seed) {
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_noiseLevel = noiseLevel;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Gaussian");
}

void GenerateImpulseNoise(std::string meshDir, double noiseLevel, double impulsiveLevel, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_noiseLevel = noiseLevel;
	t_noiseParams.m_impulsiveLevel = impulsiveLevel;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Impulsive");
}

void GenerateExponentialNoise(std::string meshDir, double lambda, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_lambda = lambda;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Exponential");
}

void GenerateExtremeValueNoise(std::string meshDir, double a, double b, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_a = a;
	t_noiseParams.m_b = b;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "ExtremeValue");
}

void GenerateGammaNoise(std::string meshDir, double alpha, double beta, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_alpha = alpha;
	t_noiseParams.m_beta = beta;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Gamma");
}

void GenerateLaplaceNoise(std::string meshDir, double mu, double b, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_mu = mu;
	t_noiseParams.m_b = b;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Laplace");
}

void GenerateLogNormalNoise(std::string meshDir, double m, double s, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_m = m;
	t_noiseParams.m_s = s;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "LogNormal");
}

void GenerateUniformNoise(std::string meshDir, double a, double b, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_a = a;
	t_noiseParams.m_b = b;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Uniform");
}

void GenerateWeibullNoise(std::string meshDir, double a, double b, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_a = a;
	t_noiseParams.m_b = b;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Weibull");
}

void GenerateCauchyNoise(std::string meshDir, double x0, double gamma, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_x0 = x0;
	t_noiseParams.m_gamma = gamma;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Cauchy");
}

void GenerateFisherNoise(std::string meshDir, double d1, double d2, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_d1 = d1;
	t_noiseParams.m_d2 = d2;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Fisher");
}

void GenerateStudentNoise(std::string meshDir, double n, double scale, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_n = n;
	t_noiseParams.m_scale = scale;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "Student");
}

void GenerateChiSquaredNoise(std::string meshDir, double n, double scale, int noiseDirection, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	NoiseParams t_noiseParams;
	t_noiseParams.m_n = n;
	t_noiseParams.m_scale = scale;
	t_meshProcessor->slotGenNoise(t_noiseParams, noiseDirection, "ChiSquared");
}



void GenerateRandomVerticesTopologyNoise(std::string meshDir, std::string meshNoisedDir, double percentage, int seed)
{
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	t_meshProcessor->SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_percentage = percentage;
	t_meshProcessor->slotGenTopologyNoise(t_noiseParams, "RandomVertices");
}

void GenerateRandomFacesTopologyNoise(std::string meshDir, std::string meshNoisedDir, double percentage, int seed) {
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	t_meshProcessor->SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_percentage = percentage;
	t_meshProcessor->slotGenTopologyNoise(t_noiseParams, "RandomFaces");
}

void GenerateOneRandomClusterTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int seed) {
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	t_meshProcessor->SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_standadDeviation = standadDeviation;
	t_noiseParams.m_maxDistance = maxDistance;
	t_noiseParams.m_divider = divider;
	t_meshProcessor->slotGenTopologyNoise(t_noiseParams, "OneRandomCluster");
}

void GenerateOneSpecifiedClusterTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int clusterCenterVertexId, int seed) {
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	t_meshProcessor->SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_standadDeviation = standadDeviation;
	t_noiseParams.m_maxDistance = maxDistance;
	t_noiseParams.m_divider = divider;
	t_noiseParams.m_clusterCenterVertexId = clusterCenterVertexId;
	t_meshProcessor->slotGenTopologyNoise(t_noiseParams, "OneClusterInSpecifiedVertex");
}

void GenerateSetOfRandomClustersTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int clustersNumber, int interClusterDistance, int makeClustersFar, int seed) {
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	t_meshProcessor->SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_standadDeviation = standadDeviation;
	t_noiseParams.m_maxDistance = maxDistance;
	t_noiseParams.m_divider = divider;
	t_noiseParams.m_clustersNumber = clustersNumber;
	t_noiseParams.m_interClusterDistance = interClusterDistance;
	t_noiseParams.m_makeClustersFar = (makeClustersFar == 0);
	t_meshProcessor->slotGenTopologyNoise(t_noiseParams, "SetOfRandomClusters");
}

void GenerateSetPatchTopologyNoise(std::string meshDir, std::string meshNoisedDir, double c1, double c2, double c3, double c4, double c5, int type, double percentage, int ringsNumber, int ringsNumberToDelete, int regionRadius, int seed) {
	MeshProcessor* t_meshProcessor = new MeshProcessor(meshDir);
	t_meshProcessor->SetSeed(seed);
	t_meshProcessor->SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;

	std::vector<double> constants;

	constants.push_back(c1);
	constants.push_back(c2);
	constants.push_back(c3);
	constants.push_back(c4);
	constants.push_back(c5);
	t_noiseParams.m_patchSelectionConstants = constants;

	if (type == 1) {
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Face);
	}
	else if (type == 2) {
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Edge);
	}
	else if (type == 3) {
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Corner);
	}
	else if (type == 4) {
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Transition);
	}
	else if (type == 5) {
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Face);
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Edge);
	}
	else if (type == 6) {
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Edge);
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Corner);
	}
	else if (type == 7) {
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Face);
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Corner);
	}
	else if (type == 8) {
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Face);
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Edge);
		t_noiseParams.m_patchDeleteTypes.insert(PatchType::Corner);
	}
	
	t_noiseParams.m_percentage = percentage;

	t_noiseParams.m_ringsNumber = ringsNumber;
	t_noiseParams.m_ringsNumberToDelete = ringsNumberToDelete;
	t_noiseParams.m_regionRadius = regionRadius;

	t_meshProcessor->slotGenTopologyNoise(t_noiseParams, "PatchModel");
}

PYBIND11_MODULE(CPPNoiseGenerator, m) {
	m.doc() = "pybind11 example plugin";

	m.def("GenerateGaussianNoisePy", &GenerateGaussianNoise, "Generates Gaussian noise", pb::arg("meshDir") = "", pb::arg("noiseLevel") = double(0.2), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateImpulsiveNoisePy", &GenerateImpulseNoise, "Generates impulsive noise", pb::arg("meshDir") = "", pb::arg("noiseLevel") = double(0.2), pb::arg("impulsiveLevel") = double(0.2), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateExponentialNoisePy", &GenerateExponentialNoise, "Generates Exponential distributed noise", pb::arg("meshDir") = "", pb::arg("lambda") = double(7.0), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateExtremeValueNoisePy", &GenerateExtremeValueNoise, "Generates ExtremeValue distributed noise", pb::arg("meshDir") = "", pb::arg("a") = double(0.0), pb::arg("b") = double(0.3), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateGammaNoisePy", &GenerateGammaNoise, "Generates Gamma distributed noise", pb::arg("meshDir") = "", pb::arg("alpha") = double(0.1), pb::arg("beta") = double(0.3), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateLaplaceNoisePy", &GenerateLaplaceNoise, "Generates Laplace distributed noise", pb::arg("meshDir") = "", pb::arg("mu") = double(0.0), pb::arg("b") = double(0.3), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateLogNormalNoisePy", &GenerateLogNormalNoise, "Generates LogNormal distributed noise", pb::arg("meshDir") = "", pb::arg("m") = double(-5.0), pb::arg("s") = double(40.0), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateUniformNoisePy", &GenerateUniformNoise, "Generates Unifromly distributed noise", pb::arg("meshDir") = "", pb::arg("a") = double(0.2), pb::arg("b") = double(0.4), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateWeibullNoisePy", &GenerateWeibullNoise, "Generates Weibull distributed noise", pb::arg("meshDir") = "", pb::arg("a") = double(1.0), pb::arg("b") = double(0.2), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateCauchyNoisePy", &GenerateCauchyNoise, "Generates Cauchy distributed noise", pb::arg("meshDir") = "", pb::arg("x0") = double(0.0), pb::arg("gamma") = double(0.5), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateFisherNoisePy", &GenerateFisherNoise, "Generates Fisher distributed noise", pb::arg("meshDir") = "", pb::arg("d1") = double(1.0), pb::arg("d2") = double(1.0), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateStudentNoisePy", &GenerateStudentNoise, "Generates Student distributed noise", pb::arg("meshDir") = "", pb::arg("n") = double(1.0), pb::arg("scale") = double(1.0), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateChiSquaredNoisePy", &GenerateChiSquaredNoise, "Generates ChiSquared distributed noise", pb::arg("meshDir") = "", pb::arg("n") = double(1.0), pb::arg("scale") = double(1.0), pb::arg("noiseDirection") = int(0), pb::arg("seed") = int(0));

	m.def("GenerateRandomVerticesTopologyNoisePy", &GenerateRandomVerticesTopologyNoise, "Removes random verteces with adjacent faces", pb::arg("meshDir") = "", pb::arg("meshNoisedDir") = "", pb::arg("percentage") = double(5.0), pb::arg("seed") = int(0));
	m.def("GenerateRandomFacesTopologyNoisePy", &GenerateRandomFacesTopologyNoise, "Removes random faces", pb::arg("meshDir") = "", pb::arg("meshNoisedDir") = "", pb::arg("percentage") = double(5.0), pb::arg("seed") = int(0));
	m.def("GenerateOneRandomClusterTopologyNoisePy", &GenerateOneRandomClusterTopologyNoise, "Removes one normally distributed cluster with randomly selected center", pb::arg("meshDir") = "", pb::arg("meshNoisedDir") = "", pb::arg("standadDeviation") = double(0.5), pb::arg("maxDistance") = int(50), pb::arg("divider") = double(40.0), pb::arg("seed") = int(0));
	m.def("GenerateOneSpecifiedClusterTopologyNoisePy", &GenerateOneSpecifiedClusterTopologyNoise, "Removes one normally distributed cluster with specified center", pb::arg("meshDir") = "", pb::arg("meshNoisedDir") = "", pb::arg("standadDeviation") = double(0.5), pb::arg("maxDistance") = int(50), pb::arg("divider") = double(40.0), pb::arg("clusterCenterVertexId") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateSetOfRandomClustersTopologyNoisePy", &GenerateSetOfRandomClustersTopologyNoise, "Removes set of normally distributed clusters", pb::arg("meshDir") = "", pb::arg("meshNoisedDir") = "", pb::arg("standadDeviation") = double(0.5), pb::arg("maxDistance") = int(50), pb::arg("divider") = double(40.0), pb::arg("clustersNumber") = int(3), pb::arg("interClusterDistance") = int(100), pb::arg("makeClustersFar") = int(0), pb::arg("seed") = int(0));
	m.def("GenerateSetPatchTopologyNoisePy", &GenerateSetPatchTopologyNoise, "Removes specified patches selected by restricting constants", pb::arg("meshDir") = "", pb::arg("meshNoisedDir") = "", pb::arg("c1") = double(0.005), pb::arg("c2") = double(0.01), pb::arg("c3") = double(0.01), pb::arg("c4") = double(0.01), pb::arg("c5") = double(0.05), pb::arg("type") = int(3), pb::arg("percentage") = double(30.0), pb::arg("ringsNumber") = int(2), pb::arg("ringsNumberToDelete") = int(4), pb::arg("regionRadius") = int(16), pb::arg("seed") = int(0));
}