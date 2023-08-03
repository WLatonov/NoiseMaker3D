#include "CPPNoiseGenerator.h"
#include "MeshProcessor.h"
#include <Python.h>
#include <pybind11/pybind11.h>

void CPPNoiseGenerator::GenerateGaussianNoise(std::string meshDir, double noise_level, int noise_direction) {
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_noiseLevel = noise_level;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "Gaussian");
}

void CPPNoiseGenerator::GenerateImpulseNoise(std::string meshDir, double noise_level, double impulsive_level, int noise_direction)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_noiseLevel = noise_level;
	t_noiseParams.m_impulsiveLevel = impulsive_level;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "Impulsive");
}

void CPPNoiseGenerator::GenerateExponentialNoise(std::string meshDir, double lambda, int noise_direction)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_lambda = lambda;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "Exponential");
}

void CPPNoiseGenerator::GenerateExtremeValueNoise(std::string meshDir, double a, double b, int noise_direction)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_a = a;
	t_noiseParams.m_b = b;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "ExtremeValue");
}

void CPPNoiseGenerator::GenerateGammaNoise(std::string meshDir, double alpha, double beta, int noise_direction)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_alpha = alpha;
	t_noiseParams.m_beta = beta;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "Gamma");
}

void CPPNoiseGenerator::GenerateLaplaceNoise(std::string meshDir, double mu, double b, int noise_direction)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_mu = mu;
	t_noiseParams.m_b = b;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "Laplace");
}

void CPPNoiseGenerator::GenerateLogNormalNoise(std::string meshDir, double m, double s, int noise_direction)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_m = m;
	t_noiseParams.m_s = s;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "LogNormal");
}

void CPPNoiseGenerator::GenerateUniformNoise(std::string meshDir, double a, double b, int noise_direction)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_a = a;
	t_noiseParams.m_b = b;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "Uniform");
}

void CPPNoiseGenerator::GenerateWeibullNoise(std::string meshDir, double a, double b, int noise_direction)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	NoiseParams t_noiseParams;
	t_noiseParams.m_a = a;
	t_noiseParams.m_b = b;
	t_meshProcessor.slotGenNoise(t_noiseParams, noise_direction, "Weibull");
}




void CPPNoiseGenerator::GenerateRandomVerticesTopologyNoise(std::string meshDir, std::string meshNoisedDir, double percentage)
{
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	t_meshProcessor.SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_percentage = percentage;
	t_meshProcessor.slotGenTopologyNoise(t_noiseParams, "RandomVertices");
}

void CPPNoiseGenerator::GenerateRandomFacesTopologyNoise(std::string meshDir, std::string meshNoisedDir, double percentage) {
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	t_meshProcessor.SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_percentage = percentage;
	t_meshProcessor.slotGenTopologyNoise(t_noiseParams, "RandomFaces");
}

void CPPNoiseGenerator::GenerateOneRandomClusterTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider) {
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	t_meshProcessor.SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_standadDeviation = standadDeviation;
	t_noiseParams.m_maxDistance = maxDistance;
	t_noiseParams.m_divider = divider;
	t_meshProcessor.slotGenTopologyNoise(t_noiseParams, "OneRandomCluster");
}

void CPPNoiseGenerator::GenerateOneSpecifiedClusterTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int clusterCenterVertexId) {
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	t_meshProcessor.SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_standadDeviation = standadDeviation;
	t_noiseParams.m_maxDistance = maxDistance;
	t_noiseParams.m_divider = divider;
	t_noiseParams.m_clusterCenterVertexId = clusterCenterVertexId;
	t_meshProcessor.slotGenTopologyNoise(t_noiseParams, "OneClusterInSpecifiedVertex");
}

void CPPNoiseGenerator::SetOfRandomClustersTopologyNoise(std::string meshDir, std::string meshNoisedDir, double standadDeviation, int maxDistance, double divider, int clustersNumber, int interClusterDistance, int makeClustersFar) {
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	t_meshProcessor.SetNoiseFileName(meshNoisedDir);
	TopoNoiseParams t_noiseParams;
	t_noiseParams.m_standadDeviation = standadDeviation;
	t_noiseParams.m_maxDistance = maxDistance;
	t_noiseParams.m_divider = divider;
	t_noiseParams.m_clustersNumber = clustersNumber;
	t_noiseParams.m_interClusterDistance = interClusterDistance;
	t_noiseParams.m_makeClustersFar = (makeClustersFar == 0);
	t_meshProcessor.slotGenTopologyNoise(t_noiseParams, "SetOfRandomClusters");
}

void CPPNoiseGenerator::SetPatchTopologyNoise(std::string meshDir, std::string meshNoisedDir, double c1, double c2, double c3, double c4, double c5, int type, double percentage, int ringsNumber, int ringsNumberToDelete, int regionRadius) {
	MeshProcessor t_meshProcessor = MeshProcessor(meshDir);
	t_meshProcessor.SetNoiseFileName(meshNoisedDir);
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
	
	t_noiseParams.m_percentage = 30;

	t_noiseParams.m_ringsNumber = 2;
	t_noiseParams.m_ringsNumberToDelete = 4;
	t_noiseParams.m_regionRadius = 16;

	t_meshProcessor.slotGenTopologyNoise(t_noiseParams, "PatchModel");
}