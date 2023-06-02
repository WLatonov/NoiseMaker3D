#include "GCNMeshDenoiser19.h"
#include "MeshViewer.h"
#include <Python.h>
#include <pybind11/pybind11.h>

void GCNMeshDenoiser19::GenerateGaussianNoise(std::string meshDir, double noise_level, int noise_direction) {
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_noiseLevel = noise_level;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "Gaussian");
}

void GCNMeshDenoiser19::GenerateImpulseNoise(std::string meshDir, double noise_level, double impulsive_level, int noise_direction)
{
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_noiseLevel = noise_level;
	m_noiseParams.m_impulsiveLevel = impulsive_level;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "Impulsive");
}

void GCNMeshDenoiser19::GenerateExponentialNoise(std::string meshDir, double lambda, int noise_direction)
{
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_lambda = lambda;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "Exponential");
}

void GCNMeshDenoiser19::GenerateExtremeValueNoise(std::string meshDir, double a, double b, int noise_direction)
{
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_a = a;
	m_noiseParams.m_b = b;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "ExtremeValue");
}

void GCNMeshDenoiser19::GenerateGammaNoise(std::string meshDir, double alpha, double beta, int noise_direction)
{
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_alpha = alpha;
	m_noiseParams.m_beta = beta;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "Gamma");
}

void GCNMeshDenoiser19::GenerateLaplaceNoise(std::string meshDir, double mu, double b, int noise_direction)
{
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_mu = mu;
	m_noiseParams.m_b = b;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "Laplace");
}

void GCNMeshDenoiser19::GenerateLogNormalNoise(std::string meshDir, double m, double s, int noise_direction)
{
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_m = m;
	m_noiseParams.m_s = s;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "LogNormal");
}

void GCNMeshDenoiser19::GenerateUniformNoise(std::string meshDir, double a, double b, int noise_direction)
{
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_a = a;
	m_noiseParams.m_b = b;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "Uniform");
}

void GCNMeshDenoiser19::GenerateWeibullNoise(std::string meshDir, double a, double b, int noise_direction)
{
	MeshViewer t_meshViewer = MeshViewer(meshDir);
	NoiseParams m_noiseParams;
	m_noiseParams.m_a = a;
	m_noiseParams.m_b = b;
	t_meshViewer.slotGenNoise(m_noiseParams, noise_direction, "Weibull");
}