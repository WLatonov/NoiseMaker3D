#pragma once
#include "PatchData.h"
#include "Mesh.h"
#include "Noise.h"
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <fstream>
#include <queue>
#include <string>

// Param structures for different noise types
struct NoiseParams
{
	//Gauss / Impulse
	double m_noiseLevel;

	//Impulse
	double m_impulsiveLevel;

	//Laplace
	double m_mu;
	//Laplace  / Uniform / Weibull
	double m_b;

	//Exponential
	double m_lambda;

	//Extreme / Uniform / Weibull
	double m_a;

	//Gamma
	double m_alpha;
	double m_beta;

	//LogNormal
	double m_m;
	double m_s;
};

// Param structures for different types of topology noise 
struct TopoNoiseParams
{
	//Random
	double m_percentage;

	//Cluster, Normal Destribution
	double m_standadDeviation;
	//Max distance of vertices to delete in cluster noise (number of edges)
	int m_maxDistance;
	//Divider for Normal destribution in cluster noise. 
	//The divider is required to adjust a density of vertices to delete in the vicinity of cluster center vertex
	double m_divider;
	//If one cluster center vertex is specified
	int m_clusterCenterVertexId;
	//The number of clusters to delete if more then one cluster removing is required
	int m_clustersNumber;
	//The distance between clusters center points if many clusters are to delete (number of edges)
	int m_interClusterDistance;
	//1 - make all clusters as close to each other as possible, 0 - as far as possible
	bool m_makeClustersFar;
	
	// Patches types to delete
	std::set<PatchType> m_patchDeleteTypes;
	//Constants for patches selection to delete
	std::vector<double> m_patchSelectionConstants;

	// Rings number in one patch
	int m_ringsNumber;
	// Rings number in one patch to delete
	int m_ringsNumberToDelete;
	// Region radius for one patch
	int m_regionRadius;
};

using CenterSearchQueuePair = std::pair<TriMesh::VertexHandle, int>;

class Comparator
{
	bool m_minFirst;
public:
	Comparator(const bool& minFirst = true)
	{
		m_minFirst = minFirst;
	}
	bool operator() (const CenterSearchQueuePair& left, const CenterSearchQueuePair& right) const
	{
		if (m_minFirst) return (left.second <= right.second);
		else return (left.second > right.second);
	}
};

using CenterSearchQueue = std::priority_queue<CenterSearchQueuePair, std::vector<CenterSearchQueuePair>, Comparator>;

class MeshProcessor
{
public:
	MeshProcessor(std::string currentFileName);

private:
	unsigned m_seed;
	std::default_random_engine m_generator;

	std::string m_gtFile;
	std::string m_noisyFile;

	OpenMesh::Vec3d m_center;
	double m_max = 0.;

	TriMesh* m_mesh;
	int m_facesNumber = 0;

	void meshInitializedGet(const char* file_name, bool is_original);
	void WriteNoisyMesh();
	void getFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, std::vector<TriMesh::FaceHandle> &face_neighbor);
	void getAllFaceNeighbor(TriMesh &mesh, std::vector< std::vector<TriMesh::FaceHandle> > &all_face_neighbor, bool include_central_face = false);
	void getFaceArea(TriMesh &mesh, std::vector<double> &area);
	void getFaceCentroid(TriMesh &mesh, std::vector<TriMesh::Point> &centroid);

public:
	void SetSeed(int seed);
	void SetGenerator(unsigned seed);
	void SetNoiseFileName(NoiseParams& noiseParams, int noiseDirection, std::string noise_type);
	void SetNoiseFileName(std::string noiseFileName);
	std::vector<double> GetNoiseNumbers(std::string& noiseType, double averageLength, NoiseParams& noiseParams, Noise& noise);
	void slotGenNoise(NoiseParams& noiseParams, int noiseDirection, std::string noise_type);

	void MakePatchNoise(TriMesh& mesh, TopoNoiseParams& noiseParams);
	bool FindInVec(int index, std::vector<int>& vertexIndexToDelete);
	TriMesh::VertexHandle FindNextClusterCenter(CenterSearchQueue& possibleNextClusterCenter);
	void DeleteSetOfClusters(TriMesh& mesh, TopoNoiseParams& noiseParams);
	void DeleteOneCluster(TriMesh& mesh, int startVertId, double standardDeviation, double divider, int maxDist);
	void slotGenTopologyNoise(TopoNoiseParams& noiseParams, std::string noise_type);
};