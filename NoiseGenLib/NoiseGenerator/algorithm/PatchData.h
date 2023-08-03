#pragma once
#include "Mesh.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <vector>
#include <glm.hpp>

enum PatchType : int { 
	Face = 1, 
	Edge = 2, 
	Corner = 3, 
	Transition = 4 
};

class PatchData
{
public:
	static double m_face1;
	static double m_face2;
	static double m_edge1;
	static double m_edge2;
	static double m_corner1;

	PatchData(TriMesh &mesh, std::vector<std::vector<TriMesh::FaceHandle>> &all_face_neighbor,
		std::vector<double> &face_area, std::vector<TriMesh::Point> &face_centroid,
		TriMesh::FaceIter &iter_face, int num_ring = 2, int radius = 16);
	~PatchData();

	void CalculateEigenValuesAndVectorsForVotingTensor();
	PatchType GetPatchType();
	bool IsFace();
	bool IsEdge();
	bool IsCorner();
	bool CheckIfOneOfTypes(std::set<PatchType>& types);
	static bool SetEigenTresholdConstants(std::vector<double>& constants);

private:
	glm::dvec3 m_sortedEigenValues;
	glm::dmat3x3 m_sortedEigenVectors;

	int m_patchFacesNumber = 0;
	std::vector<glm::dvec4> m_patchFacesCenters;
	std::vector<glm::dvec3> m_patchFacesNormals;

	void SortAndNormalizeEigen(Eigen::EigenSolver<Eigen::Matrix3d>& eigen_solver);
};