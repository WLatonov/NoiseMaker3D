#pragma once
#include "PatchData.h"
#include "DataManager.h"
#include "MeshDenoisingBase.h"
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <fstream>

struct OnePatchData {
	unsigned char* temp_graph_matrix = nullptr;
	double* patch_node_features = nullptr;
	glm::dvec3 gt_norm;
	glm::dvec3 center_norm;
};

struct OnePatchDataVec {
	std::vector<std::vector<unsigned char>> temp_graph_matrix;
	std::vector<std::vector<double>> patch_node_features;
	glm::dvec3 gt_norm;
	glm::dvec3 center_norm;
};

typedef std::vector<OnePatchData> PatchesData;
//typedef std::vector<OnePatchDataVec> PatchesDataVector;

struct PatchesDataVector {
	std::vector<OnePatchDataVec> m_patchesVector;
	int m_length = 0;

	inline int GetPatchesNumber() {
		return m_length;
	}
	inline void SetPatchesNumber(int length) {
		m_length = length;
	}

	inline OnePatchDataVec GetPatchByIndex(int index) {
		return m_patchesVector[index];
	}
};

class MeshViewer
{
public:
	MeshViewer();
	MeshViewer(std::string originalMeshDir, std::string noisedMeshDir, int var);
	MeshViewer(std::string currentFileName, std::string modelName, bool loadNoisedMesh = false);
	MeshViewer(std::string currentFileName);

private:
	int m_xRot;
	int m_yRot;
	int m_zRot;
	int m_projMatrixLoc;	// vertices' projection matrix location
	int m_modelMatrixLoc;	// vertices' world matrix location
	int m_viewMatrixLoc;	// vertices' camera matrix location
	int m_aLightPosLoc;    // frag's A light position location
	int m_bLightPosLoc;		// frag's B light position location
	int  m_viewPosLoc;	// frag's view position location

private:
	std::string m_current_file_name;
	std::string m_model_name;
	double m_current_noise_level;
	int m_current_noise_type;

	bool m_is_reload = false;

	bool m_is_noise = true;
	bool m_is_gt = false;
	bool m_is_denoised = false;

	bool m_is_have_noise = true;
	bool m_is_have_gt = true;

	OpenMesh::Vec3d m_center;
	double m_max = 0.;

	DataManager* m_data_manager = NULL;
	TriMesh m_noised_tri_mesh;
	TriMesh m_gt_tri_mesh;
	int m_num_faces = 0;

	std::vector<std::vector<TriMesh::FaceHandle>> m_all_face_neighbor;
	std::vector<double> m_face_area;
	std::vector<TriMesh::Point> m_face_centroid;

	std::vector<glm::vec3> m_vertices_kd_tree;
	shen::Geometry::EasyFlann m_flann_kd_tree;

	void meshInitializedGet(const char* file_name, bool is_original);

	// same with mesh denoising base
	enum FaceNeighborType { kVertexBased, kEdgeBased, kRadiusBased };

	void getFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, FaceNeighborType face_neighbor_type, std::vector<TriMesh::FaceHandle> &face_neighbor);
	void getAllFaceNeighbor(TriMesh &mesh, std::vector< std::vector<TriMesh::FaceHandle> > &all_face_neighbor, FaceNeighborType face_neighbor_type = kVertexBased, bool include_central_face = false);
	void getFaceArea(TriMesh &mesh, std::vector<double> &area);
	void getFaceCentroid(TriMesh &mesh, std::vector<TriMesh::Point> &centroid);

public:
	void slotLoadNoise(std::string noise_model_file);
	void slotLoadGT(std::string gt_model_file);
	void slotDelete();
	void slotGenNoise(double noise_level, std::string noise_type);
};