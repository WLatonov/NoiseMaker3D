#include "PatchData.h"
#include <set>

PatchData::PatchData(TriMesh &mesh, std::vector<std::vector<TriMesh::FaceHandle>> &all_face_neighbor,
	std::vector<double> &face_area, std::vector<TriMesh::Point> &face_centroid,
	TriMesh::FaceIter &iter_face, int num_ring, int radius)
{
	std::vector<std::vector<TriMesh::FaceHandle>> face_neighbor_ring(num_ring + 1);
	std::vector<std::vector<bool>> flag(num_ring + 1, std::vector<bool>(mesh.n_faces(), false));

	face_neighbor_ring[0].push_back(*iter_face);
	if (num_ring >= 1)
	{
		face_neighbor_ring[1] = all_face_neighbor[iter_face->idx()];
	}
	flag[0][iter_face->idx()] = true;

	if (num_ring >= 1)
	{
		for (int i = 0; i < (int)face_neighbor_ring[1].size(); i++)
		{
			flag[1][face_neighbor_ring[1][i].idx()] = true;
		}
	}

	// Collect sequentially faces to rings
	for (int ring = 1; ring < num_ring; ring++)
	{
		for (int i = 0; i < (int)face_neighbor_ring[ring].size(); i++)
		{
			std::vector<TriMesh::FaceHandle> temp_neighbor = all_face_neighbor[face_neighbor_ring[ring][i].idx()];
			for (int t = 0; t < (int)temp_neighbor.size(); t++)
			{
				if ((!flag[ring - 1][temp_neighbor[t].idx()]) && (!flag[ring][temp_neighbor[t].idx()]) && (!flag[ring + 1][temp_neighbor[t].idx()]))
				{
					face_neighbor_ring[ring + 1].push_back(temp_neighbor[t]);
					flag[ring + 1][temp_neighbor[t].idx()] = true;
				}
			}
		}
	}

	// for center point position
	TriMesh::Point centroidc = face_centroid[iter_face->idx()];
	double areac = face_area[iter_face->idx()];
	glm::dvec4 center_face_centorid(areac, centroidc[0], centroidc[1], centroidc[2]);
	m_patchFacesCenters.push_back(center_face_centorid);

	glm::dvec3 normalc(mesh.normal(*iter_face)[0], mesh.normal(*iter_face)[1], mesh.normal(*iter_face)[2]);
	m_patchFacesNormals.push_back(normalc);

	int numFaces = 1;
	for (int ring = 1; ring <= num_ring; ring++)
	{
		for (int i = 0; i < (int)face_neighbor_ring[ring].size(); i++)
		{
			TriMesh::Point centroid = face_centroid[face_neighbor_ring[ring][i].idx()];
			glm::dvec4 neighbor_faces_centorid(areac, centroid[0], centroid[1], centroid[2]);
			m_patchFacesCenters.push_back(neighbor_faces_centorid);
			glm::dvec3 normal(mesh.normal(face_neighbor_ring[ring][i])[0], mesh.normal(face_neighbor_ring[ring][i])[1], mesh.normal(face_neighbor_ring[ring][i])[2]);
			m_patchFacesNormals.push_back(normal);
			numFaces++;
		}
	}
	m_patchFacesNumber = numFaces;
}

PatchData::~PatchData()
{
	m_patchFacesCenters.clear();
	m_patchFacesNormals.clear();
}

void PatchData::SortAndNormalizeEigen(Eigen::EigenSolver<Eigen::Matrix3d>& eigen_solver) {
	double max_eigen_value = eigen_solver.eigenvalues()[0].real();
	double min_eigen_value = eigen_solver.eigenvalues()[0].real();
	int max_ev_idx = 0;
	int min_ev_idx = 0;
	for (int i = 0; i < 3; i++)
	{
		if (max_eigen_value < eigen_solver.eigenvalues()[i].real())
		{
			max_eigen_value = eigen_solver.eigenvalues()[i].real();
			max_ev_idx = i;
		}
		if (min_eigen_value > eigen_solver.eigenvalues()[i].real())
		{
			min_eigen_value = eigen_solver.eigenvalues()[i].real();
			min_ev_idx = i;
		}
	}
	int middle_ev_idx = 3 - max_ev_idx - min_ev_idx;

	// only one face
	if (middle_ev_idx < 0 || middle_ev_idx > 2)
	{
		min_ev_idx = 0;
		middle_ev_idx = 1;
		max_ev_idx = 2;
	}

	glm::dvec3 sorted_eigen_value(eigen_solver.eigenvalues()[max_ev_idx].real(),
		eigen_solver.eigenvalues()[middle_ev_idx].real(),
		eigen_solver.eigenvalues()[min_ev_idx].real());

	glm::dvec3 max_eigen_vec(eigen_solver.eigenvectors().col(max_ev_idx).row(0).value().real(),
		eigen_solver.eigenvectors().col(max_ev_idx).row(1).value().real(),
		eigen_solver.eigenvectors().col(max_ev_idx).row(2).value().real());
	glm::dvec3 middel_eigen_vec(eigen_solver.eigenvectors().col(middle_ev_idx).row(0).value().real(),
		eigen_solver.eigenvectors().col(middle_ev_idx).row(1).value().real(),
		eigen_solver.eigenvectors().col(middle_ev_idx).row(2).value().real());
	glm::dvec3 min_eigen_vec(eigen_solver.eigenvectors().col(min_ev_idx).row(0).value().real(),
		eigen_solver.eigenvectors().col(min_ev_idx).row(1).value().real(),
		eigen_solver.eigenvectors().col(min_ev_idx).row(2).value().real());

	glm::dmat3x3 sorted_eigen_vec(max_eigen_vec, middel_eigen_vec, min_eigen_vec);

	if (sorted_eigen_vec[0][0] * m_patchFacesNormals[0][0]
		+ sorted_eigen_vec[0][1] * m_patchFacesNormals[0][1]
		+ sorted_eigen_vec[0][2] * m_patchFacesNormals[0][2] < 0.)
	{
		sorted_eigen_vec[0] *= -1.;
		sorted_eigen_vec[1] *= -1.;
		sorted_eigen_vec[2] *= -1.;
	}

	m_sortedEigenValues = sorted_eigen_value / sqrt(pow(sorted_eigen_value[0], 2) + pow(sorted_eigen_value[1], 2) + pow(sorted_eigen_value[2], 2));

	// for normal base change
	m_sortedEigenVectors = glm::inverse(sorted_eigen_vec);
}

void PatchData::CalculateEigenValuesAndVectorsForVotingTensor()
{
	double area_max = 0;
	for (int i = 0; i < m_patchFacesNumber; i++)
	{
		if (area_max < m_patchFacesCenters[i][0])
		{
			area_max = m_patchFacesCenters[i][0];
		}
	}

	Eigen::Matrix3d tensor_feature_mat;
	Eigen::Matrix3d temp_mat;
	tensor_feature_mat = Eigen::Matrix3d::Zero();

	for (int i = 1; i < m_patchFacesNumber; i++)
	{
		temp_mat = Eigen::Matrix3d::Zero();
		glm::dvec3 temp_center(m_patchFacesCenters[i][1] - m_patchFacesCenters[0][1],
			m_patchFacesCenters[i][2] - m_patchFacesCenters[0][2],
			m_patchFacesCenters[i][3] - m_patchFacesCenters[0][3]);
		double temp_area = m_patchFacesCenters[i][0] / area_max;
		double  temp_center_norm_para = exp(-3 * glm::length(temp_center));
		glm::dvec3 temp_weight = glm::cross(glm::cross(temp_center, m_patchFacesNormals[i]), temp_center);
		temp_weight = glm::normalize(temp_weight);
		glm::dvec3 temp_normal = 2. * glm::dot(m_patchFacesNormals[i], temp_weight) * temp_weight - m_patchFacesNormals[i];

		for (int r = 0; r < 3; r++)
		{
			for (int c = 0; c < 3; c++)
			{
				temp_mat(r, c) = temp_normal[r] * temp_normal[c];
			}
		}
		temp_mat *= temp_area * temp_center_norm_para;
		tensor_feature_mat += temp_mat;
	}

	Eigen::EigenSolver<Eigen::Matrix3d> eigen_solver(tensor_feature_mat);
	SortAndNormalizeEigen(eigen_solver);
}

double PatchData::m_face1;
double PatchData::m_face2;
double PatchData::m_edge1;
double PatchData::m_edge2;
double PatchData::m_corner1;

PatchType PatchData::GetPatchType() {
	if (m_sortedEigenValues[1] < m_face1 && m_sortedEigenValues[2] < m_face2) {
		return PatchType::Face;
	}
	else if (m_sortedEigenValues[1] > m_edge1 && m_sortedEigenValues[2] < m_edge2) {
		return PatchType::Edge;
	}
	else if (m_sortedEigenValues[2] > m_corner1) {
		return PatchType::Corner;
	}
	else {
		return PatchType::Transition;
	}
}
bool PatchData::IsFace() {
	return m_sortedEigenValues[1] < m_face1 && m_sortedEigenValues[2] < m_face2;
}

bool PatchData::IsEdge() {
	return m_sortedEigenValues[1] > m_edge1 && m_sortedEigenValues[2] < m_edge2;
}

bool PatchData::IsCorner() {
	return m_sortedEigenValues[2] > m_corner1;
}

bool PatchData::CheckIfOneOfTypes(std::set<PatchType>& types) {
	bool isRequired = false;
	for (auto& type : types) {
		if (type == PatchType::Face) {
			isRequired |= IsFace();
		}
		if (type == PatchType::Edge) {
			isRequired |= IsEdge();
		}
		if (type == PatchType::Corner) {
			isRequired |= IsCorner();
		}
		if (type == PatchType::Transition) {
			isRequired |= (!IsFace() && !IsEdge() && !IsCorner());
		}
	}
	return isRequired;
}

bool PatchData::SetEigenTresholdConstants(std::vector<double>& constants) {
	if (constants.size() < 5) {
		return false;
	}
	m_face1 = constants[0];
	m_face2 = constants[1];
	m_edge1 = constants[2];
	m_edge2 = constants[3];
	m_corner1 = constants[4];
	return true;
}