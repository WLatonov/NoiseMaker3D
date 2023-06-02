#include <memory>

#include "MeshViewer.h"
#include "FlannKDTree.h"
#include "MeshNormalFiltering.h"

#include <omp.h>
#include <cstdio>
#include <ctime>

MeshViewer::MeshViewer(std::string currentFileName) : m_xRot(0), m_yRot(0), m_zRot(0)
{
	m_data_manager = new DataManager();
	m_gtFile = currentFileName;
	meshInitializedGet(currentFileName.c_str(), true);
	m_current_noise_level = -1.0;
	m_current_noise_type = -1;

	m_noisyFile = m_gtFile;
	m_noisyFile = m_noisyFile.substr(0, m_noisyFile.size() - 4);
	m_noisyFile = m_noisyFile + "_noise.obj";
}

void MeshViewer::meshInitializedGet(const char* file_name, bool is_original)
{
	std::cout << "Loading... mesh dir: " << file_name << std::endl;

	if (!is_original)
	{
		if (!OpenMesh::IO::read_mesh(m_noised_tri_mesh, file_name)) {
			std::cout << "Reading error!" << std::endl;
		}
		
		int new_num_faces = m_noised_tri_mesh.n_faces();
		if (m_num_faces == 0)
		{
			m_num_faces = new_num_faces;
		}
		else
		{
			if (new_num_faces != m_num_faces)
			{
				std::cout << "Not Match!" << std::endl;
				return;
			}
		}
		std::cout << "\tNumber of faces: " << m_num_faces << std::endl;

		m_is_have_noise = true;

		m_is_noise = true;
		m_is_gt = false;
		m_is_denoised = false;

		m_data_manager->setNoisyMesh(m_noised_tri_mesh);
		m_data_manager->setDenoisedMesh(m_noised_tri_mesh);

		m_noised_tri_mesh.request_face_normals();
		m_noised_tri_mesh.request_vertex_normals();
		m_noised_tri_mesh.update_normals();

		m_center = OpenMesh::Vec3d(0., 0., 0.);
		for (TriMesh::VertexIter v_it = m_noised_tri_mesh.vertices_begin(); v_it != m_noised_tri_mesh.vertices_end(); v_it++)
		{
			OpenMesh::Vec3d center_position_temp = m_noised_tri_mesh.point(*v_it);
			m_center += center_position_temp;
		}
		m_center /= int(m_noised_tri_mesh.n_vertices());

		m_max = 0;
		for (TriMesh::VertexIter v_it = m_noised_tri_mesh.vertices_begin(); v_it != m_noised_tri_mesh.vertices_end(); v_it++)
		{
			OpenMesh::Vec3d position_temp = m_noised_tri_mesh.point(*v_it) - m_center;
			for (int i = 0; i < 3; i++)
			{
				double temp_max = abs(position_temp[i]);
				if (temp_max > m_max)
				{
					m_max = temp_max;
				}
			}
		}

		for (TriMesh::VertexIter v_it = m_noised_tri_mesh.vertices_begin(); v_it != m_noised_tri_mesh.vertices_end(); v_it++)
		{
			m_noised_tri_mesh.point(*v_it) -= m_center;
			m_noised_tri_mesh.point(*v_it) /= m_max / 1.0;
			OpenMesh::Vec3d center_position_temp = m_noised_tri_mesh.point(*v_it);
			m_vertices_kd_tree.push_back(glm::vec3(center_position_temp[0], center_position_temp[1], center_position_temp[2]));
		}

		m_flann_kd_tree = shen::Geometry::EasyFlann(m_vertices_kd_tree);

		FaceNeighborType face_neighbor_type = kVertexBased;

		getAllFaceNeighbor(m_noised_tri_mesh, m_all_face_neighbor, face_neighbor_type, false);
		getFaceArea(m_noised_tri_mesh, m_face_area);
		getFaceCentroid(m_noised_tri_mesh, m_face_centroid);
	}
	else
	{
		OpenMesh::IO::read_mesh(m_gt_tri_mesh, file_name);

		int new_num_faces = m_gt_tri_mesh.n_faces();
		std::cout << "Faces = " << new_num_faces << std::endl;
		if (m_num_faces == 0)
		{
			m_num_faces = new_num_faces;
		}
		else
		{
			if (new_num_faces != m_num_faces)
			{
				std::cout << "Not Match!" << std::endl;
				return;
			}
		}
		std::cout << "\tNumber of faces: " << m_num_faces << std::endl;

		m_is_have_gt = true;

		m_is_gt = true;
		m_is_noise = false;
		m_is_denoised = false;

		m_data_manager->setOriginalMesh(m_gt_tri_mesh);

		m_gt_tri_mesh.request_face_normals();
		m_gt_tri_mesh.request_vertex_normals();
		m_gt_tri_mesh.update_normals();

		m_center = OpenMesh::Vec3d(0., 0., 0.);
		for (TriMesh::VertexIter v_it = m_gt_tri_mesh.vertices_begin(); v_it != m_gt_tri_mesh.vertices_end(); v_it++)
		{
			OpenMesh::Vec3d center_position_temp = m_gt_tri_mesh.point(*v_it);
			m_center += center_position_temp;
		}
		m_center /= int(m_gt_tri_mesh.n_vertices());

		m_max = 0;
		for (TriMesh::VertexIter v_it = m_gt_tri_mesh.vertices_begin(); v_it != m_gt_tri_mesh.vertices_end(); v_it++)
		{
			OpenMesh::Vec3d position_temp = m_gt_tri_mesh.point(*v_it) - m_center;
			for (int i = 0; i < 3; i++)
			{
				double temp_max = abs(position_temp[i]);
				if (temp_max > m_max)
				{
					m_max = temp_max;
				}
			}
		}

		for (TriMesh::VertexIter v_it = m_gt_tri_mesh.vertices_begin(); v_it != m_gt_tri_mesh.vertices_end(); v_it++)
		{
			m_gt_tri_mesh.point(*v_it) -= m_center;
			m_gt_tri_mesh.point(*v_it) /= m_max / 1.0;
		}

		m_is_have_gt = true;
	}
}

void MeshViewer::getFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, FaceNeighborType face_neighbor_type, std::vector<TriMesh::FaceHandle> &face_neighbor)
{
	face_neighbor.clear();
	if (face_neighbor_type == kEdgeBased)
	{
		for (TriMesh::FaceFaceIter ff_it = mesh.ff_iter(fh); ff_it.is_valid(); ff_it++)
			face_neighbor.push_back(*ff_it);
	}
	else if (face_neighbor_type == kVertexBased)
	{
		std::set<int> neighbor_face_index; neighbor_face_index.clear();

		for (TriMesh::FaceVertexIter fv_it = mesh.fv_begin(fh); fv_it.is_valid(); fv_it++)
		{
			for (TriMesh::VertexFaceIter vf_it = mesh.vf_iter(*fv_it); vf_it.is_valid(); vf_it++)
			{
				if ((*vf_it) != fh)
					neighbor_face_index.insert(vf_it->idx());
			}
		}

		for (std::set<int>::iterator iter = neighbor_face_index.begin(); iter != neighbor_face_index.end(); ++iter)
		{
			face_neighbor.push_back(TriMesh::FaceHandle(*iter));
		}
	}
}

void MeshViewer::getAllFaceNeighbor(TriMesh &mesh, std::vector<std::vector<TriMesh::FaceHandle> > &all_face_neighbor, FaceNeighborType face_neighbor_type, bool include_central_face)
{
	all_face_neighbor.resize(mesh.n_faces());
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		std::vector<TriMesh::FaceHandle> face_neighbor;
		getFaceNeighbor(mesh, *f_it, face_neighbor_type, face_neighbor);
		if (include_central_face) face_neighbor.push_back(*f_it);
		all_face_neighbor[f_it->idx()] = face_neighbor;
	}
}

void MeshViewer::getFaceArea(TriMesh &mesh, std::vector<double> &area)
{
	area.resize(mesh.n_faces());

	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		std::vector<TriMesh::Point> point;
		point.resize(3); int index = 0;
		for (TriMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
		{
			point[index] = mesh.point(*fv_it);
			index++;
		}
		TriMesh::Point edge1 = point[1] - point[0];
		TriMesh::Point edge2 = point[1] - point[2];
		double S = 0.5 * (edge1 % edge2).length();
		area[(*f_it).idx()] = S;
	}
}

void MeshViewer::getFaceCentroid(TriMesh &mesh, std::vector<TriMesh::Point> &centroid)
{
	centroid.resize(mesh.n_faces(), TriMesh::Point(0.0, 0.0, 0.0));
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		TriMesh::Point pt = mesh.calc_face_centroid(*f_it);
		centroid[(*f_it).idx()] = pt;
	}
}

void MeshViewer::slotLoadNoise(std::string noise_model_file)
{
	if (m_is_have_noise)
	{
		std::cout << "Existing noisy model, delete first." << std::endl;
		return;
	}

	if (noise_model_file == "")
	{
		return;
	}

	meshInitializedGet(noise_model_file.c_str(), false);
}

void MeshViewer::slotLoadGT(std::string gt_model_file)
{
	if (m_is_have_gt)
	{
		std::cout << "Existing gt model, delete first." << std::endl;
		return;
	}

	if (gt_model_file == "" || m_is_have_gt)
	{
		return;
	}

	meshInitializedGet(gt_model_file.c_str(), true);
}

void MeshViewer::slotDelete()
{
	if (!m_is_have_gt && !m_is_have_noise)
	{
		return;
	}

	if (m_is_have_noise)
	{
		m_noised_tri_mesh.clean();
		m_vertices_kd_tree.clear();
		m_all_face_neighbor.clear();
		m_face_area.clear();
		m_face_centroid.clear();

		m_is_have_noise = false;
	}

	if (m_is_have_gt)
	{
		m_gt_tri_mesh.clean();
		m_is_have_gt = false;
	}

	delete m_data_manager;
	m_data_manager = new DataManager();
	
	m_num_faces = 0;
	m_is_reload = true;
}

void MeshViewer::slotGenNoise(NoiseParams& noiseParams, int noiseDirection, std::string noise_type)
{
	m_is_have_noise = true;
	m_data_manager->MeshToOriginalMesh();

	Noise* noise;
	noise = new Noise(m_data_manager);

	//START
	auto meshtest = m_data_manager->getMesh();
	int vertNum = (int)meshtest.n_vertices();
	if (vertNum == 0) {
		return;
	}

	Noise::NoiseDirection noise_direction = (noiseDirection == 0) ? Noise::NoiseDirection::kNormal : Noise::NoiseDirection::kRandom;
	m_data_manager->MeshToOriginalMesh();
	TriMesh mesh = m_data_manager->getMesh();
	// compute average length of mesh
	double average_length = 0.0;
	for (TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
		average_length += mesh.calc_edge_length(*e_it);
	double edge_numbers = (double)mesh.n_edges();
	average_length /= edge_numbers;

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();

	if (noise_type == "Gaussian") {
		// add noise
		double standard_derivation = average_length * noiseParams.m_noiseLevel;

		std::vector<double> GaussianNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomGaussianNumbers(0, standard_derivation, (int)mesh.n_vertices(), GaussianNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

				TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * GaussianNumbers[v_it->idx()];
				mesh.set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh.n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh.point(*v_it) + RandomDirections[index] * GaussianNumbers[index];
				mesh.set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Impulsive") {
		// add noise
		double standard_derivation = average_length * noiseParams.m_noiseLevel;
		int impulsive_vertex_number = (int)(mesh.n_vertices() * noiseParams.m_impulsiveLevel);

		std::vector<double> GaussianNumbers;
		std::vector<TriMesh::Normal> RandomDirections;
		std::vector<std::pair<int, double> > VertexListAndGaussianNumbers;

		noise->randomImpulsiveNumbers(0, (int)mesh.n_vertices() - 1, impulsive_vertex_number, 0, standard_derivation, VertexListAndGaussianNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			for (int i = 0; i < (int)VertexListAndGaussianNumbers.size(); i++) {
				int index = VertexListAndGaussianNumbers[i].first;
				TriMesh::VertexHandle vh = mesh.vertex_handle(index);
				TriMesh::Point p = mesh.point(vh) + mesh.normal(vh) * VertexListAndGaussianNumbers[i].second;
				mesh.set_point(vh, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections(impulsive_vertex_number, RandomDirections);
			for (int i = 0; i < (int)VertexListAndGaussianNumbers.size(); i++) {
				int index = VertexListAndGaussianNumbers[i].first;
				TriMesh::VertexHandle vh = mesh.vertex_handle(index);
				TriMesh::Point p = mesh.point(vh) + RandomDirections[i] * VertexListAndGaussianNumbers[i].second;
				mesh.set_point(vh, p);
			}
		}
	}
	else if (noise_type == "Laplace") {
		std::vector<double> LaplaceNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomLaplaceNumbers(noiseParams.m_mu, noiseParams.m_b, (int)mesh.n_vertices(), LaplaceNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

				TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * LaplaceNumbers[v_it->idx()];
				mesh.set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh.n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh.point(*v_it) + RandomDirections[index] * LaplaceNumbers[index];
				mesh.set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Exponential") {
		std::vector<double> ExponentialNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomExponentialNumbers(noiseParams.m_lambda, (int)mesh.n_vertices(), ExponentialNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

				TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * ExponentialNumbers[v_it->idx()];
				mesh.set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh.n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh.point(*v_it) + RandomDirections[index] * ExponentialNumbers[index];
				mesh.set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "ExtremeValue") {
		std::vector<double> ExtremeValueNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomExtremeValueNumbers(noiseParams.m_a, noiseParams.m_b, (int)mesh.n_vertices(), ExtremeValueNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

				TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * ExtremeValueNumbers[v_it->idx()];
				mesh.set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh.n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh.point(*v_it) + RandomDirections[index] * ExtremeValueNumbers[index];
				mesh.set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Gamma") {
		std::vector<double> GammaNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomGammaNumbers(noiseParams.m_alpha, noiseParams.m_beta, (int)mesh.n_vertices(), GammaNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

				TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * GammaNumbers[v_it->idx()];
				mesh.set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh.n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh.point(*v_it) + RandomDirections[index] * GammaNumbers[index];
				mesh.set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "LogNormal") {
		std::vector<double> LogNormalNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomLogNormalNumbers(noiseParams.m_m, noiseParams.m_s, (int)mesh.n_vertices(), LogNormalNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

				TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * LogNormalNumbers[v_it->idx()];
				mesh.set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh.n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh.point(*v_it) + RandomDirections[index] * LogNormalNumbers[index];
				mesh.set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Uniform") {
		std::vector<double> UniformNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomUniformNumbers(noiseParams.m_a, noiseParams.m_b, (int)mesh.n_vertices(), UniformNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

				TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * UniformNumbers[v_it->idx()];
				mesh.set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh.n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh.point(*v_it) + RandomDirections[index] * UniformNumbers[index];
				mesh.set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Weibull") {
		std::vector<double> WeibullNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomWeibullNumbers(noiseParams.m_a, noiseParams.m_b, (int)mesh.n_vertices(), WeibullNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

				TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * WeibullNumbers[v_it->idx()];
				mesh.set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh.n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh.point(*v_it) + RandomDirections[index] * WeibullNumbers[index];
				mesh.set_point(*v_it, p);
			}
		}
	}
	
	//FINISH

	m_data_manager->setMesh(mesh);
	m_data_manager->setNoisyMesh(mesh);
	m_data_manager->setDenoisedMesh(mesh);

	m_noised_tri_mesh = mesh;

	m_noised_tri_mesh.request_face_normals();
	m_noised_tri_mesh.request_vertex_normals();
	m_noised_tri_mesh.update_normals();

	m_center = OpenMesh::Vec3d(0., 0., 0.);
	for (TriMesh::VertexIter v_it = m_noised_tri_mesh.vertices_begin(); v_it != m_noised_tri_mesh.vertices_end(); v_it++)
	{
		OpenMesh::Vec3d center_position_temp = m_noised_tri_mesh.point(*v_it);
		m_center += center_position_temp;
	}
	m_center /= int(m_noised_tri_mesh.n_vertices());

	m_max = 0;
	for (TriMesh::VertexIter v_it = m_noised_tri_mesh.vertices_begin(); v_it != m_noised_tri_mesh.vertices_end(); v_it++)
	{
		OpenMesh::Vec3d position_temp = m_noised_tri_mesh.point(*v_it) - m_center;
		for (int i = 0; i < 3; i++)
		{
			double temp_max = abs(position_temp[i]);
			if (temp_max > m_max)
			{
				m_max = temp_max;
			}
		}
	}

	for (TriMesh::VertexIter v_it = m_noised_tri_mesh.vertices_begin(); v_it != m_noised_tri_mesh.vertices_end(); v_it++)
	{
		m_noised_tri_mesh.point(*v_it) -= m_center;
		m_noised_tri_mesh.point(*v_it) /= m_max / 1.0;
		OpenMesh::Vec3d center_position_temp = m_noised_tri_mesh.point(*v_it);
		m_vertices_kd_tree.push_back(glm::vec3(center_position_temp[0], center_position_temp[1], center_position_temp[2]));
	}

	m_flann_kd_tree = shen::Geometry::EasyFlann(m_vertices_kd_tree);

	FaceNeighborType face_neighbor_type = kVertexBased;

	getAllFaceNeighbor(m_noised_tri_mesh, m_all_face_neighbor, face_neighbor_type, false);
	getFaceArea(m_noised_tri_mesh, m_face_area);
	getFaceCentroid(m_noised_tri_mesh, m_face_centroid);

	m_is_noise = true;
	m_is_denoised = false;
	m_is_gt = false;

	try {
		if (!OpenMesh::IO::write_mesh(m_noised_tri_mesh, m_noisyFile)) {
			std::cout << "Writing error!" << std::endl;
		}
	} catch(...) {
		std::cout << "External execption!" << std::endl;
	}
}