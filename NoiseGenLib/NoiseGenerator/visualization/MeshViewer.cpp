#include <memory>

#include "MeshViewer.h"
#include "FlannKDTree.h"
#include "Noise.h"
#include "MeshNormalFiltering.h"

#include <omp.h>
#include <cstdio>
#include <ctime>

MeshViewer::MeshViewer() : m_xRot(0), m_yRot(0), m_zRot(0)
{
	m_data_manager = new DataManager();
	m_current_file_name = "C:\\Users\\20962270\\source\\repos\\GCNMeshDenoiser19\\GCNMeshDenoiser19\\examples\\example_noisy.obj";
	meshInitializedGet(m_current_file_name.c_str(), false);
	meshInitializedGet("C:\\Users\\20962270\\source\\repos\\GCNMeshDenoiser19\\GCNMeshDenoiser19\\examples\\example.obj", true);
	int start = m_current_file_name.find_last_of('/');
	int end = m_current_file_name.find_last_of('.');
	m_model_name = m_current_file_name.substr(start + 1, end - start - 1);
	m_current_noise_level = -1.0;
	m_current_noise_type = -1;
}

MeshViewer::MeshViewer(std::string originalMeshDir, std::string noisedMeshDir, int var) : m_xRot(0), m_yRot(0), m_zRot(0)
{
	m_data_manager = new DataManager();
	meshInitializedGet(noisedMeshDir.c_str(), false);
	meshInitializedGet(originalMeshDir.c_str(), true);
	m_current_noise_level = -1.0;
	m_current_noise_type = -1;
	var = 0;
}

MeshViewer::MeshViewer(std::string currentFileName, std::string modelName, bool loadNoisedMesh) : m_xRot(0), m_yRot(0), m_zRot(0)
{
	m_data_manager = new DataManager();
	m_current_file_name = currentFileName;
	m_model_name = modelName;
	std::string fullMeshPath = m_current_file_name + "\\" + m_model_name + ".obj";
	meshInitializedGet(fullMeshPath.c_str(), true);
	if (loadNoisedMesh) {
		std::string fullNoisedMeshPath = m_current_file_name + "\\" + m_model_name + "_noisy.obj";
		meshInitializedGet(fullNoisedMeshPath.c_str(), false);
	}
	m_current_noise_level = -1.0;
	m_current_noise_type = -1;
}

MeshViewer::MeshViewer(std::string currentFileName) : m_xRot(0), m_yRot(0), m_zRot(0)
{
	m_data_manager = new DataManager();
	meshInitializedGet(currentFileName.c_str(), true);
	m_current_noise_level = -1.0;
	m_current_noise_type = -1;
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

		m_current_file_name = file_name;
		int start = m_current_file_name.find_last_of('/');
		int end = m_current_file_name.find_last_of('.');
		m_model_name = m_current_file_name.substr(start + 1, end - start - 1);

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

		m_current_file_name = file_name;
		int start = m_current_file_name.find_last_of('/');
		int end = m_current_file_name.find_last_of('.');
		m_model_name = m_current_file_name.substr(start + 1, end - start - 1);

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

void MeshViewer::slotGenNoise(double noise_level, std::string noise_type)
{
	if (!m_is_have_gt)
	{
		std::cout << "Load GT first!" << std::endl;
		return;
	}

	/*if (m_is_have_noise)
	{
		m_noised_tri_mesh.clean();
		m_vertices_kd_tree.clear();
		m_all_face_neighbor.clear();
		m_face_area.clear();
		m_face_centroid.clear();
	}*/

	m_is_have_noise = true;
	m_data_manager->MeshToOriginalMesh();

	Noise* noise;
	if (noise_type == "Gaussian")
	{
		noise = new Noise(m_data_manager, noise_level, 0);
	}
	else
	{
		noise = new Noise(m_data_manager, noise_level, 1);
	}
	//noise->addNoise();

	//START


	auto meshtest = m_data_manager->getMesh();
	int vertNum = (int)meshtest.n_vertices();
	if (vertNum == 0) {
		return;
	}

	int noise_type_index = 0;
	int noise_direction_index = 1;
	double impulsive_level = noise_level;
	m_data_manager->MeshToOriginalMesh();

	// compute average length of mesh
	double average_length = 0.0;
	TriMesh mesh = m_data_manager->getMesh();
	for (TriMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++)
		average_length += mesh.calc_edge_length(*e_it);
	double edge_numbers = (double)mesh.n_edges();
	average_length /= edge_numbers;

	// add noise
	double standard_derivation = average_length * noise_level;
	int impulsive_vertex_number = (int)(mesh.n_vertices() * impulsive_level);

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();

	std::vector<double> GaussianNumbers;
	std::vector<TriMesh::Normal> RandomDirections;
	std::vector<std::pair<int, double> > VertexListAndGaussianNumbers;
	std::cout << "Noise important numbers = " << edge_numbers << "\t" << standard_derivation << "\t" << impulsive_vertex_number << "\t" << std::distance(mesh.vertices_begin(), mesh.vertices_end()) << std::endl;

	noise->randomGaussianNumbers(0, standard_derivation, (int)mesh.n_vertices(), GaussianNumbers);
	int i = 0;
	for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {

		TriMesh::Point p = mesh.point(*v_it) + mesh.normal(*v_it) * GaussianNumbers[v_it->idx()];
		mesh.set_point(*v_it, p);
		++i;
		std::cout << "\t test = " << i << "\n";
	}

	m_data_manager->setMesh(mesh);
	m_data_manager->setNoisyMesh(mesh);
	m_data_manager->setDenoisedMesh(mesh);
	std::cout << "\t Test 1 " << "\n";

	int new_num_faces = mesh.n_faces();
	std::cout << "Faces = " << new_num_faces << std::endl;


	//FINISH

	//m_noised_tri_mesh = m_data_manager->getNoisyMesh();
	m_noised_tri_mesh = mesh;

	std::cout << "\t Test 2 " << "\n";

	new_num_faces = m_noised_tri_mesh.n_faces();
	std::cout << "Faces = " << new_num_faces << std::endl;

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

	//std::string fullNoisedMeshPath = m_current_file_name + "\\" + m_model_name + "_noisy.obj";
	std::string fullNoisedMeshPath = "C:\\Users\\20962270\\source\\repos\\GCNMeshDenoiser19\\GCNMeshDenoiser19\\examples\\example_noisy.obj";

	new_num_faces = m_noised_tri_mesh.n_faces();
	std::cout << "Faces = " << new_num_faces << std::endl;

	if (!OpenMesh::IO::write_mesh(m_noised_tri_mesh, "C:\\Users\\20962270\\source\\repos\\GCNMeshDenoiser19\\GCNMeshDenoiser19\\examples\\example_noisy1.obj" /*fullNoisedMeshPath*/)) {
		std::cout << "Writing error!" << std::endl;
	}
}