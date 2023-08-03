#include <memory>
#include <algorithm>
#include "MeshProcessor.h"

#include <omp.h>
#include <cstdio>
#include <ctime>

MeshProcessor::MeshProcessor(std::string currentFileName)
{
	m_gtFile = currentFileName;
	meshInitializedGet(currentFileName.c_str(), true);
	m_noisyFile = m_gtFile;
	m_noisyFile = m_noisyFile.substr(0, m_noisyFile.size() - 4);
}

void MeshProcessor::meshInitializedGet(const char* file_name, bool is_original)
{
	TriMesh* m_gtTriMesh = new TriMesh();
	OpenMesh::IO::read_mesh(*m_gtTriMesh, file_name);
	m_facesNumber = m_gtTriMesh->n_faces();
	std::cout << "Faces = " << m_facesNumber << std::endl;
	m_mesh = m_gtTriMesh;
}

void MeshProcessor::getFaceNeighbor(TriMesh &mesh, TriMesh::FaceHandle fh, std::vector<TriMesh::FaceHandle> &face_neighbor)
{
	face_neighbor.clear();
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

void MeshProcessor::getAllFaceNeighbor(TriMesh &mesh, std::vector<std::vector<TriMesh::FaceHandle> > &all_face_neighbor, bool include_central_face)
{
	all_face_neighbor.resize(mesh.n_faces());
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		std::vector<TriMesh::FaceHandle> face_neighbor;
		getFaceNeighbor(mesh, *f_it, face_neighbor);
		if (include_central_face) face_neighbor.push_back(*f_it);
		all_face_neighbor[f_it->idx()] = face_neighbor;
	}
}

void MeshProcessor::getFaceArea(TriMesh &mesh, std::vector<double> &area)
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

void MeshProcessor::getFaceCentroid(TriMesh &mesh, std::vector<TriMesh::Point> &centroid)
{
	centroid.resize(mesh.n_faces(), TriMesh::Point(0.0, 0.0, 0.0));
	for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		TriMesh::Point pt = mesh.calc_face_centroid(*f_it);
		centroid[(*f_it).idx()] = pt;
	}
}

void MeshProcessor::SetNoiseFileName(NoiseParams& noiseParams, int noiseDirection, std::string noise_type) {
	m_noisyFile += "_noise_" + noise_type;
	if (noise_type == "Gaussian") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_noiseLevel);
	}
	else if (noise_type == "Impulsive") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_noiseLevel)+ "_" + std::to_string(noiseParams.m_impulsiveLevel);
	}
	else if (noise_type == "Exponential") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_lambda);
	}
	else if (noise_type == "ExtremeValue") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_a) + "_" + std::to_string(noiseParams.m_b);
	}
	else if (noise_type == "Gamma") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_alpha) + "_" + std::to_string(noiseParams.m_beta);
	}
	else if (noise_type == "Laplace") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_mu) + "_" + std::to_string(noiseParams.m_b);
	}
	else if (noise_type == "LogNormal") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_m) + "_" + std::to_string(noiseParams.m_s);
	}
	else if (noise_type == "Uniform") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_a) + "_" + std::to_string(noiseParams.m_b);
	}
	else if (noise_type == "Weibull") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_a) + "_" + std::to_string(noiseParams.m_b);
	}

	if (noiseDirection == 1) {
		m_noisyFile += "_1";
	}
	m_noisyFile += ".obj";
}

void MeshProcessor::SetNoiseFileName(std::string noiseFileName) {
	m_noisyFile = noiseFileName;
}

void MeshProcessor::slotGenNoise(NoiseParams& noiseParams, int noiseDirection, std::string noise_type)
{
	SetNoiseFileName(noiseParams, noiseDirection, noise_type);
	Noise* noise;
	noise = new Noise();

	//START
	auto meshtest = m_mesh;
	int vertNum = (int)meshtest->n_vertices();
	if (vertNum == 0) {
		return;
	}

	Noise::NoiseDirection noise_direction = (noiseDirection == 0) ? Noise::NoiseDirection::kNormal : Noise::NoiseDirection::kRandom;
	TriMesh* mesh = meshtest;

	double average_length = 0.0;
	for (TriMesh::EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++)
		average_length += mesh->calc_edge_length(*e_it);
	double edge_numbers = (double)mesh->n_edges();
	average_length /= edge_numbers;

	mesh->request_face_normals();
	mesh->request_vertex_normals();
	mesh->update_normals();

	if (noise_type == "Gaussian") {
		// add noise
		double standard_derivation = average_length * noiseParams.m_noiseLevel;

		std::vector<double> GaussianNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomGaussianNumbers(0, standard_derivation, (int)mesh->n_vertices(), GaussianNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {

				TriMesh::Point p = mesh->point(*v_it) + mesh->normal(*v_it) * GaussianNumbers[v_it->idx()];
				mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh->n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh->point(*v_it) + RandomDirections[index] * GaussianNumbers[index];
				mesh->set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Impulsive") {
		// add noise
		double standard_derivation = average_length * noiseParams.m_noiseLevel;
		int impulsive_vertex_number = (int)(mesh->n_vertices() * noiseParams.m_impulsiveLevel);

		std::vector<double> GaussianNumbers;
		std::vector<TriMesh::Normal> RandomDirections;
		std::vector<std::pair<int, double> > VertexListAndGaussianNumbers;

		noise->randomImpulsiveNumbers(0, (int)mesh->n_vertices() - 1, impulsive_vertex_number, 0, standard_derivation, VertexListAndGaussianNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			for (int i = 0; i < (int)VertexListAndGaussianNumbers.size(); i++) {
				int index = VertexListAndGaussianNumbers[i].first;
				TriMesh::VertexHandle vh = mesh->vertex_handle(index);
				TriMesh::Point p = mesh->point(vh) + mesh->normal(vh) * VertexListAndGaussianNumbers[i].second;
				mesh->set_point(vh, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections(impulsive_vertex_number, RandomDirections);
			for (int i = 0; i < (int)VertexListAndGaussianNumbers.size(); i++) {
				int index = VertexListAndGaussianNumbers[i].first;
				TriMesh::VertexHandle vh = mesh->vertex_handle(index);
				TriMesh::Point p = mesh->point(vh) + RandomDirections[i] * VertexListAndGaussianNumbers[i].second;
				mesh->set_point(vh, p);
			}
		}
	}
	else if (noise_type == "Laplace") {
		std::vector<double> LaplaceNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomLaplaceNumbers(noiseParams.m_mu, average_length * noiseParams.m_b, (int)mesh->n_vertices(), LaplaceNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {

				TriMesh::Point p = mesh->point(*v_it) + mesh->normal(*v_it) * LaplaceNumbers[v_it->idx()];
				mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh->n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh->point(*v_it) + RandomDirections[index] * LaplaceNumbers[index];
				mesh->set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Exponential") {
		std::vector<double> ExponentialNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomExponentialNumbers(noiseParams.m_lambda / average_length, (int)mesh->n_vertices(), ExponentialNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {

				TriMesh::Point p = mesh->point(*v_it) + mesh->normal(*v_it) * ExponentialNumbers[v_it->idx()];
				mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh->n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh->point(*v_it) + RandomDirections[index] * ExponentialNumbers[index];
				mesh->set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "ExtremeValue") {
		std::vector<double> ExtremeValueNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomExtremeValueNumbers(noiseParams.m_a, average_length * noiseParams.m_b, (int)mesh->n_vertices(), ExtremeValueNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {

				TriMesh::Point p = mesh->point(*v_it) + mesh->normal(*v_it) * ExtremeValueNumbers[v_it->idx()];
				mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh->n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh->point(*v_it) + RandomDirections[index] * ExtremeValueNumbers[index];
				mesh->set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Gamma") {
		std::vector<double> GammaNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomGammaNumbers(noiseParams.m_alpha, average_length * noiseParams.m_beta, (int)mesh->n_vertices(), GammaNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {

				TriMesh::Point p = mesh->point(*v_it) + mesh->normal(*v_it) * GammaNumbers[v_it->idx()];
				mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh->n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh->point(*v_it) + RandomDirections[index] * GammaNumbers[index];
				mesh->set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "LogNormal") {
		std::vector<double> LogNormalNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomLogNormalNumbers(noiseParams.m_m, noiseParams.m_s, (int)mesh->n_vertices(), LogNormalNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {

				TriMesh::Point p = mesh->point(*v_it) + mesh->normal(*v_it) * LogNormalNumbers[v_it->idx()];
				mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh->n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh->point(*v_it) + RandomDirections[index] * LogNormalNumbers[index];
				mesh->set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Uniform") {
		std::vector<double> UniformNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomUniformNumbers(average_length * noiseParams.m_a, average_length * noiseParams.m_b, (int)mesh->n_vertices(), UniformNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {

				TriMesh::Point p = mesh->point(*v_it) + mesh->normal(*v_it) * UniformNumbers[v_it->idx()];
				mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh->n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh->point(*v_it) + RandomDirections[index] * UniformNumbers[index];
				mesh->set_point(*v_it, p);
			}
		}
	}
	else if (noise_type == "Weibull") {
		std::vector<double> WeibullNumbers;
		std::vector<TriMesh::Normal> RandomDirections;

		noise->randomWeibullNumbers(noiseParams.m_a, average_length * noiseParams.m_b, (int)mesh->n_vertices(), WeibullNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {

				TriMesh::Point p = mesh->point(*v_it) + mesh->normal(*v_it) * WeibullNumbers[v_it->idx()];
				mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections((int)mesh->n_vertices(), RandomDirections);
			for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = mesh->point(*v_it) + RandomDirections[index] * WeibullNumbers[index];
				mesh->set_point(*v_it, p);
			}
		}
	}
	
	//FINISH

	m_mesh->request_face_normals();
	m_mesh->request_vertex_normals();
	m_mesh->update_normals();

	try {
		if (!OpenMesh::IO::write_mesh(*m_mesh, m_noisyFile)) {
			std::cout << "Writing error!" << std::endl;
		}
	} catch(...) {
		std::cout << "External execption!" << std::endl;
	}
}

void MeshProcessor::MakeEpidemicNoise(TriMesh& mesh, TopoNoiseParams& noiseParams) {
	int ringsNumber = noiseParams.m_ringsNumber;
	int regionRadius = noiseParams.m_regionRadius;

	std::vector<std::vector<TriMesh::FaceHandle>> m_all_face_neighbor;
	std::vector<double> m_face_area;
	std::vector<TriMesh::Point> m_face_centroid;
	getAllFaceNeighbor(mesh, m_all_face_neighbor, false);
	getFaceArea(mesh, m_face_area);
	getFaceCentroid(mesh, m_face_centroid);

	std::set<PatchType> types = noiseParams.m_patchDeleteTypes;
	std::set<TriMesh::FaceHandle> deleteHandle;
	PatchData::SetEigenTresholdConstants(noiseParams.m_patchSelectionConstants);
	int vicinity = noiseParams.m_ringsNumberToDelete;

#pragma omp parallel for
	for (TriMesh::FaceIter faceIter = mesh.faces_begin(); faceIter != mesh.faces_end(); ++faceIter) {
		PatchData* patchData = new PatchData(mesh, m_all_face_neighbor, m_face_area, m_face_centroid, faceIter, ringsNumber, regionRadius);
		patchData->CalculateEigenValuesAndVectorsForVotingTensor();
		bool isTooDeleteByType = patchData->CheckIfOneOfTypes(types);
		if (!isTooDeleteByType) {
			continue;
		}
		deleteHandle.insert(*faceIter);
		delete patchData;
	}

	auto percentage = noiseParams.m_percentage;
	double probability = static_cast<double>(percentage / 100);
	std::vector<bool> isToDelete = GetBoolVectorWithSpecifiedProbabiliy(deleteHandle.size(), probability);

	int ind = 0;
	std::set<TriMesh::FaceHandle> deleteHandleNew;
	for (auto& handle : deleteHandle) {
		if (isToDelete[ind]) {
			deleteHandleNew.emplace(handle);
		}
		++ind;
	}
	deleteHandle = deleteHandleNew;

	for (int vic = 0; vic < vicinity; ++vic) {
		std::set<TriMesh::FaceHandle> deleteHandleNew;
		for (auto& handle : deleteHandle) {
			for (auto faceSecondIter = mesh.ff_iter(handle); faceSecondIter.is_valid(); ++faceSecondIter) {
				deleteHandleNew.insert(*faceSecondIter);
			}
			deleteHandleNew.insert(handle);
		}
		deleteHandle = deleteHandleNew;
	}

	for (auto& handle : deleteHandle) {
		mesh.delete_face(handle, true);
	}
	mesh.garbage_collection();
}

bool MeshProcessor::FindInVec(int index, std::vector<int>& vertexIndexToDelete) {
	for (int& vertIndex : vertexIndexToDelete) {
		if (vertIndex == index) {
			return true;
		}
	}
	return false;
}

TriMesh::VertexHandle MeshProcessor::FindNextClusterCenter(CenterSearchQueue& possibleNextClusterCenter) {
	auto pairFirst = possibleNextClusterCenter.top();
	int bestDistance = pairFirst.second;
	std::vector<TriMesh::VertexHandle> allBestHandles;
	while (!possibleNextClusterCenter.empty()) {
		auto pairCurrent = possibleNextClusterCenter.top();
		possibleNextClusterCenter.pop();
		int dist = pairCurrent.second;
		if (dist == bestDistance) {
			allBestHandles.push_back(pairCurrent.first);
		}
		else {
			break;
		}
	}
	auto itemsNumber = allBestHandles.size();
	auto resId = GetRandomVertexId(itemsNumber);
	return allBestHandles[resId];
}

void MeshProcessor::DeleteSetOfClusters(TriMesh& mesh, TopoNoiseParams& noiseParams) {
	std::vector<int> result;
	double standDev = noiseParams.m_standadDeviation;
	int clustersNumber = noiseParams.m_clustersNumber;
	int interClusterDist = noiseParams.m_interClusterDistance;
	bool makeClustersConcentrate = !noiseParams.m_makeClustersFar;
	double divider = noiseParams.m_divider;
	int maxDist = noiseParams.m_maxDistance;
	int maxSearchDist = interClusterDist * clustersNumber;
	auto vertexNumber = mesh.n_vertices();
	int startVertId = GetRandomVertexId(vertexNumber);
	TriMesh::VertexHandle startVertHandle;
	for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
		if (v_it->idx() == startVertId) {
			startVertHandle = *v_it;
		}
	}
	std::set<TriMesh::VertexHandle> deleteHandle;
	std::unordered_map<int, CenterSearchQueue> verticesAndDistances; // The distance from clusters centers to each vertex
	std::unordered_map<TriMesh::VertexHandle, int> totalDistMap;
	std::unordered_set<TriMesh::VertexHandle> visitedVertices;
	for (int clusterInd = 0; clusterInd < clustersNumber; ++clusterInd) {
		std::queue<CenterSearchQueuePair> bfsQueue;
		std::unordered_set<TriMesh::VertexHandle> currentVisitedVertices;
		bfsQueue.push(std::make_pair(startVertHandle, 0));
		deleteHandle.insert(startVertHandle);
		visitedVertices.insert(startVertHandle);

		CenterSearchQueue possibleNextClusterCenter = CenterSearchQueue(Comparator(makeClustersConcentrate));
		while (!bfsQueue.empty()) {
			auto pairCurrent = bfsQueue.front();
			bfsQueue.pop();
			auto currentVertexIter = pairCurrent.first;
			int distanceCurrent = pairCurrent.second;
			for (auto vertIter = mesh.vv_iter(currentVertexIter); vertIter.is_valid(); ++vertIter) {
				if (/*visitedVertices.find(*vertIter) != visitedVertices.end() ||*/
					currentVisitedVertices.find(*vertIter) != currentVisitedVertices.end()) {
					continue;
				}
				int distanceNew = distanceCurrent + 1;
				if (distanceNew > maxSearchDist) {
					continue;
				}
				if (distanceNew < interClusterDist) {
					visitedVertices.insert(*vertIter);
					totalDistMap.erase(*vertIter);
				}
				else {
					auto& totalDist = totalDistMap[*vertIter];
					totalDist += distanceNew;
					if (distanceNew == interClusterDist) {
						possibleNextClusterCenter.emplace(*vertIter, totalDist);
					}		
				}
				//Forming cluster
				double nonIntDistance = static_cast<double>(distanceNew) / divider;
				double probabilityNew = GetNormalDistributedProbability(nonIntDistance, standDev);
				bool needRemove = GetBoolWithSpecifiedProbabiliy(probabilityNew);
				if (needRemove && (distanceNew < maxDist)) {
					deleteHandle.insert(*vertIter);
				}
				bfsQueue.push(std::make_pair(*vertIter, distanceNew));
				currentVisitedVertices.insert(*vertIter);
			}
		}
		startVertHandle = FindNextClusterCenter(possibleNextClusterCenter);
	}

	for (auto& vertHand : deleteHandle) {
		mesh.delete_vertex(vertHand, true);
	}
	mesh.garbage_collection();
}

void MeshProcessor::DeleteOneCluster(TriMesh& mesh, int startVertId, double standardDeviation, double divider, int maxDist) {
	std::vector<TriMesh::VertexHandle> deleteHandle;
	std::queue<CenterSearchQueuePair> bfsQueue;
	std::unordered_set<int> visitedVertices;
	for (TriMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
		if (v_it->idx() == startVertId) {
			bfsQueue.push(std::make_pair(*v_it, 0));
			deleteHandle.push_back(*v_it);
			visitedVertices.insert(startVertId);
		}
	}

	while (!bfsQueue.empty()) {
		auto pairCurrent = bfsQueue.front();
		bfsQueue.pop();
		auto currentVertexIter = pairCurrent.first;
		int distanceCurrent = pairCurrent.second;
		for (auto vertIter = mesh.vv_iter(currentVertexIter); vertIter.is_valid(); ++vertIter) {
			auto vertId = vertIter->idx();
			if (visitedVertices.find(vertId) != visitedVertices.end()) {
				continue;
			}
			visitedVertices.insert(vertIter->idx());
			int distanceNew = distanceCurrent + 1;
			if (distanceNew >= maxDist) {
				continue;
			}
			double nonIntDistance = static_cast<double>(distanceNew) / divider;
			double probabilityNew = GetNormalDistributedProbability(nonIntDistance, standardDeviation);
			bool needRemove = GetBoolWithSpecifiedProbabiliy(probabilityNew);
			if (needRemove) {
				deleteHandle.push_back(*vertIter);
			}
			bfsQueue.push(std::make_pair(*vertIter, distanceNew));
			visitedVertices.insert(vertIter->idx());
		}
	}

	for (auto& vertHand : deleteHandle) {
		mesh.delete_vertex(vertHand, true);
	}
	mesh.garbage_collection();
}

void MeshProcessor::slotGenTopologyNoise(TopoNoiseParams& noiseParams, std::string noise_type)
{
	auto meshtest = m_mesh;
	int vertNum = (int)meshtest->n_vertices();
	if (vertNum == 0) {
		return;
	}

	TriMesh* mesh = meshtest;
	mesh->request_vertex_normals();
	mesh->update_normals();

	if (noise_type == "RandomVertices") {
		double percentage = noiseParams.m_percentage;
		auto vertexNumber = mesh->n_vertices();
		double probability = static_cast<double>(percentage / 100);
		std::vector<int> vertexIndexToDelete = GetRandomBoolVector(vertexNumber, probability);
		std::vector<TriMesh::VertexHandle> deleteHandle;

		for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
		{
			bool findSuccessfully = FindInVec(v_it->idx(), vertexIndexToDelete);
			if (findSuccessfully) {
				deleteHandle.push_back(*v_it);
			}
		}

		for (auto& vertHand : deleteHandle) {
			mesh->delete_vertex(vertHand, true);
		}
		mesh->garbage_collection();
	}
	else if (noise_type == "RandomFaces") {
		double percentage = noiseParams.m_percentage;
		auto facesNumber = mesh->n_faces();
		double probability = static_cast<double>(percentage / 100);
		std::vector<int> facesIndexToDelete = GetRandomBoolVector(facesNumber, probability);
		std::vector<TriMesh::FaceHandle> deleteHandle;

		for (TriMesh::FaceIter v_it = mesh->faces_begin(); v_it != mesh->faces_end(); v_it++)
		{
			bool findSuccessfully = FindInVec(v_it->idx(), facesIndexToDelete);
			if (findSuccessfully) {
				deleteHandle.push_back(*v_it);
			}
		}

		for (auto& faceHand : deleteHandle) {
			mesh->delete_face(faceHand, true);
		}
		mesh->garbage_collection();
	}
	else if (noise_type == "OneRandomCluster") {
		auto vertexNumber = mesh->n_vertices();
		int startVertId = GetRandomVertexId(vertexNumber);
		DeleteOneCluster(*mesh, startVertId, noiseParams.m_standadDeviation, noiseParams.m_divider, noiseParams.m_maxDistance);
	}
	else if (noise_type == "OneClusterInSpecifiedVertex") {
		DeleteOneCluster(*mesh, noiseParams.m_clusterCenterVertexId, noiseParams.m_standadDeviation, noiseParams.m_divider, noiseParams.m_maxDistance);
	}
	else if (noise_type == "SetOfRandomClusters") {
		DeleteSetOfClusters(*mesh, noiseParams);
	}
	else if (noise_type == "PatchModel") {
		MakeEpidemicNoise(*mesh, noiseParams);
	}

	try {
		if (!OpenMesh::IO::write_mesh(*mesh, m_noisyFile)) {
			std::cout << "Writing error!" << std::endl;
		}
	}
	catch (...) {
		std::cout << "External execption!" << std::endl;
	}
}