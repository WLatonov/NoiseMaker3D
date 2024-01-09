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

void MeshProcessor::WriteNoisyMesh() {
	std::ofstream myfile;
	myfile.open(m_noisyFile);
	auto facesNum = m_mesh->n_faces();
	auto verticesNum = m_mesh->n_vertices();
	myfile << "# " + std::to_string(verticesNum) + " vertices"  + ", " + std::to_string(facesNum) + " faces\n";

	for (TriMesh::VertexIter vertexIter = m_mesh->vertices_begin(); vertexIter != m_mesh->vertices_end(); ++vertexIter) {
		auto myPoint = m_mesh->point(*vertexIter);
		float x = myPoint[0];
		float y = myPoint[1];
		float z = myPoint[2];
		myfile << "v " + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
	}

	for (TriMesh::FaceIter faceIter = m_mesh->faces_begin(); faceIter != m_mesh->faces_end(); ++faceIter) {
		myfile << "f";
		for (TriMesh::FaceVertexIter fv_it = m_mesh->fv_iter(*faceIter); fv_it.is_valid(); fv_it++)
		{
			auto vertexId = fv_it->idx();
			myfile << " " << std::to_string(vertexId + 1);
		}
		myfile << "\n";
	}

	myfile.close();
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

void MeshProcessor::SetSeed(int seed) { 
	m_seed = static_cast<unsigned>(seed);
	SetGenerator(m_seed);
};

void MeshProcessor::SetGenerator(unsigned seed) {
	m_generator = std::default_random_engine(seed);
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
	else if (noise_type == "Cauchy") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_x0) + "_" + std::to_string(noiseParams.m_gamma);
	}
	else if (noise_type == "Fisher") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_d1) + "_" + std::to_string(noiseParams.m_d2);
	}
	else if (noise_type == "Student") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_n);
	}
	else if (noise_type == "ChiSquared") {
		m_noisyFile += "_" + std::to_string(noiseParams.m_n);
	}

	if (noiseDirection == 1) {
		m_noisyFile += "_1";
	}
	m_noisyFile += ".obj";
}

void MeshProcessor::SetNoiseFileName(std::string noiseFileName) {
	m_noisyFile = noiseFileName;
}

std::vector<double> MeshProcessor::GetNoiseNumbers(std::string& noiseType, double averageLength, NoiseParams& noiseParams, Noise& noise) {
	std::vector<double> res;
	if (noiseType == "Gaussian") {
		noise.randomGaussianNumbers(m_seed, 0, averageLength * noiseParams.m_noiseLevel, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "Laplace") {
		noise.randomLaplaceNumbers(m_seed, noiseParams.m_mu, averageLength * noiseParams.m_b, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "Exponential") {
		noise.randomExponentialNumbers(m_seed, noiseParams.m_lambda / averageLength, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "ExtremeValue") {
		noise.randomExtremeValueNumbers(m_seed, noiseParams.m_a, averageLength * noiseParams.m_b, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "Gamma"){
		noise.randomGammaNumbers(m_seed, noiseParams.m_alpha, averageLength * noiseParams.m_beta, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "LogNormal") {
		noise.randomLogNormalNumbers(m_seed, noiseParams.m_m, noiseParams.m_s, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "Uniform") {
		noise.randomUniformNumbers(m_seed, averageLength * noiseParams.m_a, averageLength * noiseParams.m_b, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "Weibull") {
		noise.randomWeibullNumbers(m_seed, noiseParams.m_a, averageLength * noiseParams.m_b, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "Cauchy") {
		noise.randomCauchyNumbers(m_seed, noiseParams.m_x0, averageLength * noiseParams.m_gamma, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "Fisher") {
		noise.randomFisherNumbers(m_seed, noiseParams.m_d1, noiseParams.m_d2, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "Student") {
		noise.randomStudentNumbers(m_seed, noiseParams.m_n, (int)m_mesh->n_vertices(), res);
	}
	else if (noiseType == "ChiSquared") {
		noise.randomChiSquaredNumbers(m_seed, noiseParams.m_n, (int)m_mesh->n_vertices(), res);
	}

	return res;
}

void MeshProcessor::slotGenNoise(NoiseParams& noiseParams, int noiseDirection, std::string noise_type)
{
	SetNoiseFileName(noiseParams, noiseDirection, noise_type);
	Noise* noise;
	noise = new Noise();

	int vertNum = (int)m_mesh->n_vertices();
	if (vertNum == 0) {
		return;
	}

	Noise::NoiseDirection noise_direction = (noiseDirection == 0) ? Noise::NoiseDirection::kNormal : Noise::NoiseDirection::kRandom;
	double averageLength = 0.0;
	for (TriMesh::EdgeIter e_it = m_mesh->edges_begin(); e_it != m_mesh->edges_end(); e_it++)
		averageLength += m_mesh->calc_edge_length(*e_it);
	double edgeNumbers = (double)m_mesh->n_edges();
	averageLength /= edgeNumbers;

	m_mesh->request_face_normals();
	m_mesh->request_vertex_normals();
	m_mesh->update_normals();

	std::vector<TriMesh::Normal> randomDirections;

	if (noise_type == "Impulsive") {
		double standard_derivation = averageLength * noiseParams.m_noiseLevel;
		int impulsive_vertex_number = (int)(m_mesh->n_vertices() * noiseParams.m_impulsiveLevel);
		std::vector<std::pair<int, double> > VertexListAndGaussianNumbers;

		noise->randomImpulsiveNumbers(m_seed, 0, (int)m_mesh->n_vertices() - 1, impulsive_vertex_number, 0, standard_derivation, VertexListAndGaussianNumbers);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			for (int i = 0; i < (int)VertexListAndGaussianNumbers.size(); i++) {
				int index = VertexListAndGaussianNumbers[i].first;
				TriMesh::VertexHandle vh = m_mesh->vertex_handle(index);
				TriMesh::Point p = m_mesh->point(vh) + m_mesh->normal(vh) * VertexListAndGaussianNumbers[i].second;
				m_mesh->set_point(vh, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections(m_seed, impulsive_vertex_number, randomDirections);
			for (int i = 0; i < (int)VertexListAndGaussianNumbers.size(); i++) {
				int index = VertexListAndGaussianNumbers[i].first;
				TriMesh::VertexHandle vh = m_mesh->vertex_handle(index);
				TriMesh::Point p = m_mesh->point(vh) + randomDirections[i] * VertexListAndGaussianNumbers[i].second;
				m_mesh->set_point(vh, p);
			}
		}
	}
	else if(noise_type == "Student" || noise_type == "ChiSquared") {
		std::vector<double> noiseNumbers = GetNoiseNumbers(noise_type, averageLength, noiseParams, *noise);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {

				TriMesh::Point p = m_mesh->point(*v_it) + averageLength * noiseParams.m_scale * m_mesh->normal(*v_it) * noiseNumbers[v_it->idx()];
				m_mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections(m_seed, (int)m_mesh->n_vertices(), randomDirections);
			for (TriMesh::VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = m_mesh->point(*v_it) + averageLength * noiseParams.m_scale * randomDirections[index] * noiseNumbers[index];
				m_mesh->set_point(*v_it, p);
			}
		}
	} else {
		std::vector<double> noiseNumbers = GetNoiseNumbers(noise_type, averageLength, noiseParams, *noise);
		if (noise_direction == Noise::NoiseDirection::kNormal) {
			int i = 0;
			for (TriMesh::VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {

				TriMesh::Point p = m_mesh->point(*v_it) + m_mesh->normal(*v_it) * noiseNumbers[v_it->idx()];
				m_mesh->set_point(*v_it, p);
			}
		}
		else if (noise_direction == Noise::NoiseDirection::kRandom) {
			noise->randomDirections(m_seed, (int)m_mesh->n_vertices(), randomDirections);
			for (TriMesh::VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
				int index = v_it->idx();
				TriMesh::Point p = m_mesh->point(*v_it) + randomDirections[index] * noiseNumbers[index];
				m_mesh->set_point(*v_it, p);
			}
		}
	}

	m_mesh->request_face_normals();
	m_mesh->request_vertex_normals();
	m_mesh->update_normals();

	WriteNoisyMesh();
}

void MeshProcessor::MakePatchNoise(TriMesh& mesh, TopoNoiseParams& noiseParams) {
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
		delete patchData;
		if (!isTooDeleteByType) {
			continue;
		}
		deleteHandle.insert(*faceIter);
	}

	auto percentage = noiseParams.m_percentage;
	double probability = static_cast<double>(percentage / 100);
	std::vector<bool> isToDelete = GetBoolVectorWithSpecifiedProbabiliy(m_seed, deleteHandle.size(), probability);

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
	auto resId = GetRandomVertexId(m_seed, itemsNumber);
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
	int startVertId = GetRandomVertexId(m_seed, vertexNumber);
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
				if (currentVisitedVertices.find(*vertIter) != currentVisitedVertices.end()) {
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
				bool needRemove = GetBoolWithSpecifiedProbabiliy(m_generator, probabilityNew);
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
			bool needRemove = GetBoolWithSpecifiedProbabiliy(m_generator, probabilityNew);
			std::cout << distanceNew << "\t" << probabilityNew << "\t" << needRemove << std::endl;
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
	int vertNum = (int)m_mesh->n_vertices();
	if (vertNum == 0) {
		return;
	}

	m_mesh->request_vertex_normals();
	m_mesh->update_normals();

	if (noise_type == "RandomVertices") {
		double percentage = noiseParams.m_percentage;
		auto vertexNumber = m_mesh->n_vertices();
		double probability = static_cast<double>(percentage / 100);
		std::vector<int> vertexIndexToDelete = GetRandomBoolVector(m_seed, vertexNumber, probability);
		std::vector<TriMesh::VertexHandle> deleteHandle;

		for (TriMesh::VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++)
		{
			bool findSuccessfully = FindInVec(v_it->idx(), vertexIndexToDelete);
			if (findSuccessfully) {
				deleteHandle.push_back(*v_it);
			}
		}

		for (auto& vertHand : deleteHandle) {
			m_mesh->delete_vertex(vertHand, true);
		}
		m_mesh->garbage_collection();
	}
	else if (noise_type == "RandomFaces") {
		double percentage = noiseParams.m_percentage;
		auto facesNumber = m_mesh->n_faces();
		double probability = static_cast<double>(percentage / 100);
		std::vector<int> facesIndexToDelete = GetRandomBoolVector(m_seed, facesNumber, probability);
		std::vector<TriMesh::FaceHandle> deleteHandle;

		for (TriMesh::FaceIter v_it = m_mesh->faces_begin(); v_it != m_mesh->faces_end(); v_it++)
		{
			bool findSuccessfully = FindInVec(v_it->idx(), facesIndexToDelete);
			if (findSuccessfully) {
				deleteHandle.push_back(*v_it);
			}
		}

		for (auto& faceHand : deleteHandle) {
			m_mesh->delete_face(faceHand, true);
		}
		m_mesh->garbage_collection();
	}
	else if (noise_type == "OneRandomCluster") {
		auto vertexNumber = m_mesh->n_vertices();
		int startVertId = GetRandomVertexId(m_seed, vertexNumber);
		DeleteOneCluster(*m_mesh, startVertId, noiseParams.m_standadDeviation, noiseParams.m_divider, noiseParams.m_maxDistance);
	}
	else if (noise_type == "OneClusterInSpecifiedVertex") {
		DeleteOneCluster(*m_mesh, noiseParams.m_clusterCenterVertexId, noiseParams.m_standadDeviation, noiseParams.m_divider, noiseParams.m_maxDistance);
	}
	else if (noise_type == "SetOfRandomClusters") {
		DeleteSetOfClusters(*m_mesh, noiseParams);
	}
	else if (noise_type == "PatchModel") {
		MakePatchNoise(*m_mesh, noiseParams);
	}

	WriteNoisyMesh();
}