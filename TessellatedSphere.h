//
// (c) 2019-2020 Takahiro Hashimoto
//
#pragma once

#include "VecMatDef.h"

class TessellatedSphere
{
public:
	struct IcosahedronFace
	{
		Idx3 vico;
		Idx3 eico;
		Bool3 edir;
		VecIdx v;
	};

	struct IcosahedronEdge
	{
		Idx2 e;
		VecIdx v;
	};

	TessellatedSphere();
	~TessellatedSphere();

	void Init(const int nt);

	const int& nt;
	const int& nv;
	const int& nf;
	const VecVec3& varr;
	const VecIdx3& farr;

	// Save 3D model of tesselated icosahedoron as VTK (ParaView) format
	void SaveAsVtk();

private:

	void CreateVertices();
	void CreateFace();
	void InitIcoVertices();
	void InitIcoEdges();
	void InitIcoFaces();

	// Spherical linear interpolation
	static Vec3 Slerp(const Vec3& s, const Vec3& e, const double ratio);

	int nt_; // Tesselate number
	int nv_; // Number of vertices
	int nf_; // Number of faces

	int ntri_; // Number of triangles in a single icosahedron face

	// Vertices and faces on tessalated sphere
	VecVec3 vsphere_;
	VecIdx3 fsphere_;

	VecVec3 vico_; // Icosahedron vertices
	std::vector<IcosahedronFace> fico_;
	std::vector<IcosahedronEdge> eico_;

	static const int IcoVert = 12;
	static const int IcoEdge = 30;
	static const int IcoFace = 20;
};
