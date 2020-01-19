//
// (c) 2019-2020 Takahiro Hashimoto
//

#include "TessellatedSphere.h"
#include <fstream>

TessellatedSphere::TessellatedSphere()
: nt(nt_), nv(nv_), nf(nf_), varr(vsphere_), farr(fsphere_)
{
	using namespace std;

	InitIcoVertices();
	InitIcoEdges();
	InitIcoFaces();
}

TessellatedSphere::~TessellatedSphere()
{
}

void TessellatedSphere::Init(const int nt)
{
	nt_ = nt;

	CreateVertices();
	CreateFace();
}

void TessellatedSphere::SaveAsVtk() const
{
	using namespace std;

	const std::string fname = "Sphere.vtk";
	std::ofstream ofs;
	ofs.open(fname);

	// write header 
	ofs << "# vtk DataFile Version 2.0" << endl;
	ofs << "3D model of sphere" << endl;
	ofs << "ASCII" << endl;
	ofs << endl;
	ofs << "DATASET POLYDATA" << endl;

	ofs << "POINTS " << nv_ << " float" << endl;
	for (const Vec3& v : vsphere_)
	{
		ofs << " "
			<< v.x() << " "
			<< v.y() << " "
			<< v.z() << endl;
	}
	
	ofs << endl;

	ofs << "POLYGONS " << nf_ << " " << 4 * nf_ << endl;
	for (const Idx3& f : fsphere_)
	{
		ofs << " 3 " 
			<< f.x() << " "
			<< f.y() << " "
			<< f.z() << endl;
	}
	ofs.close();
}

void TessellatedSphere::CreateVertices()
{
	// number of vertices
	nv_ = 10 * nt_*nt_ + 2;

	vsphere_.reserve(nv_);

	// push vertices of original icosahedron
	for(Vec3 & v : vico_)
	{
		vsphere_.push_back(v);
	}

	// push tesselated vertices in icosahedron edges
	if (nt_ >= 2)
	{
		/// iterate with edges
		for (int iedge = 0; iedge < eico_.size(); iedge++)
		{
			const Vec3 s(vico_[eico_[iedge].e[0]]);
			const Vec3 e(vico_[eico_[iedge].e[1]]);

			for (int it = 1; it < nt_; it++)
			{
				vsphere_.push_back(Slerp(s, e, (double)it / (double)nt_));
				eico_[iedge].v[it-1] = (int)vsphere_.size()-1;
			}
		}
	}

	// push tesselated vertices inside icosahedron faces
	if (nt_ >= 3)
	{
		// vertices inside tesselated icosahedron face;
		const int nvface = (nt_ + 1)*(nt_ + 2) / 2;
		VecVec3 vface(nvface);
		
		// iterate with faces
		for (int iface = 0; iface < fico_.size(); iface++)
		{
			/// vertices of original face
			const Vec3 v1(vico_[fico_[iface].vico[0]]);
			const Vec3 v2(vico_[fico_[iface].vico[1]]);
			const Vec3 v3(vico_[fico_[iface].vico[2]]);

			// create tesselated vertices on edges
			for (int it = 1; it < nt_; it++)
			{
				// index of left edge
				const int ileft = it*(it + 1) / 2;
				vface[ileft] = Slerp(v1, v2, (double)it / (double)nt_);
				// index of right edge
				const int iright = (it + 1)*(it + 2) / 2 - 1;
				vface[iright] = Slerp(v1, v3, (double)it / (double)nt_);
			}

			// create tesselated vertices inside icosahedron face
			for (int irow = 3; irow < nt_+1; irow++)
			{
				for (int icol = 1; icol <= irow -2; icol++)
				{
					const int istart = irow*(irow - 1) / 2;
					const int iend = istart + irow - 1;
					const Vec3 vstart(vface[istart]);
					const Vec3 vend(vface[iend]);
					vsphere_.push_back(Slerp(vstart, vend, icol / (double)(irow - 1)));
					fico_[iface].v[irow*(irow-1)/2 + icol] = (int)vsphere_.size() - 1;
				}
			}
		}
	}
}

void TessellatedSphere::CreateFace()
{
	ntri_ = nt_*nt_;
	nf_ = 20 * ntri_;
	fsphere_.reserve(nf_);

	//set vindices on vertices
	for (int iface = 0; iface < IcoFace; iface++)
	{
		fico_[iface].v[0] = fico_[iface].vico.x();
		fico_[iface].v[nt_*(nt_+1)/2] = fico_[iface].vico.y();
		fico_[iface].v[(nt_ +1)*(nt_ +2)/2 -1] = fico_[iface].vico.z();
	}

	// set vindices on edges
	for (int iface = 0; iface < IcoFace; iface++)
	{
		// set vindices on left edge
		const int e1 = fico_[iface].eico.x();
		if ( fico_[iface].edir[0] )
		{
			for (int it = 2; it <= nt_; it++)
			{
				fico_[iface].v[it*(it-1)/2] = eico_[e1].v[it-2];
			}
		}
		else
		{
			for (int it = 2; it <= nt_; it++)
			{
				fico_[iface].v[it*(it-1)/2] = eico_[e1].v[nt_-it];
			}
		}

		// set vindices on bottom edge
		const int e2 = fico_[iface].eico.y();
		if (fico_[iface].edir[1])
		{
			for (int it = 2; it <= nt_; it++)
			{
				fico_[iface].v[nt_*(nt_ +1)/2 +it -1] = eico_[e2].v[it-2];
			}
		}
		else
		{
			for (int it = 2; it <= nt_; it++)
			{
				fico_[iface].v[nt_*(nt_ +1)/2 + it - 1] = eico_[e2].v[nt_-it];
			}
		}

		// set vindices on right edge
		const int e3 = fico_[iface].eico.z();
		if (fico_[iface].edir[2])
		{
			for (int it = 2; it <= nt_; it++)
			{
				fico_[iface].v[it*(it+1)/2-1] = eico_[e3].v[nt_-it];
			}
		}
		else
		{
			for (int it = 2; it <= nt_; it++)
			{
				fico_[iface].v[it*(it+1)/2-1] = eico_[e3].v[it-2];
			}
		}
	}

	for (int iface = 0; iface < IcoFace; iface++)
	{
		// create upright triangle in single icosahedron face
		for (int irow = 1; irow <= nt; irow++)
		{
			for (int icol = 1; icol <= irow; icol++)
			{	
				const int v1 = fico_[iface].v[irow*(irow-1)/2 +icol -1];
				const int v2 = fico_[iface].v[irow*(irow+1)/2 +icol -1];
				const int v3 = fico_[iface].v[irow*(irow+1)/2 +icol];
				const Idx3 fidx(v1, v2, v3);
				fsphere_.push_back(fidx);
			}
		}

		// create reverse triangle in single icosahedron face
		for (int irow = 2; irow <= nt; irow++)
		{
			for (int icol = 1; icol <= irow -1; icol++)
			{
				const int v1 = fico_[iface].v[irow*(irow-1)/2 +icol -1];
				const int v2 = fico_[iface].v[irow*(irow+1)/2 +icol];
				const int v3 = fico_[iface].v[irow*(irow-1)/2 +icol];
				const Idx3 fidx(v1, v2, v3);
				fsphere_.push_back(fidx);
			}
		}

	}
}

void TessellatedSphere::InitIcoVertices()
{
	const double a = 1.0 / sqrt(5.0);
	const double b = (1.0 - a) / 2.0;
	const double c = (1.0 + a) / 2.0;
	const double d = sqrt(b);
	const double e = sqrt(c);

	vico_.resize(IcoVert);

	vico_[0] = Vec3(0, 0, 1.0);
	vico_[1] = Vec3(2*a, 0, a);
	vico_[2] = Vec3(b, e, a);
	vico_[3] = Vec3(-c, d, a);
	vico_[4] = Vec3(-c, -d, a);
	vico_[5] = Vec3(b, -e, a);
	vico_[6] = Vec3(c, d, -a);
	vico_[7] = Vec3(-b, e, -a);
	vico_[8] = Vec3(-2 * a, 0, -a);
	vico_[9] = Vec3(-b, -e, -a);
	vico_[10] = Vec3(c, -d, -a);
	vico_[11] = Vec3(0, 0, -1.0);
}

void TessellatedSphere::InitIcoEdges()
{
	eico_.resize(IcoEdge);

	eico_[0].e = Idx2(0, 1);
	eico_[1].e = Idx2(0, 2);
	eico_[2].e = Idx2(0, 3);
	eico_[3].e = Idx2(0, 4);
	eico_[4].e = Idx2(0, 5);
	eico_[5].e = Idx2(1, 2);
	eico_[6].e = Idx2(2, 3);
	eico_[7].e = Idx2(3, 4);
	eico_[8].e = Idx2(4, 5);
	eico_[9].e = Idx2(5, 1);
	eico_[10].e = Idx2(1, 6);
	eico_[11].e = Idx2(6, 2);
	eico_[12].e = Idx2(2, 7);
	eico_[13].e = Idx2(7, 3);
	eico_[14].e = Idx2(3, 8);
	eico_[15].e = Idx2(8, 4);
	eico_[16].e = Idx2(4, 9);
	eico_[17].e = Idx2(9, 5);
	eico_[18].e = Idx2(5, 10);
	eico_[19].e = Idx2(10, 1);
	eico_[20].e = Idx2(6, 7);
	eico_[21].e = Idx2(7, 8);
	eico_[22].e = Idx2(8, 9);
	eico_[23].e = Idx2(9, 10);
	eico_[24].e = Idx2(10, 6);
	eico_[25].e = Idx2(6, 11);
	eico_[26].e = Idx2(7, 11);
	eico_[27].e = Idx2(8, 11);
	eico_[28].e = Idx2(9, 11);
	eico_[29].e = Idx2(10, 11);

	for (int iedge = 0; iedge < IcoEdge; iedge++)
	{
		eico_[iedge].v.resize(nt_ - 1);
	}
}

void TessellatedSphere::InitIcoFaces()
{
	fico_.resize(IcoFace);

	// order of indices are defined, so that the normals direct outward
	fico_[0].vico  = Idx3(0, 1, 2);
	fico_[1].vico = Idx3(0, 2, 3);
	fico_[2].vico = Idx3(0, 3, 4);
	fico_[3].vico = Idx3(0, 4, 5);
	fico_[4].vico = Idx3(0, 5, 1);
	fico_[5].vico = Idx3(1, 6, 2);
	fico_[6].vico = Idx3(2, 6, 7);
	fico_[7].vico = Idx3(2, 7, 3);
	fico_[8].vico = Idx3(3, 7, 8);
	fico_[9].vico = Idx3(3, 8, 4);
	fico_[10].vico = Idx3(4, 8, 9);
	fico_[11].vico = Idx3(4, 9, 5);
	fico_[12].vico = Idx3(5, 9, 10);
	fico_[13].vico = Idx3(5, 10, 1);
	fico_[14].vico = Idx3(1, 10, 6);
	fico_[15].vico = Idx3(11, 7, 6);
	fico_[16].vico = Idx3(11, 8, 7);
	fico_[17].vico = Idx3(11, 9, 8);
	fico_[18].vico = Idx3(11, 10, 9);
	fico_[19].vico = Idx3(11, 6, 10);

	fico_[0].eico = Idx3(0, 5, 1); fico_[0].edir = Bool3(true, true, false);
	fico_[1].eico = Idx3(1, 6, 2); fico_[1].edir = Bool3(true, true, false);
	fico_[2].eico = Idx3(2, 7, 3); fico_[2].edir = Bool3(true, true, false);
	fico_[3].eico = Idx3(3, 8, 4); fico_[3].edir = Bool3(true, true, false);
	fico_[4].eico = Idx3(4, 9, 0); fico_[4].edir = Bool3(true, true, false);
	fico_[5].eico = Idx3(10, 11, 5); fico_[5].edir = Bool3(true, true, false);
	fico_[6].eico = Idx3(11, 20, 12); fico_[6].edir = Bool3(false, true, false);
	fico_[7].eico = Idx3(12, 13, 6); fico_[7].edir = Bool3(true, true, false);
	fico_[8].eico = Idx3(13, 21, 14); fico_[8].edir = Bool3(false, true, false);
	fico_[9].eico = Idx3(14, 15, 7); fico_[9].edir = Bool3(true, true, false);
	fico_[10].eico = Idx3(15, 22, 16); fico_[10].edir = Bool3(false, true, false);
	fico_[11].eico = Idx3(16, 17, 8); fico_[11].edir = Bool3(true, true, false);
	fico_[12].eico = Idx3(17, 23, 18); fico_[12].edir = Bool3(false, true, false);
	fico_[13].eico = Idx3(18, 19, 9); fico_[13].edir = Bool3(true, true, false);
	fico_[14].eico = Idx3(19, 24, 10); fico_[14].edir = Bool3(false, true, false);
	fico_[15].eico = Idx3(26, 20, 25); fico_[15].edir = Bool3(false, false, true);
	fico_[16].eico = Idx3(27, 21, 26); fico_[16].edir = Bool3(false, false, true);
	fico_[17].eico = Idx3(28, 22, 27); fico_[17].edir = Bool3(false, false, true);
	fico_[18].eico = Idx3(29, 23, 28); fico_[18].edir = Bool3(false, false, true);
	fico_[19].eico = Idx3(25, 24, 29); fico_[19].edir = Bool3(false, false, true);
	
	for (int i = 0; i < IcoFace; i++)
	{
		fico_[i].v.resize((nt_ + 1)*(nt_ + 2)/2);
	}
}

Vec3 TessellatedSphere::Slerp(const Vec3& s, const Vec3& e, const double ratio)
{
	const double cos_omega = s.dot(e);
	const double omega = acos(cos_omega);
	const double sin_omega = sqrt(1. - cos_omega * cos_omega);
	const double a = sin((1.0 - ratio) * omega) / sin_omega;
	const double b = sin(ratio * omega) / sin_omega;

	const Vec3 p = a * s + b * e;
	return p;
}
