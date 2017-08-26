#ifndef __KRF_FEMASSEMBLER_H__
#define __KRF_FEMASSEMBLER_H__
#include "noexcept_hack.h"

namespace FEM{
	using SpMat    = Eigen::SparseMatrix<double>;
	using SpTripl  = Eigen::Triplet<double>;
	using Vector3d = Eigen::Vector3d;
	using VectorXd = Eigen::VectorXd;
	using ArrayOfVector3d = std::vector<Vector3d,Eigen::aligned_allocator<Vector3d> >;
	
	class Assembler{
	public:
		static void div_k_grad(SpMat &A, CMesh &mesh, std::function<double(const CElement&)> k) noexcept{
			std::vector<SpTripl> triplets;
  			triplets.reserve(mesh.elementNodesCount);
			for(Mesh::TElement &e: mesh.elements){
				Integrator::localMatrix(*e.e);
				const int n = e.e->N;
				const double w = k(*(e.e));
				for(int i=0; i<n; ++i){
					for(int j=0; j<n; ++j){
						triplets.emplace_back(e.nodeNums[i], e.nodeNums[j], w*e.e->ae[i][j]);
					}
				}
			}

			A.setFromTriplets(triplets.begin(), triplets.end());
		};

		static void div_k_grad(SpMat &A, CMesh &mesh, const double k = 1.0) noexcept{
			std::vector<SpTripl> triplets;
  			triplets.reserve(mesh.elementNodesCount);
			for(Mesh::TElement &e: mesh.elements){
				Integrator::localMatrix(*e.e);
				const int n = e.e->N;
				for(int i=0; i<n; ++i){
					for(int j=0; j<n; ++j){
						triplets.emplace_back(e.nodeNums[i], e.nodeNums[j], k*e.e->ae[i][j]);
					}
				}
			}

			A.setFromTriplets(triplets.begin(), triplets.end());
		};

#if 0
		// This approach to the FEM-approximation of gradients must be tested
		static void grad_v(SpMat &G, const CMesh &mesh, std::function<double(const CElement&, Vector3d&)> v) noexcept{
			std::vector<SpTripl> triplets;
  			triplets.reserve(mesh.elementNodesCount);
  			const auto& elements = mesh.getElements();
  			Vector3d ve;
  			for(const auto& e: elements){
  				v(*e.e, ve);
  				_grad_v_element(triplets, mesh, e, ve);
  			}
  			G.setFromTriplets(triplets.begin(), triplets.end());
		}

		static void grad_v(SpMat &G, const CMesh &mesh, const ArrayOfVector3d& v) noexcept{
			std::vector<SpTripl> triplets;
  			triplets.reserve(mesh.elementNodesCount);
  			const auto& elements = mesh.getElements();
  			for(const auto& e: elements){
  				_grad_v_element(triplets, mesh, e, v[e.e->num]);
  			}
  			G.setFromTriplets(triplets.begin(), triplets.end());
		}

		static void grad_v(SpMat &G, const CMesh &mesh, const Vector3d &v) noexcept{
			std::vector<SpTripl> triplets;
  			triplets.reserve(mesh.elementNodesCount);
  			const auto& elements = mesh.getElements();
  			for(const auto& e: elements){
  				_grad_v_element(triplets, mesh, e, v);
  			}
  			G.setFromTriplets(triplets.begin(), triplets.end());
		}

		inline static void _grad_v_element(std::vector<SpTripl> &triplets, const CMesh &/*mesh*/, const Mesh::TElement &e, const Vector3d &v) noexcept{
			const int n    = e.e->N;
			const TDim dim = e.e->dim;
			for (int d=0; d<dim; ++d){
				const double vd = v[d];
				if (vd != 0){
					for (int i=0; i<n; ++i){
						for (int j=0; j<n; ++j){
							const double value = e.e->de[d][i][j]*vd;
							if (value != 0){
								triplets.emplace_back(e.nodeNums[i], e.nodeNums[j], value);
							}
						}
					}
				}
			}
		}
#endif

		static void mass(VectorXd &M, CMesh &mesh) noexcept{
			M.setZero(mesh.Nn);
			double *v = M.data();
			for(Mesh::TElement &e: mesh.elements){
				Integrator::localMatrix(*e.e);
				const int n = e.e->N;
				
				for(int i=0; i<n; ++i){
					v[e.nodeNums[i]] += e.e->me[i];
				}
			}
		};

		static void boundaryConditions(VectorXd &Mb, VectorXd &rhs, 
			const CBoundaryList &bList,
			const CMesh &mesh,
			bool isResetBefore = true) noexcept
		{
			if (isResetBefore){
				Mb.setZero(mesh.Nn);
				rhs.setZero(mesh.Nn);
			}

			for(const TBoundaryCond &bCnd: bList.condList){
				for(const Mesh::TGrane &g: bCnd.granes){
					const int n = g.element.N;
					for(int i=0; i<n; ++i){
						const double tmp = g.element.me[i];
						 Mb[g.nodeNums[i]] += tmp*bCnd.alpha;
						rhs[g.nodeNums[i]] += tmp*bCnd.beta;
					}
				}
			}
		};

		inline static void _FVMTriangle(std::vector<SpTripl> &triplets, 
			size_t n1, //!< global number of center node
			size_t n2, //!< global number of opposite node
			const Vector3d &A, //!< center node
			const Vector3d &B, const Vector3d &C, //!< BC - grane of finite volume assotiated with node A
			const Vector3d &v, //!< velocity vector,
			const CMesh &mesh
		) noexcept
		{
			auto CB = C-B;
			double value = (A-B).cross(CB).cross(CB).normalized().dot(v)*(CB.norm());
			if (value > 0){
				triplets.emplace_back(n1, n1,  value);
				triplets.emplace_back(n2, n1, -value);
			}
			else if (value < 0){
				triplets.emplace_back(n1, n2,  value);
				triplets.emplace_back(n2, n2, -value);
			}

			if (mesh.isBoundaryNode[n1] && mesh.isBoundaryNode[n2]){
				auto G = B-A;
				double value = (C-A).cross(G).cross(G).normalized().dot(v)*(G.norm());
				if (value != 0){
					triplets.emplace_back(n1, n1, value);
					triplets.emplace_back(n2, n2, value);
				}
			}
		}

		inline static void _FVMElement(std::vector<SpTripl> &triplets, 
			const CMesh &mesh, 
			const Mesh::TElement &e, 
			const Vector3d &v) noexcept
		{
			if (e.e->dim != dim2D){
				fprintf(stderr, "FVM is realized only in 2D\n");
				exit(1);
			}

			const size_t *nodeNums = e.nodeNums;
			const auto& nodes = mesh.getNodes();
			for(int i=0; i<e.e->N-1; ++i){
				_FVMTriangle(
					triplets,
					nodeNums[i], nodeNums[i+1],
					nodes[nodeNums[i]],
					0.5*(nodes[nodeNums[i]] + nodes[nodeNums[i+1]]),
					e.e->baricenter,
					v,
					mesh
				);
			}

			_FVMTriangle(
				triplets,
				nodeNums[e.e->N-1], nodeNums[0],
				nodes[nodeNums[e.e->N-1]],
				0.5*(nodes[nodeNums[e.e->N-1]] + nodes[nodeNums[0]]),
				e.e->baricenter,
				v,
				mesh
			);
		}

		static void FVM(SpMat &D, const CMesh &mesh, std::function<double(const CElement&, Vector3d&)> v) noexcept{
			std::vector<SpTripl> triplets;
  			triplets.reserve(mesh.elementNodesCount);
  			const auto& elements = mesh.getElements();
  			Vector3d ve;
  			for(const auto& e: elements){
  				v(*e.e, ve);
  				_FVMElement(triplets, mesh, e, ve);
  			}
  			D.setFromTriplets(triplets.begin(), triplets.end());
		}

		static void FVM(SpMat &D, const CMesh &mesh, const ArrayOfVector3d& v) noexcept{
			std::vector<SpTripl> triplets;
  			triplets.reserve(mesh.elementNodesCount);
  			const auto& elements = mesh.getElements();
  			for(const auto& e: elements){
  				_FVMElement(triplets, mesh, e, v[e.e->num]);
  			}
  			D.setFromTriplets(triplets.begin(), triplets.end());
		}

		static void FVM(SpMat &D, const CMesh &mesh, const Vector3d &v) noexcept{
			std::vector<SpTripl> triplets;
  			triplets.reserve(mesh.elementNodesCount);
  			const auto& elements = mesh.getElements();
  			for(const auto& e: elements){
  				_FVMElement(triplets, mesh, e, v);
  			}
  			D.setFromTriplets(triplets.begin(), triplets.end());
		}
	};
};
#endif
