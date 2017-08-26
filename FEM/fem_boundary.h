#ifndef _KRF_FEMBOUNDARY_H
#define _KRF_FEMBOUNDARY_H
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include "fem_mesh.h"
#include "noexcept_hack.h"
namespace FEM{
	using TBoundaryCondType = enum {
		bndCondNone = 0,
		bndCondDirichlet,  //!< u = u_0
		bndCondNeumann,    //!< -k du/dn = q
		bndCondConvection, //!< -k du/dn = h(u-u_{ext})
		bndCondRadiation,  //!< -k du/dn = \sigma\epsilon(u^4 - u^4_{ext})
		bndCondRobin       //!< -k du/dn = \alpha u - \beta
	};

	using Vector3d = Eigen::Vector3d;

	//! Boundary condition -k*du/dn = alpha*u-beta
	class TBoundaryCond{
	public:
		const size_t id;
		const TBoundaryCondType type; //!< condition type 
		const double alpha;
		const double beta;
		
		const Eigen::Vector3d n; //!< outward normal vector

		std::vector<std::reference_wrapper<const Mesh::TGrane> >granes;

		void addGrane(const Mesh::TGrane &g) noexcept{
			granes.emplace_back(g);
		}

		TBoundaryCond(size_t _id, TBoundaryCondType _type, double _alpha, double _beta, const Vector3d &_n) noexcept:
		    id(_id), type(_type), alpha(_alpha), beta(_beta), n(_n){};
	};

	//! boundary conditions
	class CBoundaryList{
	public:
		std::vector<TBoundaryCond> condList; //!< list of boundary conditions

		//!< get boundary condition for grane by outward normal vector
		/*!
			\param[in] normal - granes outward normal
			\return pair<boundary index, true> or pair<0, false> if condition was not found
		*/
		using search_res_t = std::pair<size_t, bool>;
		search_res_t boundaryIndex(const Vector3d &normal) const noexcept{
			constexpr double eps  = 1e-8;
			constexpr double eps2 = eps*eps;
			for (const TBoundaryCond &b: condList){
				if (normal.squaredNorm() > eps2 && (normal - b.n).squaredNorm() < eps2){
					return search_res_t(b.id, true);
				}
			}
			return search_res_t(0, false);
		}

		void addGraneByNormal(const Mesh::TGrane &g) noexcept{
			const auto search = boundaryIndex(g.normal);
			if (search.second){
				condList[search.first].addGrane(g);
			}
		}

		void addGraneListByNormal(const std::vector<Mesh::TGrane> &gList) noexcept{
			for(const auto& g: gList){
				addGraneByNormal(g);
			}
		}

		//! add boundary condition -k*du/dn = alpha*u - beta
		/*!
			\param[in] type - condition type
			\param[in] alpha
			\param[in] beta 
			\param[in] n - outward normal vector

			\return boundary index or 0 on error 
		*/
		size_t add(TBoundaryCondType type, double alpha, double beta, const Vector3d &n) noexcept{
			condList.emplace_back(condList.size(), type, alpha, beta, n);

			return condList.size();
		}

		size_t addDirichlet(double u0, const Vector3d &n) noexcept{
			return add(bndCondDirichlet, 1e+16, -1e+16*u0, n);
		}

		size_t addNeumann(double q, const Vector3d &n) noexcept{
			return add(bndCondNeumann, 0.0, q, n);
		}

		size_t addConvection(double h, double uExt, const Vector3d &n) noexcept{
			return add(bndCondConvection, h, -h*uExt, n);
		}

		size_t addCondRadiation(double sigmaEpsilon, double uExt, const Vector3d &n) noexcept{
			return add(bndCondRadiation, sigmaEpsilon, -sigmaEpsilon*uExt, n);
		}
	};
}; //namespace FEM
#endif
