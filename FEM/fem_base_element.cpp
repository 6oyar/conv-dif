#ifndef _KRF_FEMBASE_ELEMENT_REALIZATIONS_H_
#define _KRF_FEMBASE_ELEMENT_REALIZATIONS_H_

#include <cstdlib>
#include <cstdio>
#include "fem_base_element.h"
#include "fem_basis_functions.h"

namespace FEM{
	namespace BaseElement{
		constexpr t_phi1D l1D[2]     = {l1D_1, l1D_2}; //!< linear basis functions array for the 1D-element
		constexpr t_phi1D dl1D[2][1] = {{dl1D_11}, {dl1D_21}}; //!< linear basis functions derivatives array for the 1D-element
		constexpr t_phi2D q2D[4]     = {q2D_1, q2D_2, q2D_3, q2D_4}; //!< polylinear basis functions array for the quadrangle element
		constexpr t_phi2D t2D[3]     = {t2D_1, t2D_2, t2D_3}; //!< polylinear basis functions array for the triangle element
		
		/*! polylinear basis functions partial derivatives array for the quadrangle element */
		constexpr t_phi2D dq2D[4][2] = {
			{dq2D_11, dq2D_12},
			{dq2D_21, dq2D_22},
			{dq2D_31, dq2D_32},
			{dq2D_41, dq2D_42} 
		};
		
		/*! polylinear basis functions partial derivatives array for the triangle element */
		constexpr t_phi2D dt2D[3][2] = {
			{dt2D_11, dt2D_12},
			{dt2D_21, dt2D_22},
			{dt2D_31, dt2D_32} 
		};

		constexpr double _nodesL[2][1] = {{0.0}, {1.0}}; //!< nodes of 1D element
		
		/*! nodes of triangle */
		constexpr double _nodesT[3][2] = {
			{0.0, 0.0},
			{1.0, 0.0},
			{0.0, 1.0}
		};
		
		/*! nodes of quadrangle */
		constexpr double _nodesQ[4][2] = {
			{0.0, 0.0},
			{1.0, 0.0},
			{1.0, 1.0},
			{0.0, 1.0}
		};
		
		/*! nodes of the 3D hexagonal element */
		constexpr double _nodesPrizmT[6][3] = {
			{0.0, 0.0, 0.0},
			{1.0, 0.0, 0.0},
			{0.0, 1.0, 0.0},
			{0.0, 0.0, 1.0},
			{1.0, 0.0, 1.0},
			{0.0, 1.0, 1.0} 
		};

		/*! nodes of the 3D octagonal element */
		constexpr double _nodesPrizmQ[8][3] = {
			{0.0, 0.0, 0.0},
			{1.0, 0.0, 0.0},
			{1.0, 1.0, 0.0},
			{0.0, 1.0, 0.0},
			{0.0, 0.0, 1.0},
			{1.0, 0.0, 1.0},
			{1.0, 1.0, 1.0},
			{0.0, 1.0, 1.0} 
		};
		
		const double* const nodesL[2] = {_nodesL[0], _nodesL[1]};
		const double* const nodesT[3] = {_nodesT[0], _nodesT[1], _nodesT[2]};
		const double* const nodesQ[4] = {_nodesQ[0], _nodesQ[1], _nodesQ[2], _nodesQ[3]};
		const double* const nodesPrizmT[6] = {_nodesPrizmT[0], _nodesPrizmT[1], _nodesPrizmT[2], _nodesPrizmT[3], _nodesPrizmT[4], _nodesPrizmT[5]};
		const double* const nodesPrizmQ[8] = {_nodesPrizmQ[0], _nodesPrizmQ[1], _nodesPrizmQ[2], _nodesPrizmQ[3], _nodesPrizmQ[4], _nodesPrizmQ[5], _nodesPrizmQ[6], _nodesPrizmQ[7]};


		double lPhi(int i, const double *point) noexcept{
			return l1D[i](point[0]);
		}
		double ldPhi(int i, int j, const double *point) noexcept{
			return dl1D[i][j](point[0]);
		}
		double tPhi(int i, const double *point) noexcept{
			return t2D[i](point[0], point[1]);
		}
		double tdPhi(int i, int j, const double *point) noexcept{
			return dt2D[i][j](point[0], point[1]);
		}
		double qPhi(int i, const double *point) noexcept{
			return q2D[i](point[0], point[1]);
		}
		double qdPhi(int i, int j, const double *point) noexcept{
			return dq2D[i][j](point[0], point[1]);
		}
		double prizmQPhi(int i, const double *point) noexcept{
			const int k = (i > 3) ? 1 : 0;
			const int j = i - 4*k;

			return q2D[j](point[0], point[1])*l1D[k](point[2]);
		}
		double prizmQdPhi(int i, int j, const double *point) noexcept{
			const int k = (i > 3) ? 1 : 0;
			const int l = i - 4*k;

			return (j < 2) ? dq2D[l][j](point[0], point[1])*l1D[k](point[2]) : q2D[l](point[0], point[1])*dl1D[k][0](point[2]);
		}
		double prizmTPhi(int i, const double *point) noexcept{
			const int k = (i > 2) ? 1 : 0;
			const int j = i - 3*k;

			return t2D[j](point[0], point[1])*l1D[k](point[2]);
		}
		double prizmTdPhi(int i, int j, const double *point) noexcept{
			const int k = (i > 2) ? 1 : 0;
			const int l = i - 3*k;

			return (j < 2) ? dt2D[l][j](point[0], point[1])*l1D[k](point[2]) : t2D[l](point[0], point[1])*dl1D[k][0](point[2]);
		}

		double const* const* elementNodes (TElementType et) noexcept{
			switch (et) {
				case etL: return nodesL;
				case etT: return nodesT;
				case etQ: return nodesQ;
				case etPrizmT: return nodesPrizmT;
				case etPrizmQ: return nodesPrizmQ;

				default: return NULL;
			}
		}

		bool elementBasis(TElementType et, t_phi &p, t_dphi &dp) noexcept{
			const TDim dim = elementDim(et);
			switch (dim) {
				case dim1D: { // 1D linear element
					p  = lPhi;
					dp = ldPhi;

					return true;
				}

				case dim2D: { // 2D polylinear element
					switch (et) {
						case etT: { // Triangle
							p  = tPhi;
							dp = tdPhi;
					
							return true;					
						}

						case etQ: { // Quadrangle
							p  = qPhi;
							dp = qdPhi;

							return true;
						}
						default: return false;
					}
				}

				case dim3D: { // 3D polylinear element
					switch (et) {
						case etPrizmT: { // Hexagonal
							p  = prizmTPhi;
							dp = prizmTdPhi;

							return true;
						}

						case etPrizmQ: { // Octagonal
							p  = prizmQPhi;
							dp = prizmQdPhi;

							return true;
						}
                        default: return false;
					}
				}

				default: return false;
			}
			
			return false;
		};

		int elementGraneNodeNumbers(TElementType et, int pos, size_t *nums) noexcept{
			const int n = elementGranesCount(et);
			if (pos < 0 || pos >= n){
				return 0;
			}

			switch (et) {
				case etL:{
					nums[0] = pos;
					nums[1] = 1 - pos;
					return 1;
				}
				
				case etT:
				case etQ:{
					if (pos < n-1){
						nums[0] = pos;
						nums[1] = pos+1;
						nums[2] = (pos == 0) ? 2 : 0;
						return 2;
					}
					nums[0] = n-1;
					nums[1] = 0;
					nums[2] = 1;
					return 2;
				}
				
				case etPrizmT:{
					if (pos < n-2){
						nums[0]	= pos;
						nums[1] = pos+1;
						nums[2] = pos+4;
						nums[3] = pos+3;
						nums[4] = pos+2;
						return 4;
					}
					
					const int plus = (pos == n-1) ? 0 : 3;
					for(int i=0; i<3; ++i){
						nums[i] = i+plus;
					}
					nums[3] = (plus) ? 0 : 3;
					return 3;
				}

				case etPrizmQ:{
					if (pos < n-2){
						nums[0]	= pos;
						nums[1] = pos+1;
						nums[2] = pos+5;
						nums[3] = pos+4;
						return 4;
					}
					
					const int plus = (pos == n-1) ? 0 : 4;
					for(int i=0; i<4; ++i){
						nums[i] = i+plus;
					}
					nums[4] = (plus) ? 0 : 4;
					return 4;
				}
				
				default: return 0;
			}
		}

		TElementType elementGrane(TElementType et, int pos, size_t *nums, Eigen::Vector3d &normal) noexcept{
			int ret = elementGraneNodeNumbers(et, pos, nums);
			if (!ret){
				return etUnknown;
			}

			switch (et) {
				case etL:{
					normal << ((pos) ? 1.0 : -1.0), 0, 0;
					return etUnknown;
				}
				
				case etT:{
					switch(pos) {
						case 0: normal <<  0.0, -1.0, 0.0; return etL;
						case 1: normal <<  1.0,  1.0, 0.0; return etL;
						case 2: normal << -1.0,  0.0, 0.0; return etL;
						default: return etUnknown;
					}
				}
				
				case etQ:{
					switch(pos) {
						case 0: normal <<  0.0, -1.0, 0.0; return etL;
						case 1: normal <<  1.0,  0.0, 0.0; return etL;
						case 2: normal <<  0.0,  1.0, 0.0; return etL;
						case 3: normal << -1.0,  0.0, 0.0; return etL;
						default: return etUnknown;
					}
				}
				
				case etPrizmT:{
					switch(pos) {
						case 0: normal <<  0.0, -1.0,  0.0; return etQ;
						case 1: normal <<  1.0,  1.0,  0.0; return etQ;
						case 2: normal << -1.0,  0.0,  0.0; return etQ;
						case 3: normal <<  0.0,  0.0, -1.0; return etT;
						case 4: normal <<  0.0,  0.0,  1.0; return etT;
						default: return etUnknown;
					}
				}

				case etPrizmQ:{
					switch(pos) {
						case 0: normal <<  0.0, -1.0,  0.0; return etQ;
						case 1: normal <<  1.0,  0.0,  0.0; return etQ;
						case 2: normal <<  0.0,  1.0,  0.0; return etQ;
						case 3: normal << -1.0,  0.0,  0.0; return etQ;
						case 4: normal <<  0.0,  0.0, -1.0; return etQ;
						case 5: normal <<  0.0,  0.0,  1.0; return etQ;
						default: return etUnknown;
					}
				}
				
				default: return etUnknown;
			}
		}

		TElementType elementGrane(TElementType et, int pos, size_t *nums) noexcept{
			return (elementGraneNodeNumbers(et, pos, nums)) ? elementGraneType(et, pos) : etUnknown;
		}

		TDim elementDim(TElementType et) noexcept{
			switch (et) {
				case etL: return dim1D;
				case etQ:
				case etT: return dim2D;
				case etPrizmQ:
				case etPrizmT: return dim3D;
				case etCount:
				default: return dimUnknown;
			}
		};

		int elementNodesCount(TElementType et) noexcept{
			switch (et) {
				case etL: return 2;
				case etT: return 3;
				case etQ: return 4;
				case etPrizmT: return 6;
				case etPrizmQ: return 8;
				default: return 0;
			}
		}

		double elementMeasure(TElementType et) noexcept{
			switch (et) {
				case etPrizmT:
				case etT: return 0.5;
				default: return 1.0;
			}
		}

		int elementGranesCount(TElementType et) noexcept{
			switch (et) {
				case etL: return 2;
				case etT: return 3;
				case etQ: return 4;
				case etPrizmT: return 5;
				case etPrizmQ: return 6;
				default: return 0;
			};
		}

		TElementType elementGraneType(TElementType et, int pos) noexcept{
			switch (et){
				case etL:{
					return etUnknown;
				}
				
				case etT:{
					switch(pos) {
						case 0: 
						case 1: 
						case 2: return etL;
						default: return etUnknown;
					}
				}
				
				case etQ:{
					switch(pos) {
						case 0: 
						case 1: 
						case 2: 
						case 3: return etL;
						default: return etUnknown;
					}
				}
				
				case etPrizmT:{
					switch(pos) {
						case 0: 
						case 1: 
						case 2: return etQ;
						case 3: 
						case 4: return etT;
						default: return etUnknown;
					}
				}

				case etPrizmQ:{
					switch(pos) {
						case 0: 
						case 1: 
						case 2: 
						case 3: 
						case 4: 
						case 5: return etQ;
						default: return etUnknown;
					}
				}
				
				default: return etUnknown;
			}
		}
	}; //namespace BaseElement
}; //namespace FEM
#endif
