#ifndef __KRF_FEMBASEELEMENT_H__
#define __KRF_FEMBASEELEMENT_H__

#include "fem_typedefs.h"
#include <Eigen/Dense>

namespace FEM{
	//! FEM element basis structure
	/*!
		Contains basis functions and their partial derivatives for following polylinear elements:
		<ul>
			<li>1D element</li>
			<li>2D elements:
				<ul>
					<li>Triangle</li>
					<li>Quadrangle</li>
				</ul>
			</li>
			<li>3D elements:
				<ul>
					<li>Hexagonal</li>
					<li>Octagonal</li>
				</ul>
			</li>
		</ul>
	*/
	namespace BaseElement{
		//! Linear basis function for 1D element
		/*!
			\param i - number of basis function
			\param point - basis elements point coordinates

			\return \\phi_i(point[0])
		*/
		double lPhi(int i, const double *point) noexcept;
		
		//! Linear basis functions partial derivative for 1D element
		/*!
			\param i - number of basis function
			\param j - number of direction
			\param point - basis elements point coordinates

			\return value of \\partial{\\phi_i}{x_j}(point[0])
		*/
		double ldPhi(int i, int j, const double *point) noexcept;
		
		//! Bilinear basis function for triangle
		/*!
			\param i - number of basis function
			\param point - basis elements point coordinates

			\return \\phi_i(point[0], point[1])
		*/
		double tPhi(int i, const double *point) noexcept;
		
		//! Bilinear basis functions partial derivative for triangle
		/*!
			\param i - number of basis function
			\param j - number of direction
			\param point - basis elements point coordinates

			\return value of \\partial{\\phi_i}{x_j}(point[0], point[1])
		*/
		double tdPhi(int i, int j, const double *point) noexcept;
		
		//! Bilinear basis function for quadrangle
		/*!
			\param i - number of basis function
			\param point - basis elements point coordinates

			\return \\phi_i(point[0], point[1])
		*/	
		double qPhi(int i, const double *point) noexcept;
		
		//! Bilinear basis functions partial derivative for quadrangle
		/*!
			\param i - number of basis function
			\param point - basis elements point coordinates

			\return \\partial{\\phi_i}{x_j}(point[0], point[1])
		*/
		double qdPhi(int i, int j, const double *point) noexcept;
		
		//! Polylinear basis function for octagonal element
		/*!
			\param i - number of basis function
			\param point - basis elements point coordinates

			\return \\phi_i(point[0], point[1], point[2])
		*/
		double prizmQPhi(int i, const double *point) noexcept;
		
		//! Polylinear basis functions partial derivative for octagonal element
		/*!
			\param i - number of basis function
			\param point - basis elements point coordinates

			\return \\partial{\\phi_i}{x_j}(point[0], point[1], point[2])
		*/
		double prizmQdPhi(int i, int j, const double *point) noexcept;
		
		//! Polylinear basis function for hexagonal element
		/*!
			\param i - number of basis function
			\param point - basis elements point coordinates

			\return \\phi_i(point[0], point[1], point[2])
		*/
		double prizmTPhi(int i, const double *point) noexcept;
		
		//! Polylinear basis functions partial derivative for hexagonal element
		/*!
			\param i - number of basis function
			\param point - basis elements point coordinates

			\return \\partial{\\phi_i}{x_j}(point[0], point[1], point[2])
		*/
		double prizmTdPhi(int i, int j, const double *point) noexcept;

		//! Base element nodes ([0,1], [0,1]x[0,1], etc)
		/*!
			\param et - element type

			\return list of nodes with Decact coordinates
		*/
		const double *const *elementNodes (TElementType et) noexcept;

		//! Get pointers to the basis and partial derivatives functions by element type
		/*!
			\param et - element type
			\param p - pointer to the abstract basis function, given by reference
			\param dp - pointer to the abstract basis functions partial derivatives, given by reference

			\return true if success and false otherwise. 
			In success case parameters <i>p</i> and <i>dp</i> contains pointers to suitable realizations of abstract functions.
		*/ 
		bool elementBasis(TElementType et, t_phi &p, t_dphi &dp) noexcept;

		TDim elementDim(TElementType et) noexcept;
		int elementNodesCount(TElementType et) noexcept;
		double elementMeasure(TElementType et) noexcept;
		int elementGranesCount(TElementType et) noexcept;
		TElementType elementGraneType(TElementType et, int pos) noexcept;

		int elementGraneNodeNumbers(TElementType et, int pos, size_t *nums) noexcept;
		TElementType elementGrane(TElementType et, int pos, size_t *nums, Eigen::Vector3d &normal) noexcept;
		TElementType elementGrane(TElementType et, int pos, size_t *nums) noexcept;
	}; //namespace BaseElement
}; //namespace FEM

#endif
