#ifndef __KRF_FEMTYPEDEFS_H__
#define __KRF_FEMTYPEDEFS_H__

namespace FEM{
	using TDim = enum {
		dimUnknown = 0,
		dim1D = 1,
		dim2D = 2,
		dim3D = 3
	};

	using TElementType = enum {
		etUnknown = 0,
		etL = 1,   //!< segment
		etT,       //!< triangle
		etQ,       //!< quadrangle
		etPrizmT,  //!< 6-nodes prizm
		etPrizmQ,  //!< 8-nodes prizm
		etCount
	};

	//! 1D FEM elements basis function
	/*!
		\param x \\in [0,1] - coordinate

		\return \\phi(x)
	*/
	using t_phi1D = double (*)(double x);

	//! 2D FEM elements basis function
	/*!
		\param x \\in [0,1] - x-coordinate
		\param y \\in [0,1] - y-coordinate

		\return \\phi(x, y)
	*/
	using t_phi2D = double (*)(double x, double y);

	//! Abstract FEM elements basis function
	/*!
		\param  i - number of a function in elements basis
		\param point - vector of poin coordinates

		\return \\phi_i(point[0], point[1], ...)
	*/
	using t_phi = double (*)(int i, const double *point);

	//! Abstract FEM elements basis functions partial derivative
	/*!
		\param  i - number of a function in elements basis
		\param  j - number of derivative direction
		\param point - vector of poin coordinates

		\return \\partial{\\phi_i}{x_j}(point[0], point[1], ...)
	*/
	using t_dphi = double (*)(int i, int j, const double *point);
}; //namespace FEM

#endif
