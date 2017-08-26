#ifndef _KRF_FEMBASIS_FUNCTIONSDEFS_H_
#define _KRF_FEMBASIS_FUNCTIONSDEFS_H_
#include "noexcept_hack.h"

namespace FEM{
	namespace BaseElement{
		// 1D elements
		//! First basis function of 1D linear FEM element
		double l1D_1(double x) noexcept;

		//! Second basis function of 1D linear FEM element
		double l1D_2(double x) noexcept;

		//! Derivative of the first basis function of 1D linear FEM element
		constexpr double dl1D_11(double /*x*/) noexcept{
			return -1.0;
		}

		//! Derivative of the second basis function of 1D linear FEM element
		constexpr double dl1D_21(double /*x*/) noexcept{
			return 1.0;
		}

		// 2D elements
		// Quadrangle
		//! First basis function of quadrangle FEM element
		double q2D_1(double x, double y) noexcept;

		//! Second basis function of quadrangle FEM element
		double q2D_2(double x, double y) noexcept;

		//! 3-th basis function of quadrangle FEM element
		double q2D_3(double x, double y) noexcept;

		//! 4-th basis function of quadrangle FEM element
		double q2D_4(double x, double y) noexcept;

		//! x-derivative of the first basis function of quadrangle FEM element
		double dq2D_11(double x, double y) noexcept;

		//! y-derivative of the first basis function of quadrangle FEM element
		double dq2D_12(double x, double y) noexcept;

		//! x-derivative of the second basis function of quadrangle FEM element
		double dq2D_21(double x, double y) noexcept;

		//! y-derivative of the second basis function of quadrangle FEM element
		double dq2D_22(double x, double y) noexcept;

		//! x-derivative of the 3-th basis function of quadrangle FEM element
		double dq2D_31(double x, double y) noexcept;

		//! y-derivative of the 3-th basis function of quadrangle FEM element
		double dq2D_32(double x, double y) noexcept;

		//! x-derivative of the 4-th basis function of quadrangle FEM element
		double dq2D_41(double x, double y) noexcept;

		//! y-derivative of the 4-th basis function of quadrangle FEM element
		double dq2D_42(double x, double y) noexcept;

		// Triangle
		//! First basis function of triangle FEM element
		double t2D_1(double x, double y) noexcept;

		//! Second basis function of triangle FEM element
		double t2D_2(double x, double y) noexcept;

		//! 3-th basis function of triangle FEM element
		double t2D_3(double x, double y) noexcept;

		//! x-derivative of the first basis function of triangle FEM element
		constexpr double dt2D_11(double /*x*/, double /*y*/) noexcept{
			return -1.0;
		}

		//! y-derivative of the first basis function of triangle FEM element
		constexpr double dt2D_12(double /*x*/, double /*y*/) noexcept{
			return -1.0;
		}

		//! x-derivative of the second basis function of triangle FEM element
		constexpr double dt2D_21(double /*x*/, double /*y*/) noexcept{
			return 1.0;
		}

		//! y-derivative of the second basis function of triangle FEM element
		constexpr double dt2D_22(double /*x*/, double /*y*/) noexcept{
			return 0.0;
		}

		//! x-derivative of the 3-th basis function of triangle FEM element
		constexpr double dt2D_31(double /*x*/, double /*y*/) noexcept{
			return 0.0;
		}

		//! y-derivative of the 3-th basis function of triangle FEM element
		constexpr double dt2D_32(double /*x*/, double /*y*/) noexcept{
			return 1.0;
		}
	}; //namespace BaseElement
}; //namespace FEM
#endif
