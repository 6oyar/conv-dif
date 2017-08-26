#ifndef __KRF_FEM_LA_H__
#define __KRF_FEM_LA_H__
#include <memory.h>
#include <cmath>
#include "noexcept_hack.h"
namespace FEM{
	namespace LA{
		template<typename T>
		T det(T const* const* const a, int N) noexcept{
			if (N < 1) return 0.0;
			switch(N){
				case 3: return   a[1][1]*(a[0][0]*a[2][2] - a[2][0]*a[0][2]) 
							   + a[2][1]*(a[1][0]*a[0][2] - a[0][0]*a[1][2])
							   + a[0][1]*(a[2][0]*a[1][2] - a[2][2]*a[1][0]);

				case 2: return a[0][0]*a[1][1]-a[0][1]*a[1][0];

				case 1: return a[0][0];
			}

			T  ret = 0.0;
			T  mul = 1.0;
			const T **a1 = (const T**)alloca(sizeof(T)*(N-1));

			for (int i=N-1; i > -1; --i){
				int k = 0;
				for (int j=0;   j<i; ++j) a1[k++] = a[j];
				for (int j=i+1; j<N; ++j) a1[k++] = a[j];
				ret += mul*a[i][N-1]*det(a1, N-1);
				mul = -mul;
			}

			return ret;
		}

		template<typename T>
		bool solveAx_3x3(T const *const *const A, T *x, const T *b) noexcept{
			const size_t sz = sizeof(T)*3;
			static T Ai[3][3];
			const T dA = det(A, 3);
			
			if (std::abs(dA) < 1e-15) return false;

			for (size_t k=0; k<3; ++k){
				for (size_t i=0; i<3; ++i){
					memcpy((void *)(Ai[i]), (void *)(A[i]), sz);
					Ai[i][k] = b[i];
				}
				
				x[k] = det(Ai, 3)/dA;
			}

			return true;
		}
	}; //namespace LA
}; //namespace FEM


#endif