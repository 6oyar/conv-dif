#include <stdarg.h>
#include <assert.h>
#include <cstdlib>
#include <cstdio>

#include "fem_element.h"
#include "fem_la.h"

namespace FEM{
	using Vector3d   = Eigen::Vector3d;
	using AsVector3d = Eigen::Map<Vector3d>;

	CElement::CElement(TElementType et, size_t eNum, bool pIsExternalNodesStorage) noexcept: 
		num(eNum),
	    type(et), N(BaseElement::elementNodesCount(et)), dim(BaseElement::elementDim(et)), 
	    isExternalNodesStorage(pIsExternalNodesStorage)
	{
		double *memPtr; double **memPointersPtr;

		mem = new double[
			((isExternalNodesStorage) ? 0 : N*dim) + //node
			N         + //me
			N*N       + //ae
			N*N*dim   + //de
			N*dim*dim + //J
			N*dim*dim + //J1
			N*dim*N   + //JD
			N           //detJ
		];

		memPointers = new double*[
			N           + //node
			N           + //ae
			dim + N*dim + //de
			N   + N*dim + //J
			N   + N*dim + //J1
			N   + N*dim   //JD
		];

		memPtr         = mem;
		memPointersPtr = memPointers;

		node = (double **)memPointersPtr; memPointersPtr += N;
		if (!isExternalNodesStorage){
			for (int i=0; i<N; ++i){
				node[i] = memPtr; memPtr += dim;
			}
		}

		me = memPtr; memPtr += N;
		ae = (double **)memPointersPtr; memPointersPtr += N;
	    for (int i=0; i<N; ++i){
			ae[i] = memPtr; memPtr += N;
		}

		de = (double ***)memPointersPtr; memPointersPtr += dim;
		for (int i=0; i<dim; ++i){
			de[i] = (double **)memPointersPtr; memPointersPtr += N;
			for (int j=0; j<N; ++j){
				de[i][j] = memPtr; memPtr += N;
			}
		}

		J = (double ***)memPointersPtr; memPointersPtr += N;
		for (int i=0; i<N; ++i){
			J[i] = (double **)memPointersPtr; memPointersPtr += dim;
			for (int j=0; j<dim; ++j){
				J[i][j] = memPtr; memPtr += dim;
			}
		}

		J1 = (double ***)memPointersPtr; memPointersPtr += N;
		for (int i=0; i<N; ++i){
			J1[i] = (double **)memPointersPtr; memPointersPtr += dim;
			for (int j=0; j<dim; ++j){
				J1[i][j] = memPtr; memPtr += dim;
			}
		}

		JD = (double ***)memPointersPtr; memPointersPtr += N;
		for (int i=0; i<N; ++i){
			JD[i] = (double **)memPointersPtr; memPointersPtr += dim;
			for (int j=0; j<dim; ++j){
				JD[i][j] = memPtr; memPtr += N;
			}
		}

		detJ = memPtr;
        elNode = BaseElement::elementNodes(et);
        BaseElement::elementBasis(et, phi, dphi);

        _isAssembled = false;
	}

	void CElement::calculateBaricenter() noexcept{
		double *b = baricenter.data();
		memset(b, 0, sizeof(double)*3);

		for(int i=0; i<N; ++i){
			for(int j=0; j<dim; ++j){
				b[j] += node[i][j];
			}
		}

		for(int i=0; i<dim; ++i){
			b[i] /= N;
		}
	}

	void CElement::transformBaseVector(double *p) const noexcept{
		double *x = (double *)alloca(sizeof(double)*dim);
		memset(x, 0, sizeof(double)*dim);
		
		for (int i=0; i<N; ++i){
			double fi = phi(i, p);
			for(int j=0; j<dim; ++j){
				x[j] += node[i][j]*fi;
			}
		}

		for (int i=0; i<dim; ++i){
			x[i] -= node[0][i];
		}

		memcpy(p, x, sizeof(double)*dim);
	}

	Vector3d CElement::_normalToGrane(TElementType et, const size_t *nodeNums) const noexcept{
		Vector3d ret(0, 0, 0);
		switch(et){
			case etL:{
				AsVector3d A(node[nodeNums[0]], 3);
				AsVector3d B(node[nodeNums[1]], 3);
				AsVector3d C(node[nodeNums[2]], 3);

				auto AB = B-A;
				return (C-A).cross(AB).cross(AB).normalized();
			}
			case etT:
			case etQ:
			default:{
				fprintf(stderr, "Undefined normal for triangle and quadrangle!");
				exit(1);
			}
		}

		return ret;
	}

	Vector3d CElement::normalToGrane(int ng) const noexcept{
		size_t nodeNums[5];
		const TElementType et = BaseElement::elementGrane(type, ng, nodeNums);

		return _normalToGrane(et, nodeNums);
	}

	CElement CElement::graneElement(int ng, Vector3d &normal, size_t *nodeNums) const noexcept{
		const TElementType et = BaseElement::elementGrane(type, ng, nodeNums);
		normal = _normalToGrane(et, nodeNums);
		switch(et){
			case etL:{
				AsVector3d A(node[nodeNums[0]], 3);
				AsVector3d B(node[nodeNums[1]], 3);

				double norm = (B-A).norm();
				CElement ret(et);
				ret.setNode(0,  0.0, 0.0);
				ret.setNode(1, norm, 0.0);
				Integrator::localMatrix(ret);

				return ret;
			}
			default:{
				fprintf(stderr, "Granes for triangle and quadrangle is undefined");
				exit(1);
			}
		}
	}

	bool CElement::setNode(int i, ...) noexcept{
		assert(i < N && !isExternalNodesStorage);
		int k;
		va_list ap;
		va_start(ap, i);
		for (k=0; k<dim; k++) {
			node[i][k] = va_arg(ap, double);
		}
		va_end(ap);

		_isAssembled = false;
		return (bool)(k == dim);
	}

	bool CElement::setNode(int i, double *nodeCoordinates) noexcept{
		assert(i < N);
		node[i] = nodeCoordinates;
		_isAssembled = false;
		return true;
	}
    
    void CElement::printMatrix() const noexcept{
        for(int i=0; i<N; ++i){
            for(int j=0; j<N; ++j){
                printf("%2.5lf ", ae[i][j]);
            }
            printf("| %2.5lf\n", me[i]);
        }

        for(int d=0; d<dim; ++d){
        	for(int i=0; i<N; ++i){
            	for(int j=0; j<N; ++j){
                	printf("%2.5lf ", de[d][i][j]);
            	}	
            	printf("\n");
        	}
        	printf("\n");
        }
    }

	void Integrator::_calcJFull(CElement &e) noexcept{
		const int N     = e.N;
		const TDim dim   = e.dim;
		double **node   = e.node;
		auto elNode = e.elNode;
		t_dphi dphi     = e.dphi;
		
		double ***J  = e.J;
		double ***J1 = e.J1;
		double ***JD = e.JD;
		double *detJ = e.detJ;

		double *dphiTmp = (double*)alloca(sizeof(double)*N*dim);
        
		memset( J[0][0], 0, N*dim*dim*sizeof(double));
		memset(JD[0][0], 0, N*dim*N*sizeof(double));
		for (int k=0; k<N; k++){
			for(int i=0; i<dim; ++i){
				for (int l=0; l<N; ++l){
					double dphi_lik = dphiTmp[l*dim+i] = (*dphi)(l,i,elNode[k]);
					for(int j=0; j<dim; ++j){
						J[k][i][j] += dphi_lik*node[l][j];
					}
				}
			}
			detJ[k] = LA::det(J[k], dim);
			switch (dim){
				case dim1D:{
					J1[k][0][0] = 1.0/J[k][0][0];
					break;
				}

				case dim2D:{
					J1[k][0][0] =  J[k][1][1]/detJ[k];
					J1[k][0][1] = -J[k][0][1]/detJ[k];
					J1[k][1][0] = -J[k][1][0]/detJ[k];
					J1[k][1][1] =  J[k][0][0]/detJ[k];
					break;	
				}

				case dim3D:{
					J1[k][0][0] =  (J[k][1][1]*J[k][2][2] - J[k][1][2]*J[k][2][1])/detJ[k];
					J1[k][0][1] = -(J[k][0][1]*J[k][2][2] - J[k][0][2]*J[k][2][1])/detJ[k];
					J1[k][0][2] =  (J[k][0][1]*J[k][1][2] - J[k][0][2]*J[k][1][1])/detJ[k];

					J1[k][1][0] = -(J[k][1][0]*J[k][2][2] - J[k][1][2]*J[k][2][0])/detJ[k];
					J1[k][1][1] =  (J[k][0][0]*J[k][2][2] - J[k][0][2]*J[k][2][0])/detJ[k];
					J1[k][1][2] = -(J[k][0][0]*J[k][1][2] - J[k][0][2]*J[k][1][0])/detJ[k];

					J1[k][2][0] =  (J[k][1][0]*J[k][2][1] - J[k][1][1]*J[k][2][0])/detJ[k];
					J1[k][2][1] = -(J[k][0][0]*J[k][2][1] - J[k][0][1]*J[k][2][0])/detJ[k];
					J1[k][2][2] =  (J[k][0][0]*J[k][1][1] - J[k][0][1]*J[k][1][0])/detJ[k];
					break;
				}
                default:
                    fprintf(stderr, "Unknown dimension in Integrator::_calcJFull\n");
                    exit(1);
                
			}
			for (int i=0; i<dim; ++i){
				for (int l=0; l<N; ++l){
					for (int j=0; j<dim; ++j){
						JD[k][i][l] += J1[k][i][j]*dphiTmp[l*dim+j];
					}
				}
			}
		}

		e._isAssembled = true;
	}

	void Integrator::_calcJMass(CElement &e) noexcept{
		const int N = e.N;
		const int dim = e.dim;
		double **node = e.node;
		const auto elNode = e.elNode;
		t_dphi dphi = e.dphi;
		
		double ***J  = e.J;
		double *detJ = e.detJ;

		memset(J[0][0], 0, N*dim*dim*sizeof(double));
		for (int k=0; k<N; ++k){
			for(int i=0; i<dim; ++i){
				for (int l=0; l<N; ++l){
					const double dphi_lik = (*dphi)(l,i,elNode[k]);
					for(int j=0; j<dim; ++j){
						J[k][i][j] += dphi_lik*node[l][j];
					}
				}
			}
			detJ[k] = LA::det(J[k], dim);
		}
	}

	void Integrator::localMatrix(CElement &e) noexcept{
		if (e.isAssembled()) return;

		_calcJFull(e);

		const int N   = e.N;
		const int dim = e.dim;
		double ***JD  = e.JD;
		double *detJ  = e.detJ;
		
		double Sk = BaseElement::elementMeasure(e.type)/N;
		double *S = (double *)alloca(sizeof(double)*N);
		
		memcpy(e.me, detJ, N*sizeof(double));
		for (int i=0; i<N; ++i){
			e.me[i] *= Sk;
		}

		memset(e.ae[0], 0, N*N*sizeof(double));
		for (int i=0; i<N; ++i){
			for (int j=0; j<N; ++j){
				memset(S, 0, sizeof(double)*N);
				for(int k=0; k<N; ++k){
					for (int l=0; l<dim; ++l){
						S[k] += JD[k][l][i]*JD[k][l][j];
					}
				}
				for(int k=0; k<N; ++k){
					S[k] *= e.me[k];
				}

				for(int k=0; k<N; ++k){
					e.ae[i][j] += S[k];
				}
			}
		}


		memset(e.de[0][0], 0, N*N*dim*sizeof(double));
		for (int k=0; k<dim; ++k){
			for (int i=0; i<N; ++i){
				for (int j=0; j<N; ++j){
					//DJ[i][k][j] = (J^{-1}*grad{\\phi_j})(a_i)_k
					e.de[k][i][j] = JD[i][k][j]*e.me[i];
				}
			}
		}
		e.calculateBaricenter();
	}

	void Integrator::massMatrix(CElement &e) noexcept{
		if (e.isAssembled()) return;
		_calcJMass(e);

		double Sk = BaseElement::elementMeasure(e.type)/e.N;
		memcpy(e.me, e.detJ, sizeof(double)*e.N);
		for (int i=0; i<e.N; ++i){
			e.me[i] *= Sk;
		}
	}
}; //namespace FEM
