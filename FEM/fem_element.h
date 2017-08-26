#ifndef __KRF_LFEM_H__
#define __KRF_LFEM_H__
#include <cstdlib>
#include "fem_base_element.h"
#include "noexcept_hack.h"
namespace FEM{
	//! FEM element class
	class CElement{
	private:
		double *mem;          //!< memory for all matrices
		double **memPointers; //!< memory for all pointers

		bool _isAssembled;
		friend class Integrator;
		friend class Assembler;
		Eigen::Vector3d _normalToGrane(TElementType et, const size_t *nodeNums) const noexcept;

		Eigen::Vector3d baricenter;
		void calculateBaricenter() noexcept;
	public:
		const size_t num;        //!< element number
		const TElementType type; //!< element type (triangle, quadrangle, etc)
		const int N;             //!< number of nodes
		const TDim dim;          //!< dimension
		const bool isExternalNodesStorage; //!< is coordinates of element nodes is stored in external memory

		t_phi  phi;      //!< element basis functions
		t_dphi dphi;     //!< partial derivatives of the basis functions

		double *me;      //!< local mass matrix
		double **ae;     //!< local stiffnes matrix
		double ***de;    //!< local gradient matrices
		double ***J;     //!< transform Jacobi matrix in all nodes
		double ***J1;    //!< inversed Jacobi matrixes
		double *detJ;    //!< transform Jacobians in all nodes
		double ***JD;    //!< DJ[k][i][j] = (J^{-1}*grad{\\phi_j})(a_k)_i
		double **node;   //!< coordinates of nodes of initial elemet
		const double * const*elNode; //!< coordinates of basis element nodes

		//! Set coordinates of i-th node of element
		/*!
			\param i - number of elements node
			\param ... - list of \\x_1, \\x_2, ..., \\x_N coordinates 
			\return true if success and false othrwise
		*/
		bool setNode(int i, ...) noexcept;
		bool setNode(int i, double *nodeCoordinates) noexcept;
        
        void printMatrix() const noexcept;

        bool isAssembled() const noexcept{
        	return _isAssembled;
        };

        //!< p = \sum (\phi_i(p) x_i - \phi_i(0) x_i)
        void transformBaseVector(double *p) const noexcept;

        //!< outward ort normal to grane of element
        /*!
        	\param[in] ng - number of grane 
        	\return outward ort normal
        */
        Eigen::Vector3d normalToGrane(int ng) const noexcept;

	    //! Constructor
	    /*!
	    	\param et - element type
	    */
		CElement(TElementType et, size_t eNum = 0, bool pIsExternalNodesStorage = false) noexcept;
		CElement graneElement(int ng, Eigen::Vector3d &normal, size_t *nodeNums) const noexcept;

		CElement(const CElement&) = delete;
		CElement& operator=(const CElement&) = delete;
		CElement& operator=(CElement &&other) = delete;

		CElement(CElement &&other) noexcept:
			mem(other.mem), memPointers(other.memPointers),
			_isAssembled(other._isAssembled),
			num(other.num), type(other.type), N(other.N), dim(other.dim), isExternalNodesStorage(other.isExternalNodesStorage),
			phi(other.phi), dphi(other.dphi),
			me(other.me), ae(other.ae), de(other.de), J(other.J), J1(other.J1), detJ(other.detJ), JD(other.JD), 
			node(other.node), elNode(other.elNode)
		{
			other.mem = NULL;
			other.memPointers = NULL;
		}

		//! Destructor
		~CElement(){
			if (mem){
				delete[] mem;
			}

			if (memPointers){
				delete[] memPointers;
			}
		}
	}; //class CElement

	//! Calculation local matrixes of element
	class Integrator{
	private:
		//! Calculate Jacoby, inversed Jacoby matrixes and Jacobians of element
		/*!
			\param e - pointer to element
		*/
		static void _calcJFull(CElement &e) noexcept;

		//! Calculate Jacobians of element
		/*!
			\param e - pointer to element
		*/
		static void _calcJMass(CElement &e) noexcept;
	public:
		//! Calculate local mass and stiffnes matrixes of element
		/*!
			\param e - pointer to element
		*/
		static void localMatrix(CElement &e) noexcept;


		//! Calculate local mass matrix of element
		/*!
			\param e - pointer to element
		*/
		static void massMatrix(CElement &e) noexcept;
	}; //class Integrator
}; //namespace FEM
#endif
