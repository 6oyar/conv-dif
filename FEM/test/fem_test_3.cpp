/* test FEM+FVM for du/dt - div k grad u + v grad u = f */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <exception>
#include <functional>
#include <chrono>

#include <image/bitmap_image.hpp>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <FEM/fem_element.h>
#include <FEM/fem_mesh.h>
#include <FEM/fem_boundary.h>
#include <FEM/fem_assembler.h>

#include "fem_test_utils.h"

double source(double x, double y){
	static const double A = 1e2;
	static const double Q = 1e2;
	double rho = (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5);
	return Q*std::exp(-A*rho);
}

void outErrors(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact){
	printf("Errors: \n");
	printf("\t||Au-b||  = %le\n", (A*u-rhs).lpNorm<Eigen::Infinity>());
}

inline void solverBiCG(
	Eigen::SparseMatrix<double> &A, 
	Eigen::VectorXd &u, 
	Eigen::VectorXd &rhs,
	bool refactorize = false)
{
	static bool isInitialized = false;
	static Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
	if (!isInitialized || refactorize){
		if (!isInitialized){
			solver.analyzePattern(A);
		}
		solver.factorize(A);
		solver.setTolerance(1e-3);

		isInitialized = true;
	}
	
	u = solver.solve(rhs);
}

inline double SORIteration(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, double omega){
	double uError = 0.0;
	double *uData = u.data();
	for (int i=0; i<A.outerSize(); ++i){
		double &ud = uData[i];
		double f = rhs[i];
		double d = 1.0;
 		for (Eigen::SparseMatrix<double>::InnerIterator it(A,i); it; ++it){
    		const size_t j = it.row();
    		if (j != i){
    			f -= it.value()*uData[j];
    		}
    		else{
    			d = it.value();
    		}
    	}
    	const double ui = f/d;
    	double ue = omega*(ui - ud);

    	ud += ue;
		ue = std::abs(ue);
		if (ue > uError) uError = ue;
  	}

  	return uError;
}

void solverSOR(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, double omega){
	double uError = 100.0;
	int it = 0;
	while(uError > 1e-8 && it < 2000){
		uError = SORIteration(A, u, rhs, omega);
	}
}

int main(void){
	const size_t N   = 101;
	const double v   = 0.01;
	const double tau = 0.1/(v*(N-1));
	const double T   = 1000.0*tau;//30.0;
	const double a   = 1e-3;
    try{
    	FEM::SimpleMeshGenerator meshGen;
    	FEM::CMesh mesh = meshGen.onRectangle(
    		Eigen::Vector2d(0.0, 0.0), 
    		Eigen::Vector2d(1.0, 1.0), 
    		Eigen::Vector2i(N, N)
    	);

    	FEM::CBoundaryList bndCondList;
    	bndCondList.addDirichlet(
    		0, 
    		Eigen::Vector3d(-1.0, 0.0, 0.0)
    	);
    	bndCondList.addGraneListByNormal(mesh.getBoundaryGranes());

    	Eigen::VectorXd u(mesh.Nn); u.setZero();
    	Eigen::VectorXd q(mesh.Nn);

    	Eigen::SparseMatrix<double> A(mesh.Nn, mesh.Nn);
    	Eigen::SparseMatrix<double> V(mesh.Nn, mesh.Nn);
    	//Eigen::SparseMatrix<double> G(mesh.Nn, mesh.Nn);

    	Eigen::VectorXd M(mesh.Nn);
    	Eigen::VectorXd Mb(mesh.Nn);
		Eigen::VectorXd rhsb(mesh.Nn);


    	auto nodes = mesh.getNodes();
    	size_t i = 0;
    	for(const auto& x: nodes){
    		q[i] = tau*source(x[0], x[1]);
    		i++;
    	}

    	
		timer();
		FEM::Assembler::boundaryConditions(Mb, rhsb, bndCondList, mesh);
		Mb   *= tau;
		rhsb *= tau;

    	FEM::Assembler::div_k_grad(A, mesh, tau*a);
    	FEM::Assembler::mass(M, mesh);

    	FEM::Assembler::FVM(V, mesh, Eigen::Vector3d(v, 0.0, 0.0));
    	
    	V = tau*M.cwiseInverse().asDiagonal()*V;
    	q = M.asDiagonal()*q;
    	A += (Mb + M).asDiagonal();

		printf("Assembling done at %le s\n", timer());

		//spMatPortrait(V, "fem_test_3_conv.bmp");
		spMatPortrait(A, "fem_test_3_A.bmp");

		double time = tau;
		printf("\nTime = 0");
		
		while(time <= T){
			u -= V*u; //FVM
			Eigen::VectorXd rhs = (q + rhsb + M.asDiagonal()*u);
			solverBiCG(
				A,
				u,
				rhs
			);
			printf("\rTime = %le", time);
			time += tau;
		}
		printf("\n");

		writeSolutionTecplot("fem_test_3.dat", "BiCG", u, mesh);
    }
    catch(std::exception &e){
    	fprintf(stderr, "%s\n", e.what());
    	exit(EXIT_FAILURE);
    } 
	
	return EXIT_SUCCESS;
}
