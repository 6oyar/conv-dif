/* test Eigen solvers for Laplass */
#include <iostream>
#include <cstdio>
#include <cstdlib>
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

double solution(double x, double y){
	return 256.0*(1.0-x)*(1.0-x)*x*x*(1.0-y)*(1.0-y)*y*y;
}

double source(double x, double y){
	return -256.0*( (1.0-y)*(1.0-y)*y*y*(2.0 - 12.0*x + 12.0*x*x)
	               +(1.0-x)*(1.0-x)*x*x*(2.0 - 12.0*y + 12.0*y*y));
}

void outErrors(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact){
	printf("Errors: \n");
	printf("\t||u-u^*|| = %le\n", (u-uExact).lpNorm<Eigen::Infinity>());
	printf("\t||Au-b||  = %le\n", (A*u-rhs).lpNorm<Eigen::Infinity>());
	printf("\t||Au^*-b||  = %le\n", (A*uExact-rhs).lpNorm<Eigen::Infinity>());
}

void solverBiCG(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact){
	printf("--------- BiCG solver ----------\n");
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
	timer();
	solver.analyzePattern(A);
	printf("\t1. Analyzed at %le s\n", timer());
	solver.factorize(A);
	printf("\t2. Factorized at %le s\n", timer());
	solver.setTolerance(1e-3);
	
	u = solver.solve(rhs);
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}

void solverLLT(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact){
	printf("--------- LLT solver ----------\n");
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
	timer();
	//solver.analyzePattern(A);
	solver.compute(A);
	printf("\t1. Analyzed at %le s\n", timer());
	//solver.factorize(A);
	//printf("\t2. Factorized at %le s\n", timer());
	
	u = solver.solve(rhs);
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}

void solverLDLT(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact){
	printf("--------- LDLT solver ----------\n");
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
	timer();
	solver.analyzePattern(A);
	printf("\t1. Analyzed at %le s\n", timer());
	solver.factorize(A);
	printf("\t2. Factorized at %le s\n", timer());
	
	u = solver.solve(rhs);
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}

void solverCholesky(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact){
	printf("--------- Cholesky solver ----------\n");
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
	timer();
	solver.analyzePattern(A);
	printf("\t1. Analyzed at %le s\n", timer());
	solver.factorize(A);
	printf("\t2. Factorized at %le s\n", timer());
	
	u = solver.solve(rhs);
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
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

void solverSOR(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, double omega, Eigen::VectorXd &uExact){
	printf("--------- SOR solver ----------\n");
	timer();
	double uError = 100.0;
	int it = 0;
	while(uError > 1e-8 && it < 2000){
		uError = SORIteration(A, u, rhs, omega);
	}
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}

int main(void){
	const size_t N = 101;
    try{
    	FEM::SimpleMeshGenerator meshGen;
    	FEM::CMesh mesh = meshGen.onRectangle(Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 1.0), Eigen::Vector2i(N, N));

    	Eigen::VectorXd u(mesh.Nn);
    	Eigen::VectorXd rhs(mesh.Nn);
    	Eigen::VectorXd uExact(mesh.Nn);
    	auto nodes = mesh.getNodes();
    	size_t i = 0;
    	for(const auto& x: nodes){
    		rhs[i] = source(x[0], x[1]);
    		uExact[i] = solution(x[0], x[1]);
    		i++;
    	}
    	

    	FEM::CBoundaryList bndCondList;
    	bndCondList.addDirichlet(0, Eigen::Vector3d(-1,  0, 0));
    	bndCondList.addDirichlet(0, Eigen::Vector3d( 1,  0, 0));
    	bndCondList.addDirichlet(0, Eigen::Vector3d( 0,  1, 0));
    	bndCondList.addDirichlet(0, Eigen::Vector3d( 0, -1, 0));

    	bndCondList.addGraneListByNormal(mesh.getBoundaryGranes());
		Eigen::VectorXd Mb(mesh.Nn);
		Eigen::VectorXd rhsb(mesh.Nn);
		timer();
		FEM::Assembler::boundaryConditions(Mb, rhsb, bndCondList, mesh);

    	Eigen::SparseMatrix<double> A(mesh.Nn, mesh.Nn);
    	FEM::Assembler::div_k_grad(A, mesh);

    	Eigen::VectorXd M(mesh.Nn);
    	FEM::Assembler::mass(M, mesh);

    	A += Mb.asDiagonal();
		rhs = M.asDiagonal()*rhs + rhsb;

		/*Eigen::SparseMatrix<double> V(mesh.Nn, mesh.Nn);
		FEM::Assembler::FVM(V, mesh, Eigen::Vector3d(1.0, 0.0, 0.0));
		spMatPortrait(V, "conv.bmp");
		std::cout << V << std::endl;*/

		printf("Assembling done at %le s\n", timer());
		
		solverLDLT(A, u, rhs, uExact); 
		writeSolutionTecplot("u_ldlt.dat", "LDLT", u, mesh); 
		u.setZero(mesh.Nn);
		
		solverLLT(A, u, rhs, uExact); 
		writeSolutionTecplot("u_llt.dat", "LLT", u, mesh); 
		u.setZero(mesh.Nn);

		solverCholesky(A, u, rhs, uExact);
		writeSolutionTecplot("u_chol.dat", "Colesky", u, mesh); 
		u.setZero(mesh.Nn);

		solverBiCG(A, u, rhs, uExact);
		writeSolutionTecplot("u_bicg.dat", "BiCG", u, mesh);
		u.setZero(mesh.Nn);


		solverSOR(A, u, rhs, 1.8, uExact); 
		writeSolutionTecplot("u_sor.dat", "SOR", u, mesh);
		u.setZero(mesh.Nn);
    }
    catch(std::exception &e){
    	fprintf(stderr, "%s\n", e.what());
    	exit(EXIT_FAILURE);
    }
	
	return EXIT_SUCCESS;
}
