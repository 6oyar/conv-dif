
#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <fstream>
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

#include "FEM/test/fem_test_utils.h"
#include "test.h"

double source(double x, double y) {
	return -256.0*((1.0 - y)*(1.0 - y)*y*y*(2.0 - 12.0*x + 12.0*x*x)
		+ (1.0 - x)*(1.0 - x)*x*x*(2.0 - 12.0*y + 12.0*y*y));
	//static const double A = 1e2;
	//static const double Q = 1e2;
	//double rho = (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5);
	//return Q*std::exp(-A*rho);
}
double solution(double x, double y) {
	return 256.0*(1.0 - x)*(1.0 - x)*x*x*(1.0 - y)*(1.0 - y)*y*y;
}


void outErrors(
	Eigen::SparseMatrix<double> &A, 
	Eigen::VectorXd &u, 
	Eigen::VectorXd &rhs, 
	Eigen::VectorXd &uExact) {
	printf("Errors: \n");
	printf(
		"\t||u-u^*|| = %le\n", 
		(u - uExact).lpNorm<Eigen::Infinity>());
	printf(
		"\t||Au-b||  = %le\n", 
		(A*u - rhs).lpNorm<Eigen::Infinity>());
	printf(
		"\t||Au^*-b||  = %le\n",
		(A*uExact - rhs).lpNorm<Eigen::Infinity>());
}
void solverLLTtest(
	Eigen::SparseMatrix<double> &A,
	Eigen::VectorXd &u, Eigen::VectorXd &rhs, 
	Eigen::VectorXd &uExact) {
	printf("--------- LLT solver ----------\n");
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower > solver;
	timer();
	solver.analyzePattern(A);
	//solver.compute(A);
	printf("\t1. Analyzed at %le s\n", timer());
	solver.factorize(A);
	printf("\t2. Factorized at %le s\n", timer());

	u = solver.solve(rhs);
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}
void solverLDLTtest(
	Eigen::SparseMatrix<double> &A,
	Eigen::VectorXd &u, 
	Eigen::VectorXd &rhs, 
	Eigen::VectorXd &uExact) {
	printf("--------- LDLT solver ----------\n");
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>,
		Eigen::Lower > solver;
	timer();
	solver.analyzePattern(A);
	printf("\t1. Analyzed at %le s\n", timer());
	solver.factorize(A);
	printf("\t2. Factorized at %le s\n", timer());

	u = solver.solve(rhs);
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}
void solverCGtest(Eigen::SparseMatrix<double> &A, 
	Eigen::VectorXd &u, 
	Eigen::VectorXd &rhs, 
	Eigen::VectorXd &uExact) {
	printf("--------- CG solver ----------\n");
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
		Eigen::Lower > solver;
	timer();
	solver.analyzePattern(A);
	printf("\t1. Analyzed at %le s\n", timer());
	solver.factorize(A);
	printf("\t2. Factorized at %le s\n", timer());
	solver.setMaxIterations(999);
	solver.setTolerance(1e-3);

	u = solver.solve(rhs);
	std::cout << solver.iterations() << "\n";
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}
void solverBiCGtestIlut(Eigen::SparseMatrix<double> &A, 
	Eigen::VectorXd &u, 
	Eigen::VectorXd &rhs, 
	Eigen::VectorXd &uExact) {
	printf("--------- BiCG+Ilut solver ----------\n");
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
		Eigen::IncompleteLUT<double>> solver;
	timer();
	solver.analyzePattern(A);
	printf("\t1. Analyzed at %le s\n", timer());
	solver.factorize(A);
	printf("\t2. Factorized at %le s\n", timer());
	solver.setMaxIterations(999);
	solver.setTolerance(1e-3);

	u = solver.solve(rhs);
	std::cout << solver.iterations() << "\n";
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}
void solverBiCGtest(Eigen::SparseMatrix<double> &A, 
	Eigen::VectorXd &u, 
	Eigen::VectorXd &rhs, 
	Eigen::VectorXd &uExact) {
	printf("--------- BiCG solver ----------\n");
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
	timer();
	solver.analyzePattern(A);
	printf("\t1. Analyzed at %le s\n", timer());
	solver.factorize(A);
	printf("\t2. Factorized at %le s\n", timer());
	solver.setMaxIterations(999);
	solver.setTolerance(1e-3);

	u = solver.solve(rhs);
	std::cout<<solver.iterations()<< "\n";
	printf("\t3. Solved at %le s\n", timer());
	outErrors(A, u, rhs, uExact);
}





void test1() {
	size_t N = 101;
	try {

		printf("Calculating for N = %d\n", N);
		FEM::SimpleMeshGenerator meshGen;
		FEM::CMesh mesh = meshGen.onRectangle(
			Eigen::Vector2d(0.0, 0.0), 
			Eigen::Vector2d(1.0, 1.0), 
			Eigen::Vector2i(N, N));

		Eigen::VectorXd u(mesh.Nn);
		Eigen::VectorXd rhs(mesh.Nn);
		Eigen::VectorXd uExact(mesh.Nn);
		auto nodes = mesh.getNodes();
		size_t i = 0;
		for (const auto& x : nodes) {
			rhs[i] = source(x[0], x[1]);
			uExact[i] = solution(x[0], x[1]);
			i++;
		}


		FEM::CBoundaryList bndCondList;
		bndCondList.addDirichlet(0, Eigen::Vector3d(-1, 0, 0));
		bndCondList.addDirichlet(0, Eigen::Vector3d(1, 0, 0));
		bndCondList.addDirichlet(0, Eigen::Vector3d(0, 1, 0));
		bndCondList.addDirichlet(0, Eigen::Vector3d(0, -1, 0));

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
		spMatPortrait(V, "conv.bmp");*/


		std::cout << "Dimenshion of A - " << A.cols() << "\n";
		std::cout << "Non zero elements - " << A.nonZeros() << "\n";
		printf("Assembling done at %.5f s\n", timer());

		solverLDLT(A, u, rhs, uExact);
		writeSolutionTecplot("u_ldlt.dat", "LDLT", u, mesh);
		u.setZero(mesh.Nn);

		solverLLT(A, u, rhs, uExact);
		writeSolutionTecplot("u_llt.dat", "LLT", u, mesh);
		u.setZero(mesh.Nn);

		solverCholesky(A, u, rhs, uExact);
		writeSolutionTecplot("u_chol.dat", "Cholesky", u, mesh);
		u.setZero(mesh.Nn);

		solverCG(A, u, rhs, uExact);
		writeSolutionTecplot("u_cg.dat", "CG", u, mesh);
		u.setZero(mesh.Nn);

		solverBiCGILUT(A, u, rhs, uExact);
		writeSolutionTecplot("u_bicgILUT.dat", "BiCG+ilut", u, mesh);
		u.setZero(mesh.Nn);

		solverBiCG(A, u, rhs, uExact);
		writeSolutionTecplot("u_bicg.dat", "BiCG", u, mesh);
		u.setZero(mesh.Nn);
	}
	catch (std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
		exit(EXIT_FAILURE);
	}

}
void test2() {
	const size_t N = 11;
	const double v = 0.01;
	const double tau = 0.1 / (v*(N - 1));
	const double T = 1000.0*tau;//
	const double a = 1e-3;
	try {
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
		for (const auto& x : nodes) {
			q[i] = tau*source(x[0], x[1]);
			i++;
		}


		timer();
		FEM::Assembler::boundaryConditions(
			Mb, 
			rhsb, 
			bndCondList, 
			mesh);
		Mb *= tau;
		rhsb *= tau;

		FEM::Assembler::div_k_grad(A, mesh, tau*a);
		FEM::Assembler::mass(M, mesh);

		FEM::Assembler::FVM(
			V, 
			mesh, 
			Eigen::Vector3d(v, 0.0, 0.0));

		V = tau*M.cwiseInverse().asDiagonal()*V;
		q = M.asDiagonal()*q;
		A += (Mb + M).asDiagonal();

		printf("Assembling done at %lf s\n", timer());



		double aa[121][121];
		for (size_t i = 0; i < 121; i++)
		{
			for (size_t j = 0; j <121; j++)
			{
				aa[i][j] = 0.0;
			}
		}

		{
			for (int k = 0; k < A.outerSize(); ++k) {
				for (Eigen::SparseMatrix<double>::
					InnerIterator it(A, k);
					it; 
					++it) {
					double v = it.value();
					int j = it.row();
					int i = it.col();
					aa[i][j] = v;
				}
			}
		}
		std::ofstream fout;
		fout.open("Matrix.txt");


		for (size_t i = 0; i < 121; i++)
		{
			for (size_t j = 0; j < 121; j++)
			{
				fout << aa[i][j]  << "\n";
			}
		}
		fout.close();

		fout.open("rhs.txt");

		Eigen::VectorXd rhs = q + rhsb + M.asDiagonal()*u;
		for (size_t i = 0; i < 121; i++)
		{
			fout << rhs(i) <<"\n";
		}
		fout.close();
		getchar();
		//spMatPortrait(V, "fem_test_3_conv.bmp");
		spMatPortrait(A, "fem_test_3_A.bmp");

		double time = tau;
		printf("\nTime = 0");
		char buff[100];
		int k = 0;
		while (time <= T) {
			u -= V*u; //FVM
			Eigen::VectorXd rhs = (
				q + rhsb + M.asDiagonal()*u);
			solveSLE(
				A,
				u,
				rhs
			);
			k++;
			if (k % 10 == 0) {
				long tm = time;
				sprintf(buff, 
					"%05ld.%05ld", 
					tm, 
					(long)((time-tm)*1e5));
				writeSolutionTecplot(buff, 
					"BiCG",
					u, 
					mesh);
				printf("\rTime = %lf", time);
				k = 0;
			}
			time += tau;
		}
		std::cout <<"\n"<< k << "\n";
		printf("\n");
	/*	char buff[100];
		sprintf(buff, "%f", time);
		writeSolutionTecplot(buff, "BiCG", u, mesh);*/
	}
	catch (std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
		exit(EXIT_FAILURE);
	}

}
void test3() {
	const int N = 201;
	const double v = 0.01;
	const double tau = 0.1 / (v*(N - 1));
	const double T = 1000.0*tau;//
	const double a = 1e-3;
	try {
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
		bndCondList.addGraneListByNormal(
			mesh.getBoundaryGranes());

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
		for (const auto& x : nodes) {
			q[i] = tau*source(x[0], x[1]);
			i++;
		}


		timer();
		FEM::Assembler::boundaryConditions(
			Mb,
			rhsb, 
			bndCondList, 
			mesh);
		Mb *= tau;
		rhsb *= tau;

		FEM::Assembler::div_k_grad(A, mesh, tau*a);
		FEM::Assembler::mass(M, mesh);

		FEM::Assembler::FVM(
			V, 
			mesh,
			Eigen::Vector3d(v, 0.0, 0.0));
	
		printf("Assembling done at %lf s\n", timer());
		V = tau*M.cwiseInverse().asDiagonal()*V;
		q = M.asDiagonal()*q;
		A += (Mb + M).asDiagonal();

		//spMatPortrait(V, "fem_test_3_conv.bmp");
		//spMatPortrait(A, "fem_test_3_A.bmp");
		printf("Size of matrix - %d\n", A.innerSize());
		printf("Non zeros - %d\n", A.nonZeros());
		double time = tau;
		//printf("\nTime = 0");
		char buff[100];
		while (time <= T) {
			u -= V*u; //FVM
			Eigen::VectorXd rhs = (
				q + rhsb + M.asDiagonal()*u);
			solveSLE(
				A,
				u,
				rhs
			);
				printf("\rTime = %lf", time);

			time += tau;
		}
	//	printf("\nBest time = %lf", bestTime);
		printf("\n");
	//	char buff[100];
		sprintf(buff, "%f", time);
	//	writeSolutionTecplot(buff, "BiCG", u, mesh);
	}
	catch (std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
		exit(EXIT_FAILURE);
	}

}
void testExact() {
	const size_t N = 501;
	try {
		FEM::SimpleMeshGenerator meshGen;
		FEM::CMesh mesh = meshGen.onRectangle(
			Eigen::Vector2d(0.0, 0.0), 
			Eigen::Vector2d(1.0, 1.0), 
			Eigen::Vector2i(N, N));

		Eigen::VectorXd u(mesh.Nn);
		Eigen::VectorXd rhs(mesh.Nn);
		Eigen::VectorXd uExact(mesh.Nn);
		auto nodes = mesh.getNodes();
		size_t i = 0;
		for (const auto& x : nodes) {
			rhs[i] = source(x[0], x[1]);
			uExact[i] = solution(x[0], x[1]);
			i++;
		}


		FEM::CBoundaryList bndCondList;
		bndCondList.addDirichlet(
			0, 
			Eigen::Vector3d(-1, 0, 0));
		bndCondList.addDirichlet(
			0, 
			Eigen::Vector3d(1, 0, 0));
		bndCondList.addDirichlet(
			0, 
			Eigen::Vector3d(0, 1, 0));
		bndCondList.addDirichlet(
			0,
			Eigen::Vector3d(0, -1, 0));

		bndCondList.addGraneListByNormal(
			mesh.getBoundaryGranes());
		Eigen::VectorXd Mb(mesh.Nn);
		Eigen::VectorXd rhsb(mesh.Nn);
		timer();
		FEM::Assembler::boundaryConditions(
			Mb, 
			rhsb, 
			bndCondList, 
			mesh);

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



		writeSolutionTecplot("u_exact.dat", "Exact", uExact, mesh);

		u.setZero(mesh.Nn);
		solverBiCGtest(A, u, rhs, uExact);
		writeSolutionTecplot("u_bicg.dat", "BiCG", u, mesh);
		u.setZero(mesh.Nn);

		solverBiCGtestIlut(A, u, rhs, uExact);
		writeSolutionTecplot("u_bicg_ilut.dat", "BiCG", u, mesh);
		u.setZero(mesh.Nn);

		
		solverCGtest(A, u, rhs, uExact);
		writeSolutionTecplot("u_cg.dat", "BiCG", u, mesh);
		u.setZero(mesh.Nn);

		solverLLTtest(A, u, rhs, uExact);
		writeSolutionTecplot("u_llt.dat", "LLT", u, mesh);
		u.setZero(mesh.Nn);

		solverLDLTtest(A, u, rhs, uExact);
		writeSolutionTecplot("u_ldlt.dat", "LLT", u, mesh);
		u.setZero(mesh.Nn);



		
		
		
	}
	catch (std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
		exit(EXIT_FAILURE);
	}

}
int main(void) {
	//test1();
	//test2();
	//test3();
	//testExact();
	system("pause");
	return EXIT_SUCCESS;
}
