#pragma once
#ifndef _TEST_
#define _TEST_
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include "fem_test_utils.h"
#include <iostream>

static double bestTime = std::numeric_limits<double>::max();


void solverLLT(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact) {
	printf("--------- LLT solver ----------\n");
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower > solver;
	double time1, time2;
	timer();
	solver.compute(A);
	time1 = timer();
	if (solver.info() != 0) {
		std::cout << "There was a failure." << "\n";
	}
	else {
		printf("\t1. Computed at %lf s\n", time1);
		u = solver.solve(rhs);
		time2 = timer();
		printf("\t2. Solved at %lf s\n", time2);
		printf("\t3. Total time %lf s\n", time1 + time2);
	}
}
void solverCholesky(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact) {
	printf("--------- Cholesky solver ----------\n");
	Eigen::SimplicialCholesky< Eigen::SparseMatrix<double>, Eigen::Lower > solver;
	double time1, time2;
	timer();
	solver.compute(A);
	time1 = timer();
	if (solver.info() != 0) {
		std::cout << "There was a failure." << "\n";
	}
	else {
		printf("\t1. Computed at %lf s\n", time1);
		u = solver.solve(rhs);
		time2 = timer();
		printf("\t2. Solved at %lf s\n", time2);
		printf("\t3. Total time %lf s\n", time1 + time2);
	}

}
void solverCG(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact) {
	printf("--------- CG solver ----------\n");
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
	double time1, time2;
	timer();
	solver.compute(A);
	solver.setTolerance(1e-3);
	time1 = timer();
	if (solver.info() != 0) {
		std::cout << "There was a failure." << "\n";
	}
	else {
		printf("\t1. Computed at %lf s\n", time1);
		u = solver.solve(rhs);
		time2 = timer();
		printf("\t2. Solved at %lf s\n", time2);
		printf("\t3. Total time %lf s\n", time1 + time2);
	}
}
void solverBiCGILUT(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact) {
	printf("--------- BiCG + ILUT solver ----------\n");
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
	double time1, time2;
	timer();
	solver.analyzePattern(A);
	solver.factorize(A);
	solver.setTolerance(1e-3);
	time1 = timer();
	if (solver.info() != 0) {
		std::cout << "There was a failure." << "\n";
	}
	else {
		printf("\t1. Computed at %lf s\n", time1);
		u = solver.solve(rhs);
		time2 = timer();
		printf("\t2. Solved at %lf s\n", time2);
		printf("\t3. Total time %lf s\n", time1 + time2);
	}
}
void solverBiCG(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact) {
	printf("--------- BiCG + Jacobi preconditioner solver ----------\n");
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver; // по умолчанию диагональный предобуславливатель
	double time1, time2;
	timer();
	solver.analyzePattern(A);
	solver.factorize(A);
	solver.setTolerance(1e-3);
	time1 = timer();
	if (solver.info() != Eigen::Success) {
		std::cout << "There was a failure." << "\n";
	}
	else {
		printf("\t1. Computed at %lf s\n", time1);
		u = solver.solve(rhs);
		time2 = timer();
		printf("\t2. Solved at %lf s\n", time2);
		printf("\t3. Total time %lf s\n", time1 + time2);
	}
}
void solverLDLT(Eigen::SparseMatrix<double> &A, Eigen::VectorXd &u, Eigen::VectorXd &rhs, Eigen::VectorXd &uExact) {
	printf("--------- LDLT solver ----------\n");
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower > solver;
	double time1, time2;
	timer();
	solver.compute(A);
	time1 = timer();
	if (solver.info() != 0) {
		std::cout << "There was a failure." << "\n";
	}
	else {
		printf("\t1. Computed at %lf s\n", time1);
		u = solver.solve(rhs);
		time2 = timer();
		printf("\t2. Solved at %lf s\n", time2);
		printf("\t3. Total time %lf s\n", time1 + time2);
	}
}


//void solveSLE(
//	Eigen::SparseMatrix<double> &A,
//	Eigen::VectorXd &u,
//	Eigen::VectorXd &rhs,
//	bool refactorize = false) {
//	static bool isInitialized = false;
//	static Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
//	if (!isInitialized || refactorize) {
//		if (!isInitialized) {
//			solver.analyzePattern(A);
//		}
//		solver.factorize(A);
//		isInitialized = true;
//	}
//
//	u = solver.solve(rhs);
//}
void solveSLE(
	Eigen::SparseMatrix<double> &A,
	Eigen::VectorXd &u,
	Eigen::VectorXd &rhs,
	bool refactorize = false) {
	
	static bool isInitialized = false;
	static Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,Eigen::Lower> solver;
	solver.setTolerance(1e-3);
	if (!isInitialized || refactorize) {
		timer();
		if (!isInitialized) {
			solver.analyzePattern(A);
		}
		solver.factorize(A);
		printf("Computed at %lf s\n", timer());
		isInitialized = true;
	}
	timer();
	u = solver.solve(rhs);
	double curTime = timer();
	if (curTime < bestTime) { bestTime = curTime; 
	}
}

#endif
