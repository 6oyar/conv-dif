/* test FEM local matrices assembling */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <functional>

#include <image/bitmap_image.hpp>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include <FEM/fem_element.h>
#include <FEM/fem_mesh.h>
#include <FEM/fem_boundary.h>
#include <FEM/fem_assembler.h>

#include "fem_test_utils.h"

int main(void){
	FEM::CElement et(FEM::etT);
	et.setNode(0, 0.0, 0.0);
	et.setNode(1, 1.0, 0.0);
	et.setNode(2, 0.0, 1.0);
	FEM::Integrator::localMatrix(et);
    et.printMatrix();
	
	printf("\n");    
    
    FEM::CElement eq(FEM::etQ);
    eq.setNode(0, 0.0, 0.0);
    eq.setNode(1, 1.0, 0.0);
    eq.setNode(2, 1.0, 1.0);
    eq.setNode(3, 0.0, 1.0);
    FEM::Integrator::localMatrix(eq);
    eq.printMatrix();

    printf("\n");


    double nodes[8][3] = {
    	{0.0, 0.0, 0.0},
    	{1.0, 0.0, 0.0},
    	{1.0, 1.0, 0.0},
    	{0.0, 1.0, 0.0},
    	{0.0, 0.0, 1.0},
    	{1.0, 0.0, 1.0},
    	{1.0, 1.0, 1.0},
    	{0.0, 1.0, 1.0}
    };
    FEM::CElement eqp(FEM::etPrizmQ, true);
    for(int i=0; i<8; ++i){
    	eqp.setNode(i, nodes[i]);
    }
    FEM::Integrator::localMatrix(eqp);
    eqp.printMatrix();

    auto A = eqp.ae;
    const int N = eqp.N;
    double max, min; max = min = A[0][0];
    for(int i=0; i<N; ++i){
    	for(int j=0; j<N; ++j){
    		if (max < A[i][j]){
    			max = A[i][j];
    		}
    		if (min > A[i][j]){
    			min = A[i][j];
    		}
    	}
    }

    bitmap_image im(N, N);
    im.set_all_channels(255, 255, 255);
    for (int i=0; i<N; ++i){
    	for(int j=0; j<N; ++j){
    		if (A[i][j] > 0){
    			unsigned char cl = 255 - static_cast<unsigned char>(255.0*A[i][j]/max);
    			im.set_pixel(i, j, cl, cl, 255);
    		}
    		else if (A[i][j] < 0){
    			unsigned char cl =  255 - static_cast<unsigned char>(255.0*A[i][j]/min);
    			im.set_pixel(i, j, 255, cl, cl);
    		}
    	}
    }
    im.save_image("fem_test_qprizm_ae.bmp");

    Eigen::Vector3d test(0.5, -0.25, 1.0);
    std::cout << test.cwiseInverse() << std::endl;
	
	return EXIT_SUCCESS;
}
