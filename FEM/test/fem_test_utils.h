#ifndef __KRF_FEM_TEST_UTILS_H__
#define __KRF_FEM_TEST_UTILS_H__

#include <chrono>
#include <cstdio>

#include <Eigen/Sparse>
#include <image/bitmap_image.hpp>
#include <FEM/fem_mesh.h>

double timer(){
	static char is = 0;
	static std::chrono::time_point<std::chrono::high_resolution_clock> start;
	if (!is){
		is = 1;
		start = std::chrono::high_resolution_clock::now();
		return 0;
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end-start;
	start = end;
	return diff.count();
}

void spMatPortrait(Eigen::SparseMatrix<double> &A, const char *fname){
	bitmap_image im(A.rows(), A.cols());
	im.set_all_channels(255, 255, 255);

	const double *v = A.valuePtr();
	double min = v[0];
	double max = v[0];
	for(int i=0; i<A.nonZeros(); ++i){
		if (v[i] > max) max = v[i];
		if (v[i] < min) min = v[i];
	}

	for (int k=0; k<A.outerSize(); ++k){
		for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it){
			double v = it.value();
			int j = it.row();
			int i = it.col();
			
			if (v > 0){
    			unsigned char cl = 200 - static_cast<unsigned char>(200.0*v/max);
    			im.set_pixel(i, j, cl, cl, 255);
    		}
    		else if (v < 0){
    			unsigned char cl =  200 - static_cast<unsigned char>(200.0*v/min);
    			im.set_pixel(i, j, 255, cl, cl);
    		}
		}
	}

	im.save_image(fname);
}

void writeSolutionTecplot(const char *fname, const char *title, const Eigen::VectorXd &u, const FEM::CMesh &mesh){
	FILE *fp = fopen(fname, "w");
	if (!fp) return;
	int i, k;

	fprintf(fp, "TITLE = \"%s\"\n", title);
	fprintf(fp, "VARIABLES = \"X\"\n\"Y\"\n\"Z\"\n\"U\"\n");
	fprintf(fp, "ZONE T=\"Plate\"\nN=%ld, E=%ld, F=FEPOINT ET=TRIANGLE\nDT=(SINGLE SINGLE SINGLE SINGLE)\n", mesh.Nn, mesh.Ne);
	
	auto &nodes = mesh.getNodes();
	i = 0;
	for (const auto &node: nodes){
		fprintf(fp, "%le %le %le %le\n", node[0], node[1], node[2], u[i]);
		i++;
	}

	auto &elements = mesh.getElements();
	for (const auto &el: elements){
		const size_t *nodeNums = el.nodeNums;
		fprintf(fp, "%ld %ld %ld\n", nodeNums[0]+1, nodeNums[1]+1, nodeNums[2]+1);
	}
	
	fclose(fp);
}
#endif
