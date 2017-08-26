#ifndef _KRF_FEMMESH_H_
#define _KRF_FEMMESH_H_
#include <cstdio>
#include <stdexcept>
#include <functional>
#include <unordered_map>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Sparse>

#include "noexcept_hack.h"
#include "fem_element.h"

namespace FEM{
	using Vector3d = Eigen::Vector3d;
	namespace Mesh{
		struct TElement{
			CElement *e;
			size_t *nodeNums;

			void alloc(TElementType et, size_t eNum) noexcept{
				release();
				e = new CElement(et, eNum, true);
				nodeNums = new size_t[BaseElement::elementNodesCount(et)];
			};

			void release() noexcept{
				if (e){
					delete e;
				}
				if (nodeNums){
					delete[] nodeNums;
				}
			}

			TElement() noexcept:e(NULL), nodeNums(NULL){};

			TElement(TElement &&other) noexcept:e(other.e),nodeNums(other.nodeNums){
				other.e = NULL;
				other.nodeNums = NULL;
			};

			TElement& operator=(TElement &&other) noexcept{
				e = other.e; other.e = NULL;
				nodeNums = other.nodeNums; other.nodeNums = NULL;
				return *this;
			};

			TElement(const TElement &) = delete;
			TElement& operator=(const TElement&) = delete;


			~TElement(){
				release();
			}
		};

		struct TGrane{
			const CElement& meshElement;
			CElement element;
			Vector3d normal;
			size_t *nodeNums;

			void release() noexcept{
				if (nodeNums){
					delete[] nodeNums;
					nodeNums = NULL;
				}
			}

			TGrane(const TGrane &) = delete;
			void operator=(const TGrane&) = delete;

			TGrane(TGrane &&other) noexcept: 
				meshElement(other.meshElement), 
			    element(std::move(other.element)), normal(std::move(other.normal)),
			    nodeNums(other.nodeNums)
			{
				other.nodeNums = NULL;
			}
			TGrane(const TElement &me, CElement &e, Vector3d &n, const size_t* nNums) noexcept:
				meshElement(*me.e), 
			    element(std::move(e)), normal(std::move(n))
			{
				nodeNums = new size_t[e.N];
				for(int i=0; i<e.N; ++i){
					nodeNums[i] = me.nodeNums[nNums[i]];
				}
			}

			~TGrane(){
				release();
			}
		};
	};

	class CMesh{
	private:
		std::vector<Vector3d,Eigen::aligned_allocator<Vector3d> > nodes;
		std::vector<Mesh::TElement> elements;
		std::vector<std::vector<size_t> > nodeElements;
		std::vector<Mesh::TGrane> boundaryGranes;
		std::vector<bool> isBoundaryNode;

		int elementNodesCount;

		friend class Assembler;
		friend class SimpleMeshReader;
		
		CMesh(const CMesh&) = delete;
    	void operator=(const CMesh&) = delete;

    	void _setElementsNodes() noexcept{
			for(size_t i=0; i<Ne; ++i){
				Mesh::TElement &e = elements[i];
				for(int j=0; j<e.e->N; ++j){
					e.e->setNode(j, nodes[e.nodeNums[j]].data());
				}
			}
		};

		void _sortNodeElements() noexcept{
			for(auto &ne: nodeElements){
				std::sort(ne.begin(), ne.end());
			};
		};

		void _boundaryGranes() noexcept{
			size_t nodeNums[5];

			struct Grane{size_t n1, n2, eNum; int gNum;};
			using T = std::unordered_multimap<size_t, Grane>;

			T granes;

			for(Mesh::TElement &me: elements){
				const TElementType et = me.e->type;
				const int ng = BaseElement::elementGranesCount(et);

				for (int i=0; i<ng; ++i){
					int nn = BaseElement::elementGraneNodeNumbers(et, i, (size_t*)nodeNums);
					for(int j=0; j<nn; ++j){
						nodeNums[j] = me.nodeNums[nodeNums[j]];
					}

					std::sort(&nodeNums[0], &nodeNums[nn]);

					for(int j=nn; j<3; ++j){
						nodeNums[j] = nodeNums[nn-1];
					}

					bool add = true;
					auto range = granes.equal_range(nodeNums[0]);
					for(auto it=range.first; it != range.second; ++it){
						Grane &g = it->second;
						if (g.n1 == nodeNums[1] && g.n2 == nodeNums[2]){
							granes.erase(it);
							add = false;
							break;
						}
					}

					if (add){
						Grane g = {nodeNums[1], nodeNums[2], me.e->num, i};
						granes.emplace(nodeNums[0], g);
					}
				}
			}

			for (auto& x: granes){
				Vector3d normal;
				const Mesh::TElement& me = elements[x.second.eNum];

				CElement eg = me.e->graneElement(x.second.gNum, normal, nodeNums);
				boundaryGranes.emplace_back(me, eg, normal, nodeNums);
			}

			for (auto &g: boundaryGranes){
				for(int i=0; i<g.element.N; ++i){
					isBoundaryNode[g.nodeNums[i]] = true;
				}
			}
		};
	public:
		const TDim dim;
		const size_t Nn, Ne;
		CMesh(TDim pdim, size_t pNn, size_t pNe) noexcept: dim(pdim), Nn(pNn), Ne(pNe){
			         nodes.reserve(Nn);          nodes.resize(Nn);
			      elements.reserve(Ne);       elements.resize(Ne);
			  nodeElements.reserve(Nn);   nodeElements.resize(Nn);
			isBoundaryNode.reserve(Nn); isBoundaryNode.resize(Nn);

			for(size_t i=0; i<Nn; ++i){
				isBoundaryNode.emplace_back(false);
			}

			elementNodesCount = 0;
		};

		CMesh(CMesh &&other) noexcept:
			nodes(std::move(other.nodes)), 
		    elements(std::move(other.elements)),
		    nodeElements(std::move(other.nodeElements)),
		    boundaryGranes(std::move(other.boundaryGranes)),
		    isBoundaryNode(std::move(other.isBoundaryNode)),
		    elementNodesCount(other.elementNodesCount),
		    dim(other.dim),
		    Nn(other.Nn),
		    Ne(other.Ne)
		{
			//
		};

		const std::vector<Mesh::TGrane>& getBoundaryGranes() const noexcept{
			return boundaryGranes;
		}

		const std::vector<Vector3d,Eigen::aligned_allocator<Vector3d> >& getNodes() const noexcept{
			return nodes;
		}

		const std::vector<Mesh::TElement>& getElements() const noexcept{
			return elements;
		}

		void setNode(size_t i, double x, double y, double z) noexcept{
			assert(i < Nn);

			nodes[i] << x,y,z;
		};

		void setElement(size_t i, TElementType et, const size_t *nodeNums){
			assert(i < Ne);

			const int n = BaseElement::elementNodesCount(et);
			for(int j=0; j<n; ++j){
				if (nodeNums[j] >= Nn){
					throw std::invalid_argument("Bad node number in list of element nodes");
				}
			}
			
			Mesh::TElement &e = elements[i];
			e.alloc(et, i);
			memcpy(e.nodeNums, nodeNums, sizeof(size_t)*n);

			for(int j=0; j<n; ++j){
				nodeElements[nodeNums[j]].push_back(i);
			}
		};

		void postprocessing() noexcept{
			_setElementsNodes();
			_sortNodeElements();
			_boundaryGranes();
		};
	};

	class SimpleMeshReader{
    private:
        char line[1025];
        static constexpr const char *const empty = " \t\r\n\0";
        char* stripSpace(char *str){
        	char *start = str;
        	
			while(*start && strchr(empty, *start)){
				start++;
			}

			return (*start) ? start : NULL;
        };
        char* nextSpace(char *str){
        	char *start = str;

        	while(*start && !strchr(empty, *start)){
        		start++;
        	}

        	return (*start) ? start : NULL;
        };
        char* nextToken(char *str){
        	char *space = nextSpace(str);
        	if (!space) return NULL;

        	return stripSpace(space);
        };
        int nextSection(FILE *fp){
        	char *p = NULL;
            while((p = nextLine(fp))){
                if (*p == '$'){
                	break;
                }
            }

            if (!p){
            	return -2;
            }

            char *sectionNum = stripSpace(p+1);
        	if (!sectionNum){
        		return -1;
        	}

        	int sn;
        	if (1 != sscanf(sectionNum, "%d", &sn)){
        		return -1;
        	}

            return sn;
        };

        char* nextLine(FILE *fp){
        	while(NULL != fgets(line, 1024, fp)){
        		char *p = stripSpace(line);
        		if (p && *p && (p[0] != '/' || p[1] != '/')){
        			return p;
        		}
        	}
        	return NULL;
        };

        void throwClose(FILE *fp, const char *what){
        	if (fp){
        		fclose(fp);
        	}
        	throw std::invalid_argument(what);
        };
	public:
		CMesh read(const char *fname, TDim dim){
			memset(line, 0, 1025);
			FILE *fp = fopen(fname, "rb");
			if (!fp){
                throw std::invalid_argument("Unable open file");
			}
            
            size_t Ne, Nn; Ne = Nn = 0;
            while(!Ne || !Nn){
            	int sn = nextSection(fp);
            	if (sn == -2){
            		break;
            	}
            	if (sn < 0){
            		continue;
            	}

            	if (sn == 1 || sn == 2){
            		char *header = nextLine(fp);
            		if (!header || *header != '#'){
            			throwClose(fp, "Undefined section header");
            		}

            		header = stripSpace(header+1);
            		if (!header){
            			throwClose(fp, "Empty section header");
            		}

            		size_t *ptr = (sn == 1) ? &Nn : &Ne;
            		if (1 != sscanf(header, "%ld", ptr)){
            			throwClose(fp, "Not number in section header");
            		}
            	}
            }

            if (!Ne){
            	throwClose(fp, "Undefined nodes number");
            }
            if (!Nn){
            	throwClose(fp, "Undefined nodes number");
            }

            CMesh mesh(dim, Nn, Ne);
            bool isNodes, isElements; isNodes = isElements = false;
            fseek(fp, 0, SEEK_SET);
            while(!isNodes || !isElements){
            	int sn = nextSection(fp);
            	if (sn == -2){
            		break;
            	}
            	if (sn < 0){
            		continue;
            	}

            	if (sn == 1){ //read nodes
            		char *str = nextLine(fp); // strip header
            		for(size_t i=0; i<Nn; ++i){
            			str = nextLine(fp);
            			if (!str){
            				throwClose(fp, "Error while node reading");
            			}

            			double x, y, z; y = z = 0;
            			size_t n = sscanf(str, "%le %le %le", &x, &y, &z);
            			if (n < dim){
            				throwClose(fp, "Error while node reading: wrong coordinates number");
            			}

            			mesh.setNode(i, x, y, z);
            		}

            		isNodes = true;
            	}// read nodes
            	else if (sn == 2){ // read elements
            		char *str = nextLine(fp); // strip header
            		for(size_t i=0; i<Ne; ++i){
            			str = nextLine(fp);
            			if (!str){
            				throwClose(fp, "Error while element reading");
            			}

            			TElementType et;
            			if (!sscanf(str, "%d", &et)){
            				throwClose(fp, "Undefined element type");
            			}
            			if (et < etL || et >= etCount){
            				throwClose(fp, "Unknown element type");
            			}

            			if (dim != BaseElement::elementDim(et)){
            				throwClose(fp, "Element of illegal dimension");
            			}

            			int n = BaseElement::elementNodesCount(et);
            			size_t nodeNums[64];
            			for(int j=0; j<n; ++j){
            				str = nextToken(str);
            				if (!str){
            					throwClose(fp, "Not enougth nodes for element");
            				}

            				if (!sscanf(str, "%ld", nodeNums+j)){
            					throwClose(fp, "Can't read integer node in element");
            				}

            				if (nodeNums[j] < 1 || nodeNums[j] > Nn){
            					throwClose(fp, "Illegal node number");
            				}

            				nodeNums[j]--;
            			}
            			mesh.elementNodesCount += n;
            			mesh.setElement(i, et, nodeNums);
            		}
            		isElements = true;
            	}// read
            }

            mesh.postprocessing();

            fclose(fp);

            return mesh;
        };
	};

	class SimpleMeshGenerator{
		public:
			CMesh onRectangle(const Eigen::Vector2d &x0, const Eigen::Vector2d &L, const Eigen::Vector2i &N){
				auto x = Eigen::ArrayXd::LinSpaced(N[0], x0[0], L[0] + x0[0]);
				auto y = Eigen::ArrayXd::LinSpaced(N[1], x0[1], L[1] + x0[1]);

				const size_t Ne = (N[0]-1)*(N[1]-1)*2;
				const size_t Nn = N[0]*N[1];

				CMesh ret(dim2D, Nn, Ne);
				int k = 0;
				for(int i=0; i<N[1]; ++i){
					for(int j=0; j<N[0]; ++j){
						ret.setNode(k, x[j], y[i], 0.0);
						k++;
					}
				}

				k = 0;
				for(int i=0; i<N[1]-1; ++i){
					for(int j=0; j<N[0]-1; ++j){
						const size_t lb = i*N[0] + j;

						const int dir = rand()%2;
						if (dir){
							{
								const size_t nodeNums[] = {lb, lb+1, lb+N[0]};
								ret.setElement(k, etT, nodeNums); k++;
							}
							{
								const size_t nodeNums[] = {lb+1, lb+N[0]+1, lb+N[0]};
								ret.setElement(k, etT, nodeNums); k++;
							}
						}
						else{
							{
								const size_t nodeNums[] = {lb, lb+N[0]+1, lb+N[0]};
								ret.setElement(k, etT, nodeNums); k++;
							}
							{
								const size_t nodeNums[] = {lb, lb+1, lb+N[0]+1};
								ret.setElement(k, etT, nodeNums); k++;
							}

						}
					}
				}

				ret.postprocessing();
				return ret;
			}
	};
}; //namespace FEM
#endif
