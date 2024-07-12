/* This is the main file for 3D FEM for computing electromagnetic scattering
 * Project begun in August 2019, with team: Siddhant, Sriram, Aggraj, Uday
 */

#include <iostream>
//#include <fftw3.h>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <sstream>
#include <cstdlib>
#include <complex>
#include <time.h>
#include <chrono>
#include <queue>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <iterator>

#include<complex_bessel.h>
#include "matplotlibcpp.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;
namespace plt = matplotlibcpp;

typedef complex<double> cmp;
typedef Eigen::Triplet<double> T;

//Constants
const double PI = 3.14159265358979323846;
const double c_0 = 299792458.0; // m/s
const double eps_0 = 8.854187817e-12; // F/m
const double mu_0 = 12.566370614e-7; // kg m/s^2 / C/s
const double Z0 = 376.730313461; //Ohm
//const cmp iota = cmp(0,1.0);

#include "class_defs.h"
#include "helper_fns.h"

int main(int argc, char* argv[])
{
	ios_base::sync_with_stdio(false); cin.tie(NULL);

	auto start = high_resolution_clock::now();

	cout<<"* ========= FEM3D v1"<<endl;

	//Correct usage
	if(argc != 2)
	{
		cout<<"** Correct way to invoke program is <exe> <config_file> **"<<endl;
		return -1;
	}

	//stock variables
	int status;
	//Create simulation input variables that will be passed to the wrapper:
	double lambda=0;
	//input mesh file
	ifstream inmesh;
	//incidence angles
	vector<double> thetas;
	//vector for tissue dielectric constants
	vector<cmp> tissues;
	// output file
	ofstream  bircs;

	double thetp, phip; //Far field point angular coordinates
	double dist_fac, a; //radius of sphere=a, farfield radius=dist_fac*a

	//Call the wrapper to read input variables
	status = input_configuration(argc, argv, &lambda, &inmesh, &thetas, &tissues, &bircs, &thetp, &phip, &a, &dist_fac);
	if(status !=0 )
	{
		cout<<"** Error in reading input configuration **"<<endl;return -1;
	}
	double theta = thetas[0];
	double phi = thetas[1];

	//Empty datastructure lists declared
	vector<node> nodelist; vector<element> elementlist; vector<edge> edgelist;
	vector<surf_element> surf_elt_list;  vector<surf_element> scat_surf_elt_list;

	/**********Create datastructures: Node, Element, Edge, Surface element, & Scatterer-surface element lists******************************/
	status = create_fem_datastructures(&inmesh,lambda,&nodelist,&elementlist,&edgelist,tissues[0], &surf_elt_list, &scat_surf_elt_list);

	/*********************Matrix Assembly: A = P - Q + R - S *****************************************************************************/
	int No=signed(edgelist.size());
	Eigen::SparseMatrix<cmp,RowMajor> Amat(No,No); //Declaring Row-Major Sparse Matrix (https://eigen.tuxfamily.org/dox/group__TutorialSparse.html)
	Amat.reserve(VectorXi::Constant(No,36));  //Reserving 36 elements per row to speed up building matrix
	Eigen::SparseVector<cmp> bvec(No);
	Amat.setZero(); bvec.setZero();

	MatrixXcd Pmat=MatrixXcd::Zero(No,No); MatrixXcd Qmat=MatrixXcd::Zero(No,No);
	MatrixXcd Rmat=MatrixXcd::Zero(No,No); MatrixXcd Smat=MatrixXcd::Zero(No,No);

	for(int e=0; e < signed(elementlist.size()); e++)
	{
		set_element_matrices(e,&nodelist,&elementlist,&edgelist);
		set_matrix_entries_PQ(e,&nodelist,&elementlist,&edgelist,lambda,&Amat,&Pmat,&Qmat);
	}

	for(int e=0; e < signed(surf_elt_list.size()); e++)
		set_matrix_entries_RS(e,&nodelist,&elementlist, &surf_elt_list, &edgelist,lambda,theta,phi,&Amat,&bvec,&Rmat,&Smat);

	/*******Uncomment to print FEM Matrix and incident field vector b**********************/
	// cout<< "\n\nFEM Matrix: "<<endl<<endl;
	// cout<< MatrixXcd(Amat)<<endl;
	// cout<<"\n\nb_vec: "<<endl;
	// cout<<VectorXcd(bvec)<<endl;
	/*********************************Matrix Inversion********************************************************************************/
	 SparseLU<SparseMatrix<cmp>, COLAMDOrdering<int>> solver;
	 solver.analyzePattern(Amat);                             // Compute the ordering permutation vector from the structural pattern of A
	 solver.factorize(Amat);                                  // Compute the numerical factorization
	 Eigen::VectorXcd xvec = solver.solve(bvec);              // Use the factors to solve the linear system
	 //cout<<xvec;

	 double relative_error = ((Amat)*xvec - (bvec)).norm() / (bvec).norm(); // norm() is L2 norm
	 cout << "\nThe relative error in matrix inversion is: " << relative_error;
	 cmp epsi = tissues[0];

	 /* Important: This loop sets area, centroid and FEM fields at centroid of exposed triangles at huygen surface */
	 /* For huygen verification without fem, fields can be changed to mie fields by editing set_surf_members fn in helper_fns.h */
	 for(int e=0; e < signed(scat_surf_elt_list.size()); e++)
		 set_surf_members(e, &scat_surf_elt_list, &nodelist, &elementlist, &edgelist, &xvec, lambda, epsi, a);

	/**Few element MMS (1,2 or 3. Currently 3. Plots true vs manufactured field)*******************************************************/
	// few_element_mms(&nodelist,&elementlist,&surf_elt_list,&edgelist,&Amat,lambda);

	/**FEM vs Mie nearfield (original problem) (Plots nearfield on huygen surface. Edit the function if you want to avoid plotting)****/
	mie_vs_fem_nearfield(epsi,&scat_surf_elt_list,lambda,a);

	/**FEM vs Mie farfield (original problem) (Plots farfield at points on a circle far off @ (dist_fac,thetp,phip))*******************/
	// mie_vs_fem_farfield(epsi,&scat_surf_elt_list,lambda,a,dist_fac,thetp,phip);

	 /*****************************************End of everything***********************************************************************/
	 auto stop = high_resolution_clock::now();
	 auto duration = duration_cast<microseconds>(stop - start);
	 cout << "\n\nProgram time: " << float(duration.count()/1000000.0) << " seconds" << endl;
	 return 0;

}

/************* Printing all data structures (include inside main() and uncomment as needed) ************************************************/
//  int count=0; //Uncomment for any printing operation
//
//  cout << "\n\nNode list:\n\n";
//  count = 0;
//  for (auto it = nodelist.cbegin(); it != nodelist.cend(); ++it)
//  {
// 		++count;
// 		cout << "Node " << count << " coordinates : " << (*it).nx << " " << (*it).ny << " " << (*it).nz << endl;
//  }
//
// cout << "\n\nEdge list with shared element info:\n\n";
// count = 0;
// for (auto it = edgelist.cbegin(); it != edgelist.cend(); ++it)
// {
// 	cout << "Edge " << count << " : " << (*it).ni << ' ' << (*it).nj<<endl; // <<" type:"<<(*it).edbtype<<endl;
// 	cout << "Its shared elements : "; print((*it).els);
// 	cout << "\n\n";
// 	++count;
// }
//
// cout << "\n\nSurface element list:\n\n";
// for (unsigned int i=0; i < surf_elt_list.size(); i++)
// {
// 	cout << "Global element no : "; print(surf_elt_list[i].glob_el_num);
// 	cout << "\nFace Nodes : "; print(surf_elt_list[i].face_nodes);
// 	cout << "\nInterior Node : "; print(surf_elt_list[i].in_node);
// 	cout<<endl<<endl;
// 	// cout << "Centroid x y z : " << surf_elt_list[i].x <<" "<< surf_elt_list[i].y <<" "<< surf_elt_list[i].z<<endl;
// 	// cout << "Centroid r theta phi : " << surf_elt_list[i].r <<" "<< surf_elt_list[i].theta <<" "<< surf_elt_list[i].phi<<endl;
// 	// cout << "Area : " << surf_elt_list[i].area<<endl<<endl;
// }
// cout<<"Size of surf list: "<<surf_elt_list.size();
//
// cout << "\n\nElement list (optionally with edge info):\n\n";
// count = 0;
// for (auto it = elementlist.cbegin(); it != elementlist.cend(); ++it)
// {
//   cout << "Element " << count << " : "; print((*it).nodes); cout <<" "<<(*it).eltype<<" "<<(*it).eleps<<endl;
//   cout << "Its edges : "; print((*it).edges);
// 	cout << endl << "Local edge to node pair : ";
// 	for (int j=0; j<6; j++)
// 		cout <<(*it).num[j][0]<<" "<<(*it).num[j][1]<<", ";
//   cout<< "\n\n";
//   ++count;
// }

/********True vs Manufac soln (multi element)********/
// int Ne = signed(elementlist.size());
// std::vector<double> x(Ne),y(Ne),z(Ne),points(Ne),mag_true(Ne); int i=0;
// for (i=0; i<Ne; i++)
// {
// 	vector<int> ein = elementlist[i].nodes;
// 	x.at(i) = (nodelist[ein[0]].nx + nodelist[ein[1]].nx + nodelist[ein[2]].nx + nodelist[ein[3]].nx)/4.0;
// 	y.at(i) = (nodelist[ein[0]].ny + nodelist[ein[1]].ny + nodelist[ein[2]].ny + nodelist[ein[3]].ny)/4.0;
// 	z.at(i) = (nodelist[ein[0]].nz + nodelist[ein[1]].nz + nodelist[ein[2]].nz + nodelist[ein[3]].nz)/4.0;
// 	points.at(i) = i+1;
// 	mag_true.at(i) = sqrt(x.at(i)*x.at(i) + y.at(i)*y.at(i) + z.at(i)*z.at(i));
// }
//
// std::vector<double> Emanfx(Ne), Emanfy(Ne), Emanfz(Ne), mag_mms(Ne);
// std::fill(Emanfx.begin(), Emanfx.end(), 0); std::fill(Emanfy.begin(), Emanfy.end(), 0); std::fill(Emanfz.begin(), Emanfz.end(), 0);
// for (i=0; i<Ne; i++)
// {
// 	manuf_at_pt(i, &elementlist, &Emanfx, &Emanfy, &Emanfz, &xvec_tilde, x.at(i), y.at(i), z.at(i),i);
// 	mag_mms.at(i) = sqrt(Emanfx.at(i)*Emanfx.at(i) + Emanfy.at(i)*Emanfy.at(i) + Emanfz.at(i)*Emanfz.at(i));
// }

/**************Single Element MMS***************************/
//
// vector<double> xforced{-0.2838,0.0423,-0.1070,0.2415,-0.1768,-0.0647};
//
// fstream fin; fin.open("XYZ.csv", ios::in);	string line;
// element ei = elementlist[0];
//
// cout<<"\n\nLocal edge to node pair:"<<endl;
// for (int i=0;i<6;i++)
// {
// 	cout<<"Edge "<<i+1<<": "<<1+ei.num[i][0]<<" "<<1+ei.num[i][1]<<endl;
// }
//
// std::vector<double> x(140),y(140),z(140),points(140),mag_true(140); char c; int i=0;
//
// while (fin && i<140)
// {
// 	getline(fin, line);
// 	stringstream s(line);
// 	s>>x.at(i)>>c>>y.at(i)>>c>>z.at(i);
// 	mag_true.at(i) = sqrt(x.at(i)*x.at(i) + y.at(i)*y.at(i) + z.at(i)*z.at(i));
// 	i++;
// }
//
// for (i=0;i<140;i++)
// 	points.at(i) = i+1;
//
// std::vector<double> Emanfx(140), Emanfy(140), Emanfz(140), mag_mms(140);
// std::vector<double> E_ls_x(140), E_ls_y(140), E_ls_z(140), mag_ls(140);
// std::vector<double> err_x(140), err_y(140), err_z(140), err_mag(140);
//
// std::fill(Emanfx.begin(), Emanfx.end(), 0);
// std::fill(Emanfy.begin(), Emanfy.end(), 0);
// std::fill(Emanfz.begin(), Emanfz.end(), 0);
//
// std::fill(E_ls_x.begin(), E_ls_x.end(), 0);
// std::fill(E_ls_y.begin(), E_ls_y.end(), 0);
// std::fill(E_ls_z.begin(), E_ls_z.end(), 0);
//
// for (int ii=0;ii<140;ii++)
// {
// 	manuf_at_pt(0, &elementlist, &Emanfx, &Emanfy, &Emanfz, &E_ls_x, &E_ls_y, &E_ls_z, &xvec_tilde, x.at(ii), y.at(ii), z.at(ii), ii, &xforced);
// 	mag_mms.at(ii) = sqrt(Emanfx.at(ii)*Emanfx.at(ii) + Emanfy.at(ii)*Emanfy.at(ii) + Emanfz.at(ii)*Emanfz.at(ii));
// 	mag_ls.at(ii) = sqrt(E_ls_x.at(ii)*E_ls_x.at(ii) + E_ls_y.at(ii)*E_ls_y.at(ii) + E_ls_z.at(ii)*E_ls_z.at(ii));
// }
