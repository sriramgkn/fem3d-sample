double dtor(double d)
{
	return PI*d/180.0;
}

double rtod(double r)
{
	return 180.0*r/PI;
}

string stripe_spaces_after_equals(string line){
	if( line.find("=") == string::npos ){
		cout<<"** Found no = sign in string, possible error in config file **"<<endl;
	}
	// extract string after = sign
	string tmp = string(line,1+line.find("="),string::npos);
	// remove white spaces
	tmp.erase(remove(tmp.begin(),tmp.end(),' '),tmp.end());
	return tmp;
}

//This function is a wrapper that reads the config file for a substrate
//heterogeneous in y, heterogeneous in x
int input_configuration(int argc, char* argv[], double * lambda, ifstream * in, vector<double> * thetas,
vector<cmp> * tissues, ofstream * outp, double *thetp, double *phip, double *a, double *dist_fac)
{
	//Read the config file
    ifstream config(argv[1]);
	if( config.is_open() ){
		cout<<"Config: '"<<argv[1]<<"'";
	}
	else{
		cout<<"** Error opening config file **"<<endl;
		return -1;
	}

	bool onetime = true, real_specified = false, imag_specified = false;
	double eps_real, eps_imag;
	string line,tmp;
	stringstream buffer, sline;
	buffer<<config.rdbuf();

	while(buffer)
	{
		while(getline(buffer,line,'\n'))
		{
						if(line.size()==0){
							continue; //if there is an empty line
						}
						if( string(line,0,1) == "#"){
							if( onetime ){
								cout<<", description: '"<<string(line,1,string::npos)<<"'"<<endl;
								onetime = false;
							}
							continue; // all comments begin with #
						}
						if( line.find("lambda") != string::npos ){
							tmp = stripe_spaces_after_equals(line);
							(*lambda) = atof(tmp.c_str());
							cout<<"Lambda="<<100.0*(*lambda)<<"cm, ";
						}
						if( line.find("meshfile") != string::npos )
						{
							tmp = stripe_spaces_after_equals(line);
							cout<<"Mesh: '"<<tmp<<"'";
							(*in).open(tmp.c_str());
						}

						if( line.find("theta") != string::npos)
						{
							tmp = stripe_spaces_after_equals(line);
							(*thetas).push_back(dtor(atof(tmp.c_str())));
							cout<<", Inc angle(s)=";
							cout<< atof(tmp.c_str())<<", ";
						}
						if( line.find("phi") != string::npos)
						{
							tmp = stripe_spaces_after_equals(line);
							(*thetas).push_back(dtor(atof(tmp.c_str())));
							cout<< atof(tmp.c_str());
						}
						if( line.find("thetp") != string::npos ){
							tmp = stripe_spaces_after_equals(line);
							(*thetp) = atof(tmp.c_str());
						}
						if( line.find("php") != string::npos ){
							tmp = stripe_spaces_after_equals(line);
							(*phip) = atof(tmp.c_str());
						}
						if( line.find("radius_a") != string::npos ){
							tmp = stripe_spaces_after_equals(line);
							(*a) = atof(tmp.c_str());
						}
						if( line.find("dist_fac") != string::npos ){
							tmp = stripe_spaces_after_equals(line);
							(*dist_fac) = atof(tmp.c_str());
						}
						//Read in tissue dielectric constants
						if( line.find("real") != string::npos ){
								tmp = stripe_spaces_after_equals(line);
								//tis_eps.real() = double(atof(tmp.c_str()));
							  eps_real = double(atof(tmp.c_str()));
								cout<<", tissue"<<"=("<<eps_real; //tis_eps.real();
								real_specified = true;
							}
							else if( line.find("imag") != string::npos ){
								tmp = stripe_spaces_after_equals(line);
								//tis_eps.imag() = double(atof(tmp.c_str()));
								eps_imag = double(atof(tmp.c_str()));
								cout<<","<<eps_imag<<")"; //tis_eps.imag()<<")";
								imag_specified = true;
						}
						cmp tis_eps(eps_real, eps_imag);
						if(real_specified && imag_specified){
							(*tissues).push_back(tis_eps);
							real_specified = false; imag_specified = false;
						}
						//Add any more config file read-ins here
						//
						if( line.find("outputfile") != string::npos ){
							tmp = stripe_spaces_after_equals(line);
							cout<<",  output: '"<<tmp<<"'";
							(*outp).open(tmp.c_str());
						}
		}
		line.clear();
	}
	//Check if a required item is not specified.
	if( abs((*lambda)) <1e-3 ){
		cout<<endl<<"** Error wavelength not specified **"<<endl;
	}
	if( !(*in).is_open() ){
		cout<<endl<<"** Error in opening mesh file **"<<endl;
		return -1;
	}
	if( (*tissues).size() == 0 ){
		cout<<endl<<"** Tissues not specified **"<<endl;
		return -1;
	}
	return 0;
}

void create_edgelist(vector<node> * nodelist, vector<element> * elementlist, vector<edge> * edgelist)
{
	int n_idx[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
	unsigned int n1, n2, i, j;
	long large=1e10; long count;
	unsigned int N = (*elementlist).size();
	vector<long> indicator;
	indicator.reserve(6*N); // "indicator" vector size = 6N
	// Iterate over 6N edges to initialize "indicator"
	for(i=0; i < 6*N; i++)
	{
		n1 = Min((*elementlist)[i/6].nodes[n_idx[i%6][0]],(*elementlist)[i/6].nodes[n_idx[i%6][1]]);
		n2 = Max((*elementlist)[i/6].nodes[n_idx[i%6][0]],(*elementlist)[i/6].nodes[n_idx[i%6][1]]);
		indicator.push_back(n1 + n2*large);
	}
	/* Remove duplicates from "indicator" by assigning it
	 to an unordered set "red". O(N) operation! */
	unordered_set<long> red(indicator.begin(),indicator.end());
	unordered_map<long,long> red_map; // Define an unordered map "red_map"
	red_map.reserve(red.size()); // "red_map" has the same size as "red"

  count=0;
	for (auto it = red.cbegin(); it != red.cend(); ++it) // Iterate over "red" to initialize "edgelist" and "red_map"
	{
		edge ed((*it)/large,(*it)%large);
		(*edgelist).push_back(ed);	 //"edgelist" initialized with unique edge objects
		red_map[*it] = count; //"red_map" initialized with first value same as that of "red", while second value is an index (global edge number)
		++count;
	}
	red.clear();
	vector< vector<int> > elts(red_map.size()); // Define a vector of vector of ints named "elts"
	// Iterate over "indicator" (6N) to store element info of each edge in "elts"
	count=0;
	for (auto it = indicator.cbegin(); it != indicator.cend(); ++it)
	{
		elts[red_map.find(*it)->second].push_back(count/6);
		++count;
	}
	indicator.clear(); indicator.shrink_to_fit(); red_map.clear();
	// Iterate over "edgelist" for element-edge correspondence
	for (i=0; i < (*edgelist).size(); i++)
	{
		for(j=0; j < elts[i].size(); j++)
		(*elementlist)[elts[i][j]].edges.push_back(i); // Edge info of each element stored in "edges" : a vector in element class

		(*edgelist)[i].els.reserve(elts[i].size());
		(*edgelist)[i].els = elts[i]; // Element info of each edge stored in "els" : a vector in edge class
		elts[i].clear(); elts[i].shrink_to_fit();
	}
	int Ns=(*edgelist).size(); cout<<", Edges("<<Ns<<")."; // Print the number of unique edges

}

int create_fem_datastructures(ifstream * in, double lambda, vector<node> * nodelist, vector<element> * elementlist,\
	 		vector<edge> * edgelist, cmp diel_const, vector<surf_element> * surf_elt_list, vector<surf_element> * scat_surf_elt_list)
{
	//cout<< "Dielectric constant: "<<diel_const<<endl;
	int No=0,Ne=0;

	string line,tmp;
	stringstream buffer, sline;
	bool read_no = false, read_el = false;
	vector<int> temp;

	//First pass through the file to figure out size of lists and allocate memory.
	buffer<<(*in).rdbuf();
	while(buffer)
	{
				while(getline(buffer,line,'\n'))
				{
					if( string(line,0,5) == "*NODE")
					{read_no=true;continue;}
					else if( string(line,0,37) == "******* E L E M E N T S *************")
					{read_no=false;read_el=true;continue;}
					else if( string(line,0,19) == "*ELEMENT, type=C3D4") // "C3D4" means tetrahedral element
					{read_el=false;break;}

					if(read_no){No++;}
					else if(read_el){Ne++;}
				}
	}
	if(No == -1)
	{
		cout<<"**Error in file/no nodes read**"<<endl; return -1;
	}
	//rewind the buffer after clearing error flags
	buffer.clear();	buffer.seekg(0);
	//Allocate memory
	(*nodelist).reserve(No);
	(*elementlist).reserve(Ne);
 	while(buffer)
	{
				getline(buffer, line);
				if (line.substr(0, 5)=="*NODE") //Finding Nodes
				{
			    	 while(buffer)
			       {

						   		 getline(buffer, line);
				           if(line=="******* E L E M E N T S *************")
				           	break;

	                 istringstream iss(line);
	                 double num,x,y,z;
	                 char c;
	                 iss >> num  >> c >> x >> c >> y >> c >> z;//reading each line

	                 node nd(2*lambda*x,2*lambda*y,2*lambda*z);

	                 nd.ndtype='i'; // By default, every node is assigned to be interior (why?) anyway, this var is not needed

	                 (*nodelist).push_back(nd);//pushing into nodelist
	    			 	}
					}
					// Note: Volume 95 is scatterer.
					// Volume 52 = volume between scatterer surface and huygen surface
					// Volume 94 = volume between huygen surface and surface of computational domain
					else if (line.substr(0, 19)=="*ELEMENT, type=C3D4") //Finding Elements
					{
								 char medium = 'a';	 // Initialize each element with medium air

								 cmp epsr=1, mur=1;  // and relative permittivity as well as permaebility to be 1
								 while(buffer)  //Read until the end of file
						   	 {
							 				getline(buffer, line);
											if (line.substr(0, 35)=="*ELEMENT, type=C3D4, ELSET=Volume94")
											{
												getline(buffer, line);
											}
			           			if (line.substr(0, 35)=="*ELEMENT, type=C3D4, ELSET=Volume95")   //elements inside the scatterer
			                {
						   				 		medium = 's'; epsr = diel_const; getline(buffer, line);
											}
			 								istringstream iss(line);

										  int num,n1,n2,n3,n4;
										  char c;
										  iss >> num  >> c >> n1 >> c >> n2 >> c >> n3 >> c >> n4;
										  element el(0,0,0,0);
										  el.nodes[0]=n1-1; el.nodes[1]=n2-1; el.nodes[2]=n3-1; el.nodes[3]=n4-1;   //subtract 1 from node numbers

										  el.eltype = medium; el.eleps=epsr; el.elmu=mur;

		                  (*elementlist).push_back(el);
							 		}
 									(*elementlist).pop_back();  // Remove the blank line read from .inp file
				    }
			}

			No=(*nodelist).size(); cout<<"\nCreated: Nodes("<<No<<")";cout.flush();
			Ne=(*elementlist).size(); cout<<", elements("<<Ne<<")";cout.flush();

			create_edgelist(nodelist,elementlist,edgelist);  // Calling the function for creating edgelists

      buffer.clear(); buffer.seekg(0, ios::beg); // Clear the buffer and read from the beginning
      while(buffer)
      {
						 getline(buffer, line);
						 //Surfaces 20 to 25 are computational domain surfaces (outer cube)
						 //Surfaces 68 to 75 are scatterer surfaces (sphere)

						 //In the following if statement, numbers "14" and "20" are specific to the current types
						 //of mesh files named of the form "mesh_lc_alpha.inp" and "three_elt_mesh.inp"
						 //(the logic needs to be re-thought for other types of mesh files)
						 if (line=="*ELEMENT, type=CPS3, ELSET=Surface14" || line=="*ELEMENT, type=CPS3, ELSET=Surface20")
						 {
									 while(buffer)
									 {
												 getline(buffer,line);
												 if(line.substr(0,33)=="*ELEMENT, type=C3D4, ELSET=Volume" || line=="*ELEMENT, type=CPS3, ELSET=Surface68")
													 break;
												 if(line.substr(0,34)=="*ELEMENT, type=CPS3, ELSET=Surface")
													 getline(buffer,line);

												 istringstream iss(line);
												 int num,x[3];	 char c;
												 iss >> num  >> c >> x[0] >> c >> x[1] >> c >> x[2];
												 int n0, n1, n2, n3;
												 surf_element el;
												 el.face_nodes.push_back(x[0]-1); el.face_nodes.push_back(x[1]-1); el.face_nodes.push_back(x[2]-1);
												 for (int i=0; i<Ne; i++)
												 {
																n0 = (*elementlist)[i].nodes[0]; n1 = (*elementlist)[i].nodes[1];
																n2 = (*elementlist)[i].nodes[2]; n3 = (*elementlist)[i].nodes[3];

																if ( (x[0]-1==n0 || x[0]-1==n1 || x[0]-1==n2 || x[0]-1==n3) &&
																	 (x[1]-1==n0 || x[1]-1==n1 || x[1]-1==n2 || x[1]-1==n3) &&
																	 (x[2]-1==n0 || x[2]-1==n1 || x[2]-1==n2 || x[2]-1==n3) )
																{
																	el.glob_el_num.push_back(i);
																	el.in_node.push_back(n0 + n1 + n2 + n3 - x[0] - x[1] - x[2] + 3);
																}
													}
													(*surf_elt_list).push_back(el); //Initializing surface element list @ comp domain surface
											}
						 }
						 //Surfaces 77,79,172,174,176,777 are the huygen cube surfaces (inner cube)
						 if (line.substr(0,36)=="*ELEMENT, type=CPS3, ELSET=Surface77")	//77 is the first huygen cube surface in inp file, 777 last
						 {
										 while(buffer)
										 {
										 		 getline(buffer,line);
												 if(line=="*ELEMENT, type=C3D4, ELSET=Volume52")
													 break;
												 if(line.substr(0,34)=="*ELEMENT, type=CPS3, ELSET=Surface")
													 getline(buffer,line);

												 istringstream iss(line);
												 int num,x[3];	 char c;
												 iss >> num  >> c >> x[0] >> c >> x[1] >> c >> x[2];
												 int n0, n1, n2, n3;
												 surf_element el;
												 el.face_nodes.push_back(x[0]-1); el.face_nodes.push_back(x[1]-1); el.face_nodes.push_back(x[2]-1);
												 for (int i=0; i<Ne; i++)
												 {
																n0 = (*elementlist)[i].nodes[0]; n1 = (*elementlist)[i].nodes[1];
																n2 = (*elementlist)[i].nodes[2]; n3 = (*elementlist)[i].nodes[3];

																if ( (x[0]-1==n0 || x[0]-1==n1 || x[0]-1==n2 || x[0]-1==n3) &&
																	 (x[1]-1==n0 || x[1]-1==n1 || x[1]-1==n2 || x[1]-1==n3) &&
																	 (x[2]-1==n0 || x[2]-1==n1 || x[2]-1==n2 || x[2]-1==n3) )
																{
																	el.glob_el_num.push_back(i);
																	el.in_node.push_back(n0 + n1 + n2 + n3 - x[0] - x[1] - x[2] + 3);
																}
												  }
													(*scat_surf_elt_list).push_back(el); //convention: scat_surf_elt_list = list of elements *on huygen surface*
											}
								}
				}

				buffer.clear();buffer.str(""); // Close the buffer; File reading is done.
				(*in).close();
				if(No == 1)
				{
					cout<<"** Error in file/no nodes read **"<<endl; return -1;
				}
				return 0;
}

void print(std::vector<int> const &input)
{
	for (int i = 0; i < signed(input.size()); i++)
		std::cout << input.at(i) << ' ';
}

//This function sets the element matrices for all the elements
void set_element_matrices(int e, vector<node> * nodelist, vector<element> * elementlist, vector<edge> * edgelist)
{

		element ei = (*elementlist)[e];  // e^th element
		int vec[4][3]={{2,3,4},{3,1,4},{1,2,4},{2,1,3}};  // array to store the local node indices for finding a,b,c,d determinant
		int i, j;
		Matrix3d A,B,C,D;
		node nd;
		//ei.a.reserve(4); ei.b.reserve(4); ei.c.reserve(4); ei.d.reserve(4);//Allocating memory for a,b,c,d
		for(i=0; i<4; i++)
		{
		    for (j=0; j<3; ++j)
		    {
					nd = (*nodelist)[ei.nodes[vec[i][j]-1]];  // declaring nodes belonging to the element
			    A(0,j) = nd.nx; A(1,j) = nd.ny; A(2,j) = nd.nz;  // Defining coefficients for matrix A,B,C,D
			    B(0,j) = 1;     B(1,j) = nd.ny; B(2,j) = nd.nz;
			    C(0,j) = 1;     C(1,j) = nd.nx; C(2,j) = nd.nz;
			    D(0,j) = 1;     D(1,j) = nd.nx; D(2,j) = nd.ny;
			  }
		    //coefficients for scalar basis : L(x,y,z) = a + bx + cy + dz
		    (*elementlist)[e].a.push_back(A.determinant()); (*elementlist)[e].b.push_back(-B.determinant());
		    (*elementlist)[e].c.push_back(C.determinant()); (*elementlist)[e].d.push_back(-D.determinant());  //pushing into each element
		}

		Matrix4d V_mat;
	   // Finding volume of tetrahedron
	   for(i=0; i<4; i++)
	   {
		    V_mat(0,i) = 1;
				V_mat(1,i) = (*nodelist)[ei.nodes[i]].nx;
				V_mat(2,i) = (*nodelist)[ei.nodes[i]].ny;
				V_mat(3,i) = (*nodelist)[ei.nodes[i]].nz;
			//cout <<'\t'<< V_mat(1,i) <<' '<< V_mat(2,i) << ' ' << V_mat(3,i) << '\t' << endl;
	   }

	   (*elementlist)[e].V = V_mat.determinant()/6.0;  // Volume of tetrahedral element

	   double l[6]; node ni,nj; int num[6][2];

	   for(i=0;i<6;i++)  // Loop over 6 edges of ei
	   {
					for(j=0;j<4;j++)  // Loop over 4 nodes of ei
					{
						for(int k=0;k<4;k++)
						{
									if ((*edgelist)[ei.edges[i]].ni==ei.nodes[j] && (*edgelist)[ei.edges[i]].nj==ei.nodes[k])
									{
											if ( ei.nodes[j] > ei.nodes[k] )
												{num[i][0] = j; num[i][1] = k;}
											else
												{num[i][0] = k; num[i][1] = j;}
											goto edgedef;
									}
						}
					}
					edgedef: ni=(*nodelist)[ei.nodes[num[i][0]]], nj=(*nodelist)[ei.nodes[num[i][1]]];  // Edge from global node ni to global node nj

				  l[i]=sqrt(sq(ni.nx-nj.nx)+sq(ni.ny-nj.ny)+sq(ni.nz-nj.nz)); // Edge Length
				  (*elementlist)[e].l.push_back(l[i]);
					(*elementlist)[e].num[i][0] = num[i][0]; (*elementlist)[e].num[i][1] = num[i][1];
	    }

}

void calc_normal(surf_element *es, vector<node> *nodelist, vector<edge> *edgelist, vector<node> *surf_nd, node *int_nd, vector<double> *nlvec, double *area)
{

		  for (int i=0; i<3; i++)
			(*surf_nd).push_back( (*nodelist)[(*es).face_nodes[i]] ); // Surface nodes of the exposed elt

		  (*int_nd) = (*nodelist)[(*es).in_node[0]]; // interior node of the exposed elt


	     // Finding normal to the surface
	     double nlx, nly, nlz, e1_x, e1_y, e1_z, e2_x, e2_y, e2_z, e3_x, e3_y, e3_z;// e4_x, e4_y, e4_z; //nlx,nly,nlz define the unit normal

       // Surface edge vectors (e1_x, e1_y, e1_z) and (e2_x, e2_y, e2_z)
	     e1_x = (*surf_nd)[1].nx - (*surf_nd)[0].nx; e1_y = (*surf_nd)[1].ny - (*surf_nd)[0].ny; e1_z = (*surf_nd)[1].nz - (*surf_nd)[0].nz;
	     e2_x = (*surf_nd)[2].nx - (*surf_nd)[0].nx; e2_y = (*surf_nd)[2].ny - (*surf_nd)[0].ny; e2_z = (*surf_nd)[2].nz - (*surf_nd)[0].nz;

	     // Internal edge vector (e3_x, e3_y, e3_z)
	     e3_x = (*int_nd).nx - (*surf_nd)[0].nx; e3_y = (*int_nd).ny - (*surf_nd)[0].ny; e3_z = (*int_nd).nz - (*surf_nd)[0].nz;

	     //e3_x=1; e3_y=3; e3_z=4;

	     // Unnormalized normal vector [cross product of surface edge vectors] [e1 x e2]
	     nlx = e1_y*e2_z - e1_z*e2_y; nly = e1_z*e2_x - e1_x*e2_z; nlz = e1_x*e2_y - e1_y*e2_x;

	     //e4_x = (*surf_nd)[2].nx - (*surf_nd)[1].nx; e4_y = (*surf_nd)[2].ny - (*surf_nd)[1].ny; e4_z = (*surf_nd)[2].nz - (*surf_nd)[1].nz;


	     double nl_mag = sqrt(nlx*nlx + nly*nly + nlz*nlz);
	     nlx = nlx/nl_mag; nly = nly/nl_mag; nlz = nlz/nl_mag;  // Unit normal vector

       (*area)=0.5*nl_mag; // Area of exposed triangle

	     if (nlx*e3_x + nly*e3_y + nlz*e3_z>0) // If normal points inward...
	     {nlx = -nlx; nly = -nly; nlz = -nlz; }   // Reverse the direction of normal

	     (*nlvec).push_back(nlx); (*nlvec).push_back(nly); (*nlvec).push_back(nlz);

}

void set_matrix_entries_PQ(int e, vector<node> * nodelist, vector<element> * elementlist, vector<edge> * edgelist,\
						double lambda, SparseMatrix<cmp,RowMajor> * Amat, MatrixXcd *Pmat, MatrixXcd *Qmat)  // MatrixXcd * Amat
{
	   element ei = (*elementlist)[e];  // e^th element

	   int i, j, i1, i2, j1, j2;

	   MatrixXcd P(6,6), Q(6,6); //Complex matrices of size 6x6

	   //*****************************************************************************//
      // Calculating P
      for(i=0; i<6; i++)
      {
				   for(int j=0; j<6; j++)
				   {

					   	i1=ei.num[i][0]; i2=ei.num[i][1]; j1=ei.num[j][0]; j2=ei.num[j][1]; //Indices for nodes i1, i2, j1, j2;  Edge i: Node i1 -> Node i2

							P(i,j) = (ei.l[i]*ei.l[j]/(324.0*pow(ei.V,3.0)))* ( (ei.c[i1]*ei.d[i2] - ei.c[i2]*ei.d[i1]) * (ei.c[j1]*ei.d[j2] - ei.c[j2]*ei.d[j1])\
					                                   + (ei.d[i1]*ei.b[i2] - ei.d[i2]*ei.b[i1]) * (ei.d[j1]*ei.b[j2] - ei.d[j2]*ei.b[j1])\
					                                   + (ei.b[i1]*ei.c[i2] - ei.b[i2]*ei.c[i1]) * (ei.b[j1]*ei.c[j2] - ei.b[j2]*ei.c[j1]) );

							(*Amat).coeffRef(ei.edges[i],ei.edges[j]) += P(i,j);

			        (*Pmat)(ei.edges[i],ei.edges[j])+=P(i,j);

			    	}
	   	 }


     //*****************************************************************************//
	   // Calculating Q
	   Matrix3d Jac;	  //Jacobian Matrix
	   for(i=0; i<3; i++)
	   {
			Jac(0,i) = (*nodelist)[ei.nodes[i+1]].nx - (*nodelist)[ei.nodes[0]].nx;
			Jac(1,i) = (*nodelist)[ei.nodes[i+1]].ny - (*nodelist)[ei.nodes[0]].ny;
			Jac(2,i) = (*nodelist)[ei.nodes[i+1]].nz - (*nodelist)[ei.nodes[0]].nz;
	   }
	   double J = Jac.determinant();  // Jacobian

     node n1 = (*nodelist)[ei.nodes[0]], n2 = (*nodelist)[ei.nodes[1]], n3 = (*nodelist)[ei.nodes[2]], n4 = (*nodelist)[ei.nodes[3]];
		 double n[4],f[4],g[4],h[4];

     for(i=0; i<4; i++)
     {
					n[i] = (ei.a[i] + ei.b[i]*n1.nx + ei.c[i]*n1.ny + ei.d[i]*n1.nz)/(6.0*ei.V);
					f[i] = (ei.b[i]*(n2.nx-n1.nx) + ei.c[i]*(n2.ny-n1.ny) + ei.d[i]*(n2.nz-n1.nz))/(6.0*ei.V);
					g[i] = (ei.b[i]*(n3.nx-n1.nx) + ei.c[i]*(n3.ny-n1.ny) + ei.d[i]*(n3.nz-n1.nz))/(6.0*ei.V);
					h[i] = (ei.b[i]*(n4.nx-n1.nx) + ei.c[i]*(n4.ny-n1.ny) + ei.d[i]*(n4.nz-n1.nz))/(6.0*ei.V);
					// cout<<endl<<n[i]<<" "<<f[i]<<" "<<g[i]<<" "<<h[i];
					//cout<<endl<<ei.a[i] <<" "<< n1.nx <<" "<< n1.ny <<" "<< n1.nz;
	 		}

		  double I[4][4], th[4][4]; //Defining intermediate variables I and th
		  for(i=0; i<4; i++)
	    {
					   for(j=0; j<4; j++)
					   {
							   th[i][j]=ei.b[i]*ei.b[j]+ei.c[i]*ei.c[j]+ei.d[i]*ei.d[j];

							   I[i][j] = J*( n[i]*n[j]/6.0 + (n[i]*(f[j]+g[j]+h[j]) + n[j]*(f[i]+g[i]+h[i]))/24.0    \
				                      + (f[i]*f[j] + g[i]*g[j] + h[i]*h[j])/60.0    \
									  					+ (f[i]*(g[j]+h[j]) + f[j]*(g[i]+h[i]) + g[i]*h[j]+g[j]*h[i])/120.0 );
							}
				}

      for(i=0; i<6; i++)
      {
					  for(j=0; j<6; j++)
					  {
						  i1=ei.num[i][0]; i2=ei.num[i][1]; j1=ei.num[j][0]; j2=ei.num[j][1];
							Q(i,j) = (I[i1][j1]*th[i2][j2] - I[i1][j2]*th[i2][j1] - I[i2][j1]*th[i1][j2] + I[i2][j2]*th[i1][j1])\
							 					* pow(2*M_PI/lambda,2.0) * (ei.l[i]*ei.l[j]/36.0) * ei.eleps/pow(ei.V,2.0);

							(*Amat).coeffRef(ei.edges[i],ei.edges[j]) -= Q(i,j);
							(*Qmat)(ei.edges[i],ei.edges[j])+=Q(i,j);

					 	}
	  	 }

}

void set_matrix_entries_RS(int e, vector<node> * nodelist, vector<element> * elementlist, vector<surf_element> * surf_elt_list, vector<edge> * edgelist,\
						double lambda, double theta, double phi, SparseMatrix<cmp,RowMajor> * Amat, SparseVector<cmp> * bvec, MatrixXcd *Rmat, MatrixXcd *Smat) // MatrixXcd * Amat, VectorXcd * bvec
{
		/**********************Calculating Matrix R**************************************************************/
		int i, j, i1, i2, j1, j2;
		MatrixXcd R(6,6), S(6,6); // Complex matrices R and S
		double th[4][4];

		surf_element es = (*surf_elt_list)[e];
		element ei = (*elementlist)[ es.glob_el_num[0] ];

		vector<node> surf_nd; node int_nd;	vector<double> nlvec;
		double nlx,nly,nlz, area;

		calc_normal(&es, nodelist, edgelist, &surf_nd, &int_nd, &nlvec, &area); // Computing normal vector for the surface element

		nlx=nlvec[0]; nly=nlvec[1]; nlz=nlvec[2];
		double p[4], q[4], r[4];
    for(i=0; i<4; i++)
    {
				p[i] = (ei.a[i] + ei.b[i]*surf_nd[0].nx + ei.c[i]*surf_nd[0].ny + ei.d[i]*surf_nd[0].nz)/(6.0*ei.V);
				q[i] = (ei.b[i]*(surf_nd[1].nx-surf_nd[0].nx) + ei.c[i]*(surf_nd[1].ny-surf_nd[0].ny) + ei.d[i]*(surf_nd[1].nz-surf_nd[0].nz))/(6.0*ei.V);
				r[i] = (ei.b[i]*(surf_nd[2].nx-surf_nd[0].nx) + ei.c[i]*(surf_nd[2].ny-surf_nd[0].ny) + ei.d[i]*(surf_nd[2].nz-surf_nd[0].nz))/(6.0*ei.V);
	   }

	  double B[4][4]; //Defining intermediate variable B (two edge integral)
		for(i=0; i<4; i++)
	  {
			  for(j=0; j<4; j++)
			  {
				  th[i][j]=ei.b[i]*ei.b[j]+ei.c[i]*ei.c[j]+ei.d[i]*ei.d[j];
				  B[i][j] = area*(p[i]*p[j] + (p[i]*(q[j]+r[j]) + p[j]*(q[i]+r[i]))/3.0 + (q[i]*q[j] + r[i]*r[j])/6.0 + (q[i]*r[j] + q[j]*r[i])/12.0);
			  }
	  }
		for(i=0; i<6; i++)
    {
			  for(j=0; j<6; j++)
			  {
						  i1=ei.num[i][0]; i2=ei.num[i][1]; j1=ei.num[j][0]; j2=ei.num[j][1];
							R(i,j) = cmp(0.0,1.0)*(2.0*M_PI/lambda)*(ei.l[i]*ei.l[j]/36.0)*sqrt(ei.eleps)/pow(ei.V,2.0)\
											*(B[i1][j1]*th[i2][j2] - B[i1][j2]*th[i2][j1] - B[i2][j1]*th[i1][j2] + B[i2][j2]*th[i1][j1]);
							(*Amat).coeffRef(ei.edges[i],ei.edges[j]) += R(i,j);
							(*Rmat)(ei.edges[i],ei.edges[j])+=R(i,j);
			   }
   	}

		/*****************************Calculating Matrix S**************************************************************************/
	  double F0, F1, F2, F3, psi[4], u[6], v[6], w[6], zeta[6][6];
	  for(i=0; i<4; i++)
		  psi[i] = ei.b[i]*nlx + ei.c[i]*nly + ei.d[i]*nlz;

	  for(i=0; i<6; i++)
	  {
		  		i1=ei.num[i][0]; i2=ei.num[i][1];
	        F0 = ei.a[i1]*psi[i2] - ei.a[i2]*psi[i1]; F1 = ei.b[i1]*psi[i2] - ei.b[i2]*psi[i1];
	        F2 = ei.c[i1]*psi[i2] - ei.c[i2]*psi[i1]; F3 = ei.d[i1]*psi[i2] - ei.d[i2]*psi[i1];

	        u[i] = F0 + F1*surf_nd[0].nx + F2*surf_nd[0].ny + F3*surf_nd[0].nz;
	        v[i] = F1*(surf_nd[1].nx-surf_nd[0].nx) + F2*(surf_nd[1].ny-surf_nd[0].ny) + F3*(surf_nd[1].nz-surf_nd[0].nz);
		  		w[i] = F1*(surf_nd[2].nx-surf_nd[0].nx) + F2*(surf_nd[2].ny-surf_nd[0].ny) + F3*(surf_nd[2].nz-surf_nd[0].nz);
		}

    for(i=0; i<6; i++)
    {
			  for(j=0; j<6; j++)
			  {
				  	zeta[i][j] = area*(u[i]*u[j] + (u[i]*(v[j]+w[j]) + u[j]*(v[i]+w[i]))/3.0 + (v[i]*v[j] + w[i]*w[j])/6.0 + (v[i]*w[j] + v[j]*w[i])/12.0);

						S(i,j) = cmp(0.0,1.0)*(2.0*M_PI/lambda)*(ei.l[i]*ei.l[j]/1296.0) * (sqrt(ei.eleps)/pow(ei.V,4.0)) * zeta[i][j];
						(*Amat).coeffRef(ei.edges[i],ei.edges[j]) -= S(i,j);
						(*Smat)(ei.edges[i],ei.edges[j]) += S(i,j);
				}
  	}

	/**********************Calculating Incident Field Vector b***********************************************************************/
	double K1, K4, tau, xg, yg, zg, xi[4];
	cmp G0, G1, G2, G3, K2, K3, K5, b[6], upsilon[4];

	for(i=0; i<4; i++)
	{
		  xi[i] = ei.b[i]*cos(theta)*cos(phi) + ei.c[i]*cos(theta)*sin(phi) - ei.d[i]*sin(theta);
			upsilon[i] = ei.b[i]*(sqrt(ei.eleps)*nlx + sin(theta)*cos(phi)) + \
					   ei.c[i]*(sqrt(ei.eleps)*nly + sin(theta)*sin(phi)) + ei.d[i]*(sqrt(ei.eleps)*nlz + cos(theta));
	}

	xg = ( surf_nd[0].nx + surf_nd[1].nx + surf_nd[2].nx )/3.0;
	yg = ( surf_nd[0].ny + surf_nd[1].ny + surf_nd[2].ny )/3.0;
	zg = ( surf_nd[0].nz + surf_nd[1].nz + surf_nd[2].nz )/3.0;

	K3 = sqrt(ei.eleps) + nlx*sin(theta)*cos(phi) + nly*sin(theta)*sin(phi) + nlz*cos(theta);
	K4 = nlx*cos(theta)*cos(phi) + nly*cos(theta)*sin(phi) - nlz*sin(theta);
	tau = xg*sin(theta)*cos(phi) + yg*sin(theta)*sin(phi) + zg*cos(theta);
	K5 = (2*M_PI/lambda)*( -sin(2*M_PI*tau/lambda) + cmp(0,1.0)*cos(2*M_PI*tau/lambda) );

	for(i=0; i<6; i++)
	{
		  i1=ei.num[i][0]; i2=ei.num[i][1];
      F0 = ei.a[i1]*xi[i2] - ei.a[i2]*xi[i1]; F1 = ei.b[i1]*xi[i2] - ei.b[i2]*xi[i1];
      F2 = ei.c[i1]*xi[i2] - ei.c[i2]*xi[i1]; F3 = ei.d[i1]*xi[i2] - ei.d[i2]*xi[i1];

      G0 = ei.a[i1]*upsilon[i2] - ei.a[i2]*upsilon[i1]; G1 = ei.b[i1]*upsilon[i2] - ei.b[i2]*upsilon[i1];
      G2 = ei.c[i1]*upsilon[i2] - ei.c[i2]*upsilon[i1]; G3 = ei.d[i1]*upsilon[i2] - ei.d[i2]*upsilon[i1];

			K1 = F0 + F1*xg + F2*yg + F3*zg; K2 = G0 + G1*xg + G2*yg + G3*zg;

      b[i] = (ei.l[i]/(36.0*pow(ei.V,2.0))) * area * K5 * ( K1*K3 - K2*K4 );

      (*bvec).coeffRef(ei.edges[i]) += b[i];
	}
}

//Intermediate function for RCS_mie()
void calc_psi_xi_mie(int n, cmp x, cmp *psi, cmp *xi, cmp *Dpsi, cmp *Dxi)
{
	cmp jn, h2n;
	jn = sp_bessel::sph_besselJ(n,x);
	h2n = sp_bessel::sph_hankelH2(n,x);
	(*psi) = x*jn;
	(*xi) = x*h2n;
	(*Dpsi) = 0.5*(jn + x*(sp_bessel::sph_besselJ(n-1,x) - sp_bessel::sph_besselJ(n+1,x)));
	(*Dxi) = 0.5*(h2n + x*(sp_bessel::sph_hankelH2(n-1,x) - sp_bessel::sph_hankelH2(n+1,x)));
}

//Mie series Electric field is correct, but curl of mie series electric field is inaccurate
void RCS_mie(double theta, double phi, double dist_fac, int N, double lambda, cmp epsi, double a, VectorXcd *E_mie, VectorXcd *v_mie, double *rcs)
{
	double k = (2*M_PI)/lambda;
	double r = dist_fac*a;

	cmp m = sqrt(epsi);
	double rho = k*r;

	double x = k*a;

	double mu = cos(theta);
	double pai[N+1], tau[N+1];
	pai[0]=0; tau[0]=0;
	pai[1]=1; tau[1]=mu;
	pai[2]=3*mu; tau[2]=6*pow(mu,2)-3;
	for (int n=3; n<=N; n++)
	{
		pai[n] = ( (2*n-1)*mu*pai[n-1] - n*pai[n-2] )/(n-1);
		tau[n] = n*mu*pai[n] - (n+1)*pai[n-1];
	}
	cmp sumr,sump,sumt;
	sumr=0; sump=0; sumt=0;
	cmp vr, vt, vp;
	vr = 0; vt = 0; vp = 0;
	cmp coeff;
	cmp a_n, b_n, H, dH;
	cmp psi_x, psi_mx, Dpsi_x, Dpsi_mx;
	cmp xi_x, xi_mx, Dxi_x, Dxi_mx;
	for (double n=1; n<=N; n++)
	{
		coeff = pow(cmp(0,-1),n)*(2*n+1)/(n*(n+1));

		//Coefficients a_n = A(n,x,m) and b_n = B(n,x,m), hankel function H = h1(n,rho), and dH = (1/rho)d(rho*h1(n,rho))/drho
		calc_psi_xi_mie(n, x, &psi_x, &xi_x, &Dpsi_x, &Dxi_x);
		calc_psi_xi_mie(n, m*x, &psi_mx, &xi_mx, &Dpsi_mx, &Dxi_mx);
		H = sp_bessel::sph_hankelH2(n,rho);
		dH = 0.5*(H/rho + sp_bessel::sph_hankelH2(n-1,rho) - sp_bessel::sph_hankelH2(n+1,rho));

		a_n = (m*psi_mx*Dpsi_x - psi_x*Dpsi_mx)/(m*psi_mx*Dxi_x - xi_x*Dpsi_mx);
		b_n = (psi_mx*Dpsi_x - m*psi_x*Dpsi_mx)/(psi_mx*Dxi_x - m*xi_x*Dpsi_mx);

		sumr += -coeff*cmp(0,1)*a_n*n*(n+1)*(H/rho)*sin(theta)*cos(phi)*pai[int(n)];
		sumt += -coeff*cos(phi)*( cmp(0,1)*a_n*dH*tau[int(n)] + b_n*H*pai[int(n)] );
		sump -= -coeff*sin(phi)*( cmp(0,1)*a_n*dH*pai[int(n)] + b_n*H*tau[int(n)] );

		vr += -k*coeff*sin(phi)*sin(theta)*n*(n+1)*b_n*pai[int(n)]*(H/rho);
		vt += -k*coeff*sin(phi)*( b_n*tau[int(n)]*dH - cmp(0,1)*a_n*pai[int(n)]*H );
		vp += -k*coeff*cos(phi)*( b_n*pai[int(n)]*dH - cmp(0,1)*a_n*tau[int(n)]*H );
	}

	VectorXcd E_pol(3); E_pol(0)=sumr; E_pol(1)=sumt; E_pol(2)=sump;
	VectorXcd E_cart(3);
	VectorXcd v_pol(3); v_pol(0)=vr; v_pol(1)=vt; v_pol(2)=vp;
	VectorXcd v_cart(3);
	Matrix3d U;
	U << sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta),
			 cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta),
			 -sin(phi), cos(phi), 0;

	E_cart = U.transpose()*E_pol; //Cartesian mie field
	v_cart = U.transpose()*v_pol; //Cartesian curl of mie field
	*E_mie = E_cart;
	*v_mie = v_cart;
	*rcs = 4*M_PI*pow(r,2)*(pow(abs(sumr),2) + pow(abs(sumt),2) + pow(abs(sump),2));
}

//Coordinates of centroid in both cartesian and polar are correct.
//Normal vector has correct direction and Magnitude, area is correct.
void set_surf_members(int e, vector<surf_element> * scat_surf_elt_list, vector<node> * nodelist, vector<element> * elementlist, \
																	vector<edge> * edgelist, VectorXcd *xvec, double lambda, cmp epsi, double a)
{
	surf_element es = (*scat_surf_elt_list)[e];

	int j1, j2, j3;
	j1 = es.face_nodes[0]; j2 = es.face_nodes[1]; j3 = es.face_nodes[2]; //Global numbers of shared nodes
	node n1, n2, n3;
	n1 = (*nodelist)[j1]; n2 = (*nodelist)[j2]; n3 = (*nodelist)[j3];

	int g1, g2;
	g1 = es.glob_el_num[0]; g2 = es.glob_el_num[1];

	double x,y,z; //Cartesian coordinates of centroid of exposed triangle
	x = (n1.nx + n2.nx + n3.nx)/3.0; y = (n1.ny + n2.ny + n3.ny)/3.0; z = (n1.nz + n2.nz + n3.nz)/3.0;
	(*scat_surf_elt_list)[e].x = x; (*scat_surf_elt_list)[e].y = y; (*scat_surf_elt_list)[e].z = z;

	double r,theta,phi; //Polar coordinates of centroid of exposed triangle
	r = sqrt(x*x + y*y + z*z);
	theta = acos(z/r);
	if (y >= 0)
		phi = atan2(y,x);
	else
		phi = 2*M_PI + atan2(y,x);

	(*scat_surf_elt_list)[e].r = r; (*scat_surf_elt_list)[e].theta = theta; (*scat_surf_elt_list)[e].phi = phi;

	/*Uncomment for verifying ONLY huygen's principle by saving mie fields at centroid of surf elements*/
	// VectorXcd E_mie(3), curl_mie(3);
	// double dist_fac = r/a; int N=15;
	// double rcs;
	// RCS_mie(theta, phi, dist_fac, N, lambda, epsi, a, &E_mie, &curl_mie, &rcs);
	//
	// (*scat_surf_elt_list)[e].Ex = E_mie(0); (*scat_surf_elt_list)[e].Ey = E_mie(1); (*scat_surf_elt_list)[e].Ez = E_mie(2);
	// (*scat_surf_elt_list)[e].vx = curl_mie(0); (*scat_surf_elt_list)[e].vy = curl_mie(1); (*scat_surf_elt_list)[e].vz = curl_mie(2);

	/****************************Updating unit outward normal vector in surf element class*****************************************/
	double e1_x, e1_y, e1_z, e2_x, e2_y, e2_z;
	double nlx, nly, nlz;

	e1_x = n2.nx - n1.nx;	e1_y = n2.ny - n1.ny;	e1_z = n2.nz - n1.nz;
	e2_x = n3.nx - n1.nx; e2_y = n3.ny - n1.ny; e2_z = n3.nz - n1.nz;

	nlx = e1_y*e2_z - e1_z*e2_y; nly = e1_z*e2_x - e1_x*e2_z; nlz = e1_x*e2_y - e1_y*e2_x;
	double nl_mag = sqrt(nlx*nlx + nly*nly + nlz*nlz);
	nlx = nlx/nl_mag; nly = nly/nl_mag; nlz = nlz/nl_mag;
	double area = 0.5*nl_mag;
	(*scat_surf_elt_list)[e].area = area;		//Initializes area of the exposed triangle in surf_element class
	//cout<<area<<endl;
	if (nlx*x + nly*y + nlz*z<0) 	// If normal points inward...
	{nlx = -nlx; nly = -nly; nlz = -nlz;}   // Reverse the direction of normal

	(*scat_surf_elt_list)[e].nlx = nlx; (*scat_surf_elt_list)[e].nly = nly; (*scat_surf_elt_list)[e].nlz = nlz;
	//cout<<nlx<<" "<<nly<<" "<<nlz<<endl<<endl;
	Vector3d nhat(nlx,nly,nlz);

	/*******************Now compute electric field at (x,y,z) as (u1*T1 + ... u9*T9) and its curl***************************************/
	element el1, el2;
	el1 = (*elementlist)[g1]; el2 = (*elementlist)[g2];
	double L1, L2, LL1, LL2; Vector3d Tvec, TTvec;
	cmp coeff, ccoeff; int ii1, ii2; 	int i1,i2;

	cmp Ex, Ey, Ez; cmp vx, vy, vz; //(Ex,Ey,Ez): electric field at centroid & (vx,vy,vz):curl of electric field at centroid
	Ex=0; Ey=0; Ez=0; vx=0; vy=0; vz=0;
	int nii, njj, Nii, Njj;
	double Tdn, Cdn;

	Vector3d Ttan, Tperp, Ctan, Cperp;
	for (int i=0; i<6; i++)
	{
				coeff = (*xvec)[el1.edges[i]]; ccoeff = (*xvec)[el2.edges[i]];

				i1 = el1.num[i][0]; i2 = el1.num[i][1];
				ii1 = el2.num[i][0]; ii2 = el2.num[i][1];

				L1 = (el1.a[i1] + el1.b[i1]*x + el1.c[i1]*y + el1.d[i1]*z)/(6.0*el1.V);
				L2 = (el1.a[i2] + el1.b[i2]*x + el1.c[i2]*y + el1.d[i2]*z)/(6.0*el1.V);
				LL1 = (el2.a[ii1] + el2.b[ii1]*x + el2.c[ii1]*y + el2.d[ii1]*z)/(6.0*el2.V);
				LL2 = (el2.a[ii2] + el2.b[ii2]*x + el2.c[ii2]*y + el2.d[ii2]*z)/(6.0*el2.V);

				Vector3d grad_L1(el1.b[i1],el1.c[i1],el1.d[i1]), grad_L2(el1.b[i2],el1.c[i2],el1.d[i2]);
				Vector3d grad_LL1(el2.b[ii1],el2.c[ii1],el2.d[ii1]), grad_LL2(el2.b[ii2],el2.c[ii2],el2.d[ii2]);

				grad_L1 *= 1/(6.0*el1.V); grad_L2 *= 1/(6.0*el1.V);
				grad_LL1 *= 1/(6.0*el2.V); grad_LL2 *= 1/(6.0*el2.V);

				Tvec = el1.l[i]*(L1*grad_L2 - L2*grad_L1);	TTvec = el2.l[i]*(LL1*grad_LL2 - LL2*grad_LL1);

				Vector3d curl1((el1.c[i1]*el1.d[i2] - el1.c[i2]*el1.d[i1]), \
											 (el1.d[i1]*el1.b[i2] - el1.d[i2]*el1.b[i1]), \
											 (el1.b[i1]*el1.c[i2] - el1.b[i2]*el1.c[i1]));
				Vector3d curl2((el2.c[ii1]*el2.d[ii2] - el2.c[ii2]*el2.d[ii1]), \
				 							 (el2.d[ii1]*el2.b[ii2] - el2.d[ii2]*el2.b[ii1]), \
											 (el2.b[ii1]*el2.c[ii2] - el2.b[ii2]*el2.c[ii1]));

				curl1 *= el1.l[i]/(18*pow(el1.V,2));	curl2 *= el2.l[i]/(18*pow(el2.V,2));

				nii = (*edgelist)[el1.edges[i]].ni; njj = (*edgelist)[el1.edges[i]].nj;
				Nii = (*edgelist)[el2.edges[i]].ni; Njj = (*edgelist)[el2.edges[i]].nj;

				// If we find a shared edge: average the tangential component, sum the normal components (of Tvec/TTvec and curl1/curl2)
				// If we find a non-shared edge: sum directly
				// To determine whether shared or not, compare nii and njj with j1/j2/j3
				if ((nii==j1||nii==j2||nii==j3) && (njj==j1||njj==j2||njj==j3)) //check if el1's edge is a shared edge
				{
					Tdn = Tvec(0)*nlx + Tvec(1)*nly + Tvec(2)*nlz;
					Cdn = curl1(0)*nlx + curl1(1)*nly + curl1(2)*nlz;
					Ttan = Tvec - (Tdn)*nhat;
					Tperp = Tvec - Ttan;
					Ctan = curl1 - (Cdn)*nhat;
					Cperp = curl1 - Ctan;
					Ex+=coeff*(0.5*Ttan(0) + Tperp(0)); Ey+=coeff*(0.5*Ttan(1) + Tperp(1)); Ez+=coeff*(0.5*Ttan(2) + Tperp(2));
					vx+=coeff*(0.5*Ctan(0) + Cperp(0)); vy+=coeff*(0.5*Ctan(1) + Cperp(1)); vz+=coeff*(0.5*Ctan(2) + Cperp(2));
				}
				else //non-shared implies sum directly
				{
					Ex+=coeff*Tvec(0); Ey+=coeff*Tvec(1); Ez+=coeff*Tvec(2);
					vx+=coeff*curl1(0); vy+=coeff*curl1(1); vz+=coeff*curl1(2);
				}

				if ((Nii==j1||Nii==j2||Nii==j3) && (Njj==j1||Njj==j2||Njj==j3)) //check if el2's edge is a shared edge
				{
					Tdn = TTvec(0)*nlx + TTvec(1)*nly + TTvec(2)*nlz;
					Cdn = curl2(0)*nlx + curl2(1)*nly + curl2(2)*nlz;
					Ttan = TTvec - (Tdn)*nhat;
					Tperp = TTvec - Ttan;
					Ctan = curl2 - (Cdn)*nhat;
					Cperp = curl2 - Ctan;
					Ex+=ccoeff*(0.5*Ttan(0) + Tperp(0)); Ey+=ccoeff*(0.5*Ttan(1) + Tperp(1)); Ez+=ccoeff*(0.5*Ttan(2) + Tperp(2));
					vx+=ccoeff*(0.5*Ctan(0) + Cperp(0)); vy+=ccoeff*(0.5*Ctan(1) + Cperp(1)); vz+=ccoeff*(0.5*Ctan(2) + Cperp(2));
				}
				else //non-shared implies sum directly
				{
					Ex+=ccoeff*TTvec(0); Ey+=ccoeff*TTvec(1); Ez+=ccoeff*TTvec(2);
					vx+=ccoeff*curl2(0); vy+=ccoeff*curl2(1); vz+=ccoeff*curl2(2);
				}
	}

	(*scat_surf_elt_list)[e].Ex = Ex; (*scat_surf_elt_list)[e].Ey = Ey; (*scat_surf_elt_list)[e].Ez = Ez; //Electric field
	(*scat_surf_elt_list)[e].vx = vx; (*scat_surf_elt_list)[e].vy = vy; (*scat_surf_elt_list)[e].vz = vz; //Curl of E field

}

void RCS_huygen(int e, vector<surf_element> *scat_surf_elt_list, double rp, Vector3d rp_vec, double lambda, VectorXcd *u, VectorXcd *w)
{
	// This function RCS_huygen() is called in a loop over surface elements in fem3d.cpp
	// "u" and "w" passed to this function are intermediate complex 3d vectors; both initialized to (0,0,0) in fem3d.cpp but updated here
	// At the end of the loop in fem3d.cpp, complex vectors "u" and "w" have the required sum

	surf_element es = (*scat_surf_elt_list)[e];		// Surface element object named "es"

	cmp Ex, Ey, Ez;	//Field values at centroid of exposed triangle
	cmp vx, vy, vz; //Curl of electric field at the centroid

	Ex = (*scat_surf_elt_list)[e].Ex; Ey = (*scat_surf_elt_list)[e].Ey; Ez = (*scat_surf_elt_list)[e].Ez;
	vx = (*scat_surf_elt_list)[e].vx; vy = (*scat_surf_elt_list)[e].vy; vz = (*scat_surf_elt_list)[e].vz;

	double k = 2*M_PI/lambda;											// Wavevector of incident wave

	double x,y,z;
	x = es.x; y = es.y; z = es.z;										// Cartesian coordinates of near field pt (centroid)
	Vector3d r_vec{x, y, z};

	double nlx, nly, nlz;
	nlx = es.nlx; nly = es.nly; nlz=es.nlz;					// Unit outward normal
	double area = es.area;													// Area of exposed triangle

	double xp,yp,zp;
	xp = rp_vec(0); yp = rp_vec(1); zp = rp_vec(2);

	cmp expo = exp(-cmp(0,1)*k*(r_vec - rp_vec).norm());	//intermediate variable expo

	(*w)(0) += (nly*vz - nlz*vy)*expo*area;
	(*w)(1) += (nlz*vx - nlx*vz)*expo*area;
	(*w)(2) += (nlx*vy - nly*vx)*expo*area;

	cmp tx, ty, tz;
	tx = nly*Ez - nlz*Ey; ty = nlz*Ex - nlx*Ez; tz = nlx*Ey - nly*Ex;

	(*u)(0) += -((yp - y)*tz - (zp - z)*ty)*expo*area;
	(*u)(1) += -((zp - z)*tx - (xp - x)*tz)*expo*area;
	(*u)(2) += -((xp - x)*ty - (yp - y)*tx)*expo*area;
}

void manuf_at_pt(int e, vector<element> * elementlist, vector<double> * Emanfx, vector<double> * Emanfy, vector<double> * Emanfz,\
								vector<double> * E_ls_x, vector<double> * E_ls_y, vector<double> * E_ls_z, \
								VectorXcd *xvec_tilde, vector<double> * xforced, double x, double y, double z, int ii)
{
	element ei = (*elementlist)[e];
	double coeff_manf; double L1, L2; int i,i1,i2;
	Vector3d Tvec;
	double coeff_LS;

	for(i=0;i<6;i++)
	{
		coeff_LS = (*xforced)[ei.edges[i]];
		coeff_manf = real((*xvec_tilde)[ei.edges[i]]);

		i1 = ei.num[i][0]; i2 = ei.num[i][1];
		L1 = (ei.a[i1] + ei.b[i1]*x + ei.c[i1]*y + ei.d[i1]*z)/(6.0*ei.V);
		L2 = (ei.a[i2] + ei.b[i2]*x + ei.c[i2]*y + ei.d[i2]*z)/(6.0*ei.V);

		Vector3d grad_L1(ei.b[i1],ei.c[i1],ei.d[i1]), grad_L2(ei.b[i2],ei.c[i2],ei.d[i2]);

		grad_L1 *= 1/(6.0*ei.V); grad_L2 *= 1/(6.0*ei.V);

		Tvec = ei.l[i]*(L1*grad_L2 - L2*grad_L1);

		(*Emanfx)[ii] += coeff_manf*Tvec(0);
		(*Emanfy)[ii] += coeff_manf*Tvec(1);
		(*Emanfz)[ii] += coeff_manf*Tvec(2);

		(*E_ls_x)[ii] += coeff_LS*Tvec(0);
		(*E_ls_y)[ii] += coeff_LS*Tvec(1);
		(*E_ls_z)[ii] += coeff_LS*Tvec(2);
	}

}
void manufactured_soln_PQ(int e, vector<node> * nodelist, vector<element> * elementlist, vector<edge> * edgelist, double lambda, SparseVector<cmp> * b_tilde)
{
	   element ei = (*elementlist)[e];  // e^th element

	   int i, i1, i2;

	   VectorXcd P(6), Q(6); //Complex vectors P and Q

		/**************************************Calculating vectors P and Q************************************************/

	    Matrix3d J;	  //Jacobian Matrix
	    for(i=0; i<3; i++)
	    {
			J(0,i) = (*nodelist)[ei.nodes[i+1]].nx - (*nodelist)[ei.nodes[0]].nx;
			J(1,i) = (*nodelist)[ei.nodes[i+1]].ny - (*nodelist)[ei.nodes[0]].ny;
			J(2,i) = (*nodelist)[ei.nodes[i+1]].nz - (*nodelist)[ei.nodes[0]].nz;
	    }
	    double det_J = J.determinant();  	// Jacobian determinant
	    node n1 = (*nodelist)[ei.nodes[0]]; // node object of point r1

			Vector3d alpha; Matrix3d A;
	    double Q1a, Q1b, Q2a, Q2b, Q3a, Q3b;

		for(i=0;i<6;i++)
		{
			i1=ei.num[i][0]; i2=ei.num[i][1];

			P(i) = (ei.l[i]/(18.0*ei.V))*((ei.c[i1]*ei.d[i2] - ei.c[i2]*ei.d[i1]) + (ei.d[i1]*ei.b[i2] - ei.d[i2]*ei.b[i1]) + (ei.b[i1]*ei.c[i2] - ei.b[i2]*ei.c[i1]));

			(*b_tilde).coeffRef(ei.edges[i]) += P(i);

			alpha(0) = ei.a[i1]*ei.c[i2] - ei.a[i2]*ei.c[i1];
			alpha(1) = ei.a[i1]*ei.d[i2] - ei.a[i2]*ei.d[i1];
			alpha(2) = ei.a[i1]*ei.b[i2] - ei.a[i2]*ei.b[i1];

			A << ei.b[i1]*ei.c[i2]-ei.b[i2]*ei.c[i1], ei.b[i1]*ei.d[i2]-ei.b[i2]*ei.d[i1], ei.c[i1]*ei.b[i2]-ei.c[i2]*ei.b[i1],
                        0, ei.c[i1]*ei.d[i2]-ei.c[i2]*ei.d[i1], ei.d[i1]*ei.c[i2]-ei.d[i2]*ei.c[i1],
                        0, 0, ei.d[i1]*ei.b[i2]-ei.d[i2]*ei.b[i1];
			Vector3d one(1.0,1.0,1.0);

			Vector3d r1(n1.nx,n1.ny,n1.nz);

			Q1a=alpha.transpose()*r1;
			Q1b=r1.transpose()*A*r1;
			Q2a=alpha.transpose()*J*one;
			Q2b=r1.transpose()*(A+A.transpose())*J*one;
			Q3a=one.transpose()*J.transpose()*A*J*one;
			Q3b=(J.transpose()*A*J).trace();

			Q(i) = (det_J*pow(2*M_PI/lambda,2)*ei.eleps*ei.l[i]/(36.0*pow(ei.V,2.0)))*((Q1a+Q1b)/6.0 + (Q2a+Q2b)/24.0 + (Q3a+Q3b)/120.0);

			//cout<<endl<<"Q("<<i<<") = "<<Q(i)<<" and global index is "<<ei.edges[i]<<endl;
		  (*b_tilde).coeffRef(ei.edges[i]) -= Q(i);
		}

}

void manufactured_soln_RS(int e, vector<node> * nodelist, vector<element> * elementlist, vector<surf_element> * surf_elt_list, vector<edge> * edgelist, double lambda, SparseVector<cmp> * b_tilde)
{
	int i, i1, i2;
	VectorXcd R(6), S(6); 	// Complex vectors R and S

	// Getting the surf_element class and element class objects using the integer "e"
	surf_element es = (*surf_elt_list)[e];
	element ei = (*elementlist)[ es.glob_el_num[0] ];

	// For a typical surface element...
	vector<node> surf_nd; 			// surf_nd stores its three surface node
	node int_nd;					// int_nd stores its only interior node
	vector<double> nlvec;			// nlvec stores the outward unit normal vector of the exposed triangle
	double nlx,nly,nlz, area;		// nlx, nly and nlz are components of nlvec. And "area" stores the area of the exposed triangle

	calc_normal(&es, nodelist, edgelist, &surf_nd, &int_nd, &nlvec, &area); // Initializing nlvec, surf_nd, int_nd and area
	nlx=nlvec[0]; nly=nlvec[1]; nlz=nlvec[2];

   node n1=surf_nd[0],n2=surf_nd[1],n3=surf_nd[2];
   //~ //cout<<endl<<n1.nx<<" "<<n1.ny<<" "<<n1.nz;
   //~ //cout<<endl<<n2.nx<<" "<<n2.ny<<" "<<n2.nz;
   //~ //cout<<endl<<n3.nx<<" "<<n3.ny<<" "<<n3.nz;
   Vector3d alpha, one(1,1,1), n_tilde(nly,nlz,nlx), h, r1(n1.nx,n1.ny,n1.nz),u(n2.nx-n1.nx,n2.ny-n1.ny,n2.nz-n1.nz), v(n3.nx-n1.nx,n3.ny-n1.ny,n3.nz-n1.nz);
   double R1a,R1b,R2a,R2b,R3a,R3b,R3c,S1a,S1b,S2a,S2b,S3a,S3b,S3c,G,psi[4];
   Matrix3d A,M;

	 //cout<<endl<<v<<endl;

   for(int i=0;i<4;i++)
	   psi[i]=ei.b[i]*nlx+ei.c[i]*nly+ei.d[i]*nlz;

   for(i=0; i<6; i++)
   {
	  i1=ei.num[i][0]; i2=ei.num[i][1];
	  alpha(0) = ei.a[i1]*ei.c[i2] - ei.a[i2]*ei.c[i1];
	  alpha(1) = ei.a[i1]*ei.d[i2] - ei.a[i2]*ei.d[i1];
	  alpha(2) = ei.a[i1]*ei.b[i2] - ei.a[i2]*ei.b[i1];
	  //cout<<endl<<ei.a[i1]*ei.c[i2] - ei.a[i2]*ei.c[i1];
	  A << (ei.b[i1]*ei.c[i2]-ei.b[i2]*ei.c[i1]), ei.b[i1]*ei.d[i2], ei.c[i1]*ei.b[i2],
					-ei.b[i2]*ei.d[i1], (ei.c[i1]*ei.d[i2]-ei.c[i2]*ei.d[i1]), ei.d[i1]*ei.c[i2],
					-ei.c[i2]*ei.b[i1], -ei.d[i2]*ei.c[i1], (ei.d[i1]*ei.b[i2]-ei.d[i2]*ei.b[i1]);

	  //cout<<endl<<A<<endl;

	  R1a=alpha.transpose()*r1; R1b=r1.transpose()*A*r1;
	  R2a=r1.transpose()*(A + A.transpose())*(u+v);  R2b=alpha.transpose()*(u+v);  //Issue with R2a
	  R3a=u.transpose()*A*u; R3b=v.transpose()*A*v; R3c=0.5*u.transpose()*(A + A.transpose())*v;  //Issue with R3c

	  R(i) = cmp(0.0,1.0)*(2.0*M_PI/lambda)*sqrt(ei.eleps)*(area*ei.l[i]/(18.0*pow(ei.V,2.0)))*((R1a+R1b)/2.0 + (R2a+R2b)/6.0 + (R3a+R3b+R3c)/12.0);

		//cout<<endl<<"R("<<i<<") = "<<R(i)<<" and global index is "<<ei.edges[i]<<endl;
	  (*b_tilde).coeffRef(ei.edges[i]) += R(i);

	  /*****************Calculation of vector - S********************/
	  h(0) = ei.b[i1]*psi[i2] - ei.b[i2]*psi[i1];
	  h(1) = ei.c[i1]*psi[i2] - ei.c[i2]*psi[i1];
	  h(2) = ei.d[i1]*psi[i2] - ei.d[i2]*psi[i1];
	  G = ei.a[i1]*psi[i2] - ei.a[i2]*psi[i1];
	  M=h*n_tilde.transpose();

	  S1a=G*n_tilde.transpose()*r1; S1b=r1.transpose()*M*r1;
	  S2a=r1.transpose()*(M + M.transpose())*(u+v);  S2b=G*n_tilde.transpose()*(u+v);
	  S3a=u.transpose()*M*u; S3b=v.transpose()*M*v; S3c=0.5*u.transpose()*(M + M.transpose())*v;

	  S(i) = cmp(0.0,1.0)*(2.0*M_PI/lambda)*sqrt(ei.eleps)*(area*ei.l[i]/(18.0*pow(ei.V,2.0)))*((S1a+S1b)/2.0 + (S2a+S2b)/6.0 + (S3a+S3b+S3c)/12.0);

		//cout<<endl<<"S("<<i<<") = "<<S(i)<<" and global index is "<<ei.edges[i]<<endl;
	  (*b_tilde).coeffRef(ei.edges[i]) -= S(i);
   }
}

void few_element_mms(vector<node> * nodelist,vector<element> * elementlist,\
	vector<surf_element> * surf_elt_list,vector<edge> * edgelist,SparseMatrix<cmp,RowMajor> * Amat,double lambda)
{
	/******Method of Manufactured Solutions*****************************************************/
	int No=signed((*edgelist).size());
	Eigen::SparseVector<cmp> b_tilde(No); b_tilde.setZero(); // Sparse vector bvec for incident field

	for(int e=0; e < signed((*elementlist).size()); e++)   // loop over the entire element list
		manufactured_soln_PQ(e, nodelist, elementlist, edgelist, lambda, &b_tilde);

	for(int e=0; e < signed((*surf_elt_list).size()); e++)  // loop over the surface element list
		manufactured_soln_RS(e, nodelist, elementlist, surf_elt_list, edgelist, lambda, &b_tilde);

	cout<<endl<<endl<<"b tilde:\n"<<VectorXcd(b_tilde);

	SparseLU<SparseMatrix<cmp>, COLAMDOrdering<int>> solver2; // LU decomposition
	solver2.analyzePattern(*Amat);                             // Compute the ordering permutation vector from the structural pattern of A
	solver2.factorize(*Amat);                                  // Compute the numerical factorization
	Eigen::VectorXcd xvec_tilde = solver2.solve(b_tilde);     // Use the factors to solve the linear system.

	cout<<endl<<endl<<"x tilde:\n"<<VectorXcd(xvec_tilde);

	/*****Important: The below line is to be edited for points in a single or two element mesh************/
	fstream fin; fin.open("three_elt_pts.csv", ios::in);	string line;

	/************************************************/
	std::vector<double> x(207),y(207),z(207),points(207),mag_true(207); char c; int i=0;

	while (fin && i<207)
	{
		getline(fin, line);
		stringstream s(line);
		s>>x.at(i)>>c>>y.at(i)>>c>>z.at(i);
		mag_true.at(i) = sqrt(x.at(i)*x.at(i) + y.at(i)*y.at(i) + z.at(i)*z.at(i));
		i++;
	}

	for (i=0;i<207;i++)
		points.at(i) = i+1;

	std::vector<double> Emanfx(207), Emanfy(207), Emanfz(207), mag_mms(207);
	std::vector<double> E_ls_x(207), E_ls_y(207), E_ls_z(207), mag_ls(207);
	//std::vector<double> err_x(207), err_y(207), err_z(207), err_mag(207);

	std::vector<double> xforced0{2.8284,-1.732,-3.8636,1.732,1.0352,1.732};
	std::vector<double> xforced1{-2.8284,-1.0352,-1.732,1.732,3.8636,1.732};
	std::vector<double> xforced2{0,2.8284,-1.732,2.8284,-3.8636,-1.732};

	std::fill(Emanfx.begin(), Emanfx.end(), 0);std::fill(E_ls_x.begin(), E_ls_x.end(), 0);
	std::fill(Emanfy.begin(), Emanfy.end(), 0);std::fill(E_ls_y.begin(), E_ls_y.end(), 0);
	std::fill(Emanfz.begin(), Emanfz.end(), 0);std::fill(E_ls_z.begin(), E_ls_z.end(), 0);

	for (int ii=0;ii<207;ii++)
	{
		if(z.at(ii)==0)
		{
			manuf_at_pt(0, elementlist, &Emanfx, &Emanfy, &Emanfz, &E_ls_x, &E_ls_y, &E_ls_z, &xvec_tilde, &xforced0, x.at(ii), y.at(ii), z.at(ii), ii);
			manuf_at_pt(1, elementlist, &Emanfx, &Emanfy, &Emanfz, &E_ls_x, &E_ls_y, &E_ls_z, &xvec_tilde, &xforced1, x.at(ii), y.at(ii), z.at(ii), ii);
		}
		else if(z.at(ii)>0)
		{
			if(1.2247*x.at(ii) + 2.121*y.at(ii) + 0.866*z.at(ii) - 1.2247 < 0)
			{manuf_at_pt(1, elementlist, &Emanfx, &Emanfy, &Emanfz, &E_ls_x, &E_ls_y, &E_ls_z, &xvec_tilde, &xforced1, x.at(ii), y.at(ii), z.at(ii), ii);}
			else if(1.2247*x.at(ii) + 2.121*y.at(ii) + 0.866*z.at(ii) - 1.2247 > 0)
			{manuf_at_pt(2, elementlist, &Emanfx, &Emanfy, &Emanfz, &E_ls_x, &E_ls_y, &E_ls_z, &xvec_tilde, &xforced2, x.at(ii), y.at(ii), z.at(ii), ii);}
		}
		else if(z.at(ii)<0)
		{manuf_at_pt(0, elementlist, &Emanfx, &Emanfy, &Emanfz, &E_ls_x, &E_ls_y, &E_ls_z, &xvec_tilde, &xforced0, x.at(ii), y.at(ii), z.at(ii), ii);}

		mag_mms.at(ii) = sqrt(Emanfx.at(ii)*Emanfx.at(ii) + Emanfy.at(ii)*Emanfy.at(ii) + Emanfz.at(ii)*Emanfz.at(ii));
		mag_ls.at(ii) = sqrt(E_ls_x.at(ii)*E_ls_x.at(ii) + E_ls_y.at(ii)*E_ls_y.at(ii) + E_ls_z.at(ii)*E_ls_z.at(ii));
		//err_x.at(ii) = abs(Emanfx.at(ii)-z.at(ii)); err_y.at(ii) = abs(Emanfy.at(ii)-x.at(ii));
		//err_z.at(ii) = abs(Emanfz.at(ii)-y.at(ii));	err_mag.at(ii) = abs(mag_mms.at(ii)-mag_true.at(ii));
	}

	/*******Plotting true vs manufactured field***********************/
	plt::suptitle("True vs Manufactured Field");
	plt::subplot(2, 2, 1);
	plt::named_plot("True",points,z);
	//plt::named_plot("Least Squares",points,E_ls_x);
	plt::named_plot("Manufactured",points,Emanfx);
	plt::xlabel("Grid points");
	plt::ylabel("Re(E_x)");
	plt::legend();

	plt::subplot(2, 2, 2);
	plt::named_plot("True",points,x);
	//plt::named_plot("Least Squares",points,E_ls_y);
	plt::named_plot("Manufactured",points,Emanfy);
	plt::xlabel("Grid points");
	plt::ylabel("Re(E_y)");
	plt::legend();

	plt::subplot(2, 2, 3);
	plt::named_plot("True",points,y);
	//plt::named_plot("Least Squares",points,E_ls_z);
	plt::named_plot("Manufactured",points,Emanfz);
	plt::xlabel("Grid points");
	plt::ylabel("Re(E_z)");
	plt::legend();

	plt::subplot(2, 2, 4);
	plt::named_plot("True",points,mag_true);
	//plt::named_plot("Least Squares",points,mag_ls);
	plt::named_plot("Manufactured",points,mag_mms);
	plt::xlabel("Grid points");
	plt::ylabel("Magnitude");
	plt::legend();

	plt::show();
/*******************************************************************************************************************************************/
}

void mie_vs_fem_nearfield(cmp epsi,vector<surf_element> *scat_surf_elt_list, double lambda, double a)
{
	 int N=15; //Number of terms in mie series

	 VectorXcd E_mie(3), v_mie(3); double rcss;
	 int huy_size = signed((*scat_surf_elt_list).size());
	 //cout<<endl<<huy_size<<" "<<signed(surf_elt_list.size())<<endl;

	 double thet_near, phi_near, dista_fac;
	 vector<double> points(huy_size);

	 /*********Mie Series near field ***********/
	 vector<double> Exx_mie(huy_size), Eyy_mie(huy_size), Ezz_mie(huy_size), rcss_mie(huy_size);
	 for (int i=0; i<huy_size; i++)
	 {
		 points.at(i) = i;
		 surf_element es = (*scat_surf_elt_list)[i];
		 dista_fac = es.r/a; thet_near = es.theta; phi_near = es.phi;
		 RCS_mie(thet_near, phi_near, dista_fac, N, lambda, epsi, a, &E_mie, &v_mie, &rcss);
		 Exx_mie.at(i) = abs(E_mie(0)); Eyy_mie.at(i) = abs(E_mie(1)); Ezz_mie.at(i) = abs(E_mie(2)); //Be careful here. abs or real part or im part?
		 rcss_mie.at(i) = 10*log10(rcss/(lambda*lambda));
	 }

	 /*****************FEM near field*********************/
	 vector<double> Ex_fem(huy_size), Ey_fem(huy_size), Ez_fem(huy_size), rcs_fem(huy_size);
	 for (int i=0; i<huy_size; i++)
	 {
		 surf_element es = (*scat_surf_elt_list)[i];
		 Ex_fem.at(i) = abs(es.Ex); Ey_fem.at(i) = abs(es.Ey); Ez_fem.at(i) = abs(es.Ez); //Be careful here. abs or real part or im part?
		 rcss = 4*M_PI*pow(es.r,2)*( pow(abs(es.Ex),2) + pow(abs(es.Ey),2) + pow(abs(es.Ez),2) );
		 rcs_fem.at(i) = 10*log10(rcss/(lambda*lambda));
	 }

	 /*****Plotting Mie vs FEM in Near Field*****/
	 plt::suptitle("Mie vs FEM at near-field");
	 plt::subplot(2,2,1);
	 plt::named_plot("Mie near-field",points,Exx_mie);
	 plt::named_plot("FEM near-field",points,Ex_fem);
	 plt::xlabel("Points on surface"); plt::ylabel("Ex"); plt::legend();

	 plt::subplot(2,2,2);
	 plt::named_plot("Mie near-field",points,Eyy_mie);
	 plt::named_plot("FEM near-field",points,Ey_fem);
	 plt::xlabel("Points on surface"); plt::ylabel("Ey"); plt::legend();

	 plt::subplot(2,2,3);
	 plt::named_plot("Mie near-field",points,Ezz_mie);
	 plt::named_plot("FEM near-field",points,Ez_fem);
	 plt::xlabel("Points on surface"); plt::ylabel("Ez"); plt::legend();

	 plt::subplot(2,2,4);
	 plt::named_plot("Mie near-field",points,rcss_mie);
	 plt::named_plot("FEM near-field",points,rcs_fem);
	 plt::xlabel("Points on surface"); plt::ylabel("10*log(rcs/lam^2)"); plt::legend();

	 plt::show();
}

void mie_vs_fem_farfield(cmp epsi,vector<surf_element> *scat_surf_elt_list,double lambda,double a,double dist_fac,double thetp,double phip)
{
	std::vector<double> thetaas(78);		//far-field angles
	int N=15; 													//number of mie series terms
	VectorXcd E_mie(3), v_mie(3); 			//electric field and its curl
	double rcss;

	/**************Mie series far field **************/
	vector<double> Ex_mie(78), Ey_mie(78), Ez_mie(78), rcs_mie(78);
	for(int i=0;i<78;i++)
	{
		thetp+=0.04;
		RCS_mie(thetp, phip, dist_fac, N, lambda, epsi, a, &E_mie, &v_mie, &rcss);
		rcs_mie.at(i) = 10*log10(rcss/(lambda*lambda));
		Ex_mie.at(i) = abs(E_mie(0)); Ey_mie.at(i) = abs(E_mie(1)); Ez_mie.at(i) = abs(E_mie(2));
		thetaas.at(i) = thetp*180/M_PI;
	}

	/*****************FEM far field*********************/
	// Spherical coordinates of far-field point: rp and phip are fixed, and thetp is varied
	double rp = dist_fac*a;

	double xp,yp,zp;
	double k = 2*M_PI/lambda;
	VectorXcd Escat(3); std::vector<double> rcs_huygen(78); double rcs;
	vector<double> Ex_huygen(78), Ey_huygen(78), Ez_huygen(78);

	thetp=0;
	for (int i=0; i<78; i++)									// Loop over farfield points
	{
	 thetp += 0.04;
	 xp=rp*sin(thetp)*cos(phip);
	 yp=rp*sin(thetp)*sin(phip);
	 zp=rp*cos(thetp);
	 Vector3d rp_vec{xp,yp,zp};								// Farfield pt as a 3d cartesian vector
	 MatrixXcd M(3,3);												// Intermediate matrix variable M (purely a function of the far field point)
	 M << (yp*yp + zp*zp), -xp*yp, -xp*zp,
				-xp*yp, (xp*xp + zp*zp), -yp*zp,
				-xp*zp, -yp*zp, (xp*xp + yp*yp);

	 VectorXcd u(3), w(3); 										// These vectors are filled up in RCS_huygen()
	 u(0) = cmp(0,0); u(1) = cmp(0,0); u(2) = cmp(0,0);
	 w(0) = cmp(0,0); w(1) = cmp(0,0); w(2) = cmp(0,0);

	 for (int e=0; e<signed((*scat_surf_elt_list).size()); e++)				// Loop over all nearfield pts (for each farfield pt)
		 RCS_huygen(e, scat_surf_elt_list, rp, rp_vec, lambda, &u, &w); // Populates vectors u and w necessary for farfield expression

	 Escat = ( M*w/rp + cmp(0,1)*k*u )/(4*M_PI*rp*rp); 								// Farfield electric field expression using Huygen's principle

	 Ex_huygen.at(i) = abs(Escat(0)); Ey_huygen.at(i) = abs(Escat(1)); Ez_huygen.at(i) = abs(Escat(2));

	 rcs = 4*M_PI*rp*rp*(pow(abs(Escat(0)),2) + pow(abs(Escat(1)),2) + pow(abs(Escat(2)),2));

	 rcs_huygen.at(i) = 10*log10(rcs/(lambda*lambda));
	}

	/**********Plotting Mie vs propogated FEM at far-field**************/
	plt::suptitle("Mie vs FEM at far-field");
	plt::subplot(2,2,1);
	plt::named_plot("Mie at far-field",thetaas,Ex_mie);
	plt::named_plot("FEM prop to far-field",thetaas,Ex_huygen);
	plt::xlabel("theta (degrees)"); plt::ylabel("Ex"); plt::legend();

	plt::subplot(2,2,2);
	plt::named_plot("Mie at far-field",thetaas,Ey_mie);
	plt::named_plot("FEM prop to far-field",thetaas,Ey_huygen);
	plt::xlabel("theta (degrees)"); plt::ylabel("Ey"); plt::legend();

	plt::subplot(2,2,3);
	plt::named_plot("Mie at far-field",thetaas,Ez_mie);
	plt::named_plot("FEM prop to far-field",thetaas,Ez_huygen);
	plt::xlabel("theta (degrees)"); plt::ylabel("Ez"); plt::legend();

	plt::subplot(2,2,4);
	plt::named_plot("Mie at far-field",thetaas,rcs_mie);
	plt::named_plot("FEM prop to far-field",thetaas,rcs_huygen);
	plt::xlabel("theta (degrees)"); plt::ylabel("10*log(RCS/lam^2)"); plt::legend();

	plt::show();
}
