template <class T>
T sq(const T &a){
	return a*a;
}

template <class T>
T Max(T a, T b){
	return (a>b?a:b);
}

template <class T>
T Min(T a, T b){
	return (a<b?a:b);
}

//Node class
class node{
	public:
	double nx,ny,nz; //coordinates
	char ndtype; //s - surface, i - interior
	vector<int> surfno;
	node(double,double,double);
	node(){};
};

node::node(double x, double y, double z){
	nx=x;ny=y;nz=z;}

//Edge class
class edge{
public:
    int ni,nj; //node numbers of edge
	//double edlength;//edge length
    //char edbtype; //s - surface, i - interior
	vector <int> els; //elements that share this edge
	edge(int,int);
	edge(){};
	//void setlength(const vector<node> *);
    bool operator== (const edge &);
};

edge::edge(int n1,int n2){
	ni=n1; nj=n2;
	}


// void edge::setlength(const vector<node> * nlist){
// 	double nix = (*nlist)[ni].nx, niy = (*nlist)[ni].ny, niz = (*nlist)[ni].nz,
// 	       njx = (*nlist)[nj].nx, njy = (*nlist)[nj].ny, njz = (*nlist)[nj].nz;
// 	edlength=sqrt(sq(nix-njx)+sq(niy-njy)+sq(niz-njz));
// }

bool edge::operator== (const edge &e1){
	if(((ni==e1.ni)&&(nj==e1.nj))||((nj==e1.ni)&&(ni==e1.nj))) //safe because comparing integers
	{return true;}
	else
	{return false;}
}


class element{
public:

	vector<int> nodes;//node numbers (in nodelist)
	vector<int> edges;//edge numbers (in edgelist)
	char eltype;//Type: a - air, s - scatterer
	cmp eleps;//element epsilon
	cmp elmu;//element mu
	double V; //element volume
	vector<double> l; //length of edges of element
	int num[6][2]; //local edge to node pair
	vector<double> a,b,c,d;
	// vector<Matrix3d> A
	//
	// Vector3d T(int i, Vector3d rbar) // Method/function defined inside the class
	// {
	//
	// 	return (something*rbar + something)
	// }

	//Functions
	element(){};
	element(int,int,int,int);

	};

element::element(int n1,int n2,int n3, int n4){
    nodes.resize(4); edges.reserve(6);
    nodes[0]=n1;nodes[1]=n2;nodes[2]=n3;nodes[3]=n4;
}


class surf_element{
public:
	vector<int> glob_el_num;// element number (in elementlist)
	vector<int> face_nodes; // nodes on the exposed face of the element
	vector<int> in_node; 		// interior node of the element

	//double norx, nory, norz;
	double x,y,z,area; // coordinates of centroid and area of triangle
	double r, theta, phi; // polar coordinates of centroid
	cmp Ex, Ey, Ez; 	 // electric field at centroid
	cmp vx, vy, vz; 	 // curl of electric field at centroid

	double nlx,nly,nlz;

	surf_element(){};
	surf_element(int);
	};

surf_element::surf_element(int n1)
{
	face_nodes.reserve(3);
	glob_el_num.push_back(n1);
}
