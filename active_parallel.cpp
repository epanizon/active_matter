#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <fstream>
#include <unordered_set>
#include <omp.h>

using namespace std;

// GAUSSIAN RANDOM GENERATOR FROM FRANCO
double gaussrand(){
  static double V1, V2, S;
  static int phase = 0;
  double XX;
    if(phase == 0) {
	do {
		double U1 = (double)rand() / RAND_MAX;
		double U2 = (double)rand() / RAND_MAX;
		V1 = 2 * U1 - 1;
		V2 = 2 * U2 - 1;
		S = V1 * V1 + V2 * V2;
	} while(S >= 1 || S == 0);
	XX = V1 * sqrt(-2 * log(S) / S);
    } else
	XX = V2 * sqrt(-2 * log(S) / S);
  phase = 1 - phase;
  return XX;
};

// GLOBAL VARIABLES
double vel0{1.0};
double omega0{0.0};
double kT{0.5};
double gam{2};
double lx{20};
double ly{12*sqrt(3)};
double qq{2*M_PI/1}; 
double U0{0.1};
double dt{0.01};
double m{1.0};
double V0{1};
double R0{0.7};
double Rlist[2]{R0,R0};
double R2{};
double Rflok{3};
double Rflok2{};
int N{1};
int Nmax{};
int Nlist{};
int Nlistx{};
int nstep{1000};
int seed_number{0};
int nprint{10};
ofstream fout ("conf.xyz", ofstream::out);


// READ PARAMETERS FUNCTION FROM ROBERTO
int read_parameters(int argc, char *argv[]){
  char parameter,*endp;
  //initialize random number generator
  srand (time(NULL));
  for(int i=1;i<=argc-1;i++){
    parameter=*argv[i];
    switch(parameter){
    case 'N':N     =strtod(argv[i]+1,&endp);break;
    case 'n':nstep =strtod(argv[i]+1,&endp);break;
    case 'v':vel0  =strtod(argv[i]+1,&endp);break;
    case 'g':gam   =strtod(argv[i]+1,&endp);break;
    case 'm':m     =strtod(argv[i]+1,&endp);break;
    case 'T':kT    =strtod(argv[i]+1,&endp);break;
    case 'U':U0    =strtod(argv[i]+1,&endp);break;
    case 'V':V0    =strtod(argv[i]+1,&endp);break;
    case 'o':omega0=strtod(argv[i]+1,&endp);break;
    case 'd':dt    =strtod(argv[i]+1,&endp);break;
    case 'q':qq    =strtod(argv[i]+1,&endp);break;
    case 's':seed_number = strtod(argv[i]+1,&endp);break;
    case 'x':lx    =strtod(argv[i]+1,&endp);break;
    case 'y':ly    =strtod(argv[i]+1,&endp);break;
    case 'r':R0    =strtod(argv[i]+1,&endp);break;
    case 'f':Rflok =strtod(argv[i]+1,&endp);break;
    case 'p':nprint=strtod(argv[i]+1,&endp);break;
    }
  }
  R2=R0*R0;
  Rflok2=Rflok*Rflok;
  if (seed_number){srand(seed_number);}
  // PBC corrected for underlying substrate
  lx = ((int)(lx*qq/(2*M_PI)))*2*M_PI/qq;
  ly = ((int)(ly*qq/(2*M_PI*sqrt(3))))*2*M_PI/qq*sqrt(3.);
  Nmax = floor(lx/R0)*floor(ly/R0);
  if (Nmax < N) {
    cout << "TOO MANY PARTICLES IN TOO SMALL A BOX";
    return 1;
  }

  // LIST CUTOFF LARGEST BETWEEN VICSEK AND LJ
  Rlist[0]=max(R0,Rflok);
  // CALCULATION NUMBERS FOR CELL LISTS
  Nlistx = floor(lx/max(R0,Rflok));
  Rlist[0] = lx/Nlistx;
  Nlist    = Nlistx * floor(ly/max(R0,Rflok));  
  Rlist[1] = ly/floor(ly/max(R0,Rflok));
 
  return 0;
}

// PRINT PARAMETER FOR CHECK
void print_parameters(){
  cout << "#N  " << N     << endl;
  cout << "#n  " << nstep << endl; 
  cout << "#v  " << vel0  << endl; 
  cout << "#g  " << gam   << endl; 
  cout << "#m  " << m     << endl; 
  cout << "#T  " << kT    << endl; 
  cout << "#U  " << U0    << endl; 
  cout << "#o  " << omega0<< endl; 
  cout << "#dt " << dt    << endl; 
  cout << "#q  " << qq    << endl;
  cout << "#s  " << seed_number << endl;
  cout << "#lx " << lx    << endl; 
  cout << "#ly " << ly    << endl; 
  cout << "#R0 " << R0    << endl;
  cout << "#Rf " << Rflok << endl; 
  cout << "#p  " << nprint<< endl;
  cout << "#Reynold's Number (should be << 1) = " << m*vel0/R0/gam/gam << endl;
  
}

// struct of pair of double for force calculations
struct ForcesXY{
  double fx;
  double fy;

  // OVERLOAD OPERATOR- FOR WRITING -F
ForcesXY operator-() const {
  ForcesXY Y;
  Y.fx = -fx;
  Y.fy = -fy;
  return Y;
}
};

// random init on a (sparce) square lattice of size R0.
void initialize_positions(double* X, double* Y, double* Theta){
  vector<int> ran;
  int trial;
  int nx, ny;
  int nxmax;

  nxmax = floor(lx/R0);
  for (int i=0;i<Nmax;i++){ran.push_back(i);}  

  for (int i=0; i<N; i++){
    trial = rand()%(Nmax-i-1);
    nx = ran[trial]%nxmax;
    ny = ran[trial]/nxmax;
    X[i] = (nx + 0.5)*R0;
    Y[i] = (ny + 0.5)*R0;
    Theta[i] = (double)rand()/RAND_MAX*M_PI;
    ran.erase(ran.begin()+trial);
  } 
}

// functions for keeping particle in BOX.
void pbc(double* X, double* Y){
  for (int i=0; i<N; i++){
    X[i] -= floor(X[i]/lx)*lx;
    Y[i] -= floor(Y[i]/ly)*ly;
  }
}

// PBC TO GET REAL VECTOR BETWEEN TWO MOTHERFUCKING PARTICLES.
void pbc(double P1x, double P2x, double P1y, double P2y, double& x, double& y){
  x = P2x-P1x;
  x-=floor(x/lx+0.5)*lx;
  y = P2y-P1y;
  y-=floor(y/ly+0.5)*ly;
}

// two body interaction (repulsive LJ) 
ForcesXY fint(double x, double y){
  ForcesXY F;
  double rr2;
  double r0surr2;
  double r0surr6;
  double ff;
  rr2 = {x*x + y*y};
  r0surr2 = R2/rr2;
  r0surr6 = r0surr2*r0surr2*r0surr2;
  ff = V0*12./rr2*(r0surr6*r0surr6 - r0surr6);
  F.fx = ff*x;
  F.fy = ff*y;
  return F;
}

// substrate interaction and random kicks
ForcesXY fsub(double x, double y){
  ForcesXY sub;
  sub.fx  = U0*qq*2*sin(qq*x)*cos(qq*y/sqrt(3));
  sub.fy  = U0*qq*( 2/sqrt(3.)*cos(qq*x)*sin(qq*y/sqrt(3))+sin(2*y*qq/sqrt(3))/sqrt(3.)*2);
  return sub;
}

// PRINT DATA
void print_data(int it, int N, double* X, double* Y, double* Vx, double* Vy, double* Theta){
  cout << it*dt << " " ;
  for (int i=0; i < N; i++){
    cout <<X[i]<<" "<<Y[i]<< " "<<Vx[i]<< " "<<Vy[i]<<" "<<cos(Theta[i])<<" "<<sin(Theta[i])<< " ";
  }
  cout << endl;
}

// PRINT CONFIG
void print_config(int it, int N, double* X, double* Y, double* Vx, double* Vy, double* Theta){
  fout << N << endl<<"# t=" << it*dt << endl;
  for (int i=0; i < N; i++){
    fout << "P"<< i <<" "<<X[i]<<" "<<Y[i]<<" 0.0 "<<cos(Theta[i])*R0<<" "<<sin(Theta[i])*R0<<" 0.0 "<<Vx[i]<<" "<<Vy[i]<< " 0.0 "<<endl;
  }
  // fout << endl << endl;
}

// LIST ASSIGN
void list_update(double* X, double* Y, vector<int>& head, vector<int>& lscl){
  int cellnum;
  //performs list update
  // -1 MEANS EMPTY.
  for (int i=0; i<Nlist;i++) head[i] = -1;
  for (int i=0; i<N; i++){
    cellnum = floor(X[i]/Rlist[0])+floor(Y[i]/Rlist[1])*Nlistx;
    lscl[i] = head[cellnum];
    head[cellnum] = i;
  }
}

// CORRECT MODULO OF NEIGHBOUR CELL
int neighcell(int c1, int ix, int iy, int Nx, int Ntot){
 int newix;
 int newiy;
 newix = (c1+ix)%Nx;
 if (newix<0) newix+=Nx;
 newiy = (c1/Nx+iy)%(Ntot/Nx);
 if (newiy<0) newiy+=Ntot/Nx;
 return newix + newiy*Nx;
}


// MAIN 
int main(int argc, char** argv){
  // one day maybe parameters form file
  if (argc){
    int error{};
    // cout << "# INPUT FROM LINE.\n";
    error = read_parameters(argc, argv);
    print_parameters();

    if (error){
      cerr << "TOO MANY PARTICLES. RANDOM INIT FAILS " << N << " over " << Nmax;
      return 0;
    }

  }

  // C VERSION. 
  ForcesXY F{0};
  double* Xpos = new double[N]{};
  double* Ypos = new double[N]{};
  double* Xvel = new double[N]{};
  double* Yvel = new double[N]{};
  double* Theta = new double[N]{};
  double* Omega = new double[N]{};
  double* avgtheta = new double[N]{};
  int* thetacount = new int[N]{};
  vector<int> head;
  vector<int> lscl;  

  //debug move 
  double maxmove{};
  double move{};
  int   imaxmove{};
  int  imaxmove2;


  // initialization of the N particles
  initialize_positions(Xpos, Ypos, Theta);
  print_config(0, N, Xpos, Ypos, Xvel, Yvel, Theta);

  head.resize(Nlist);
  lscl.resize(N);

  list_update(Xpos, Ypos, head, lscl);

  for (int it=1; it<=nstep;it++){
    // SUBSTRATE
    for (int i=0; i<N; i++){
      F = fsub(Xpos[i],Ypos[i]);
      Xvel[i]  = F.fx;
      Yvel[i]  = F.fy;   
    }


    for (int i=0; i<N; i++){
      double ranx;
      double rany;
      ranx = gaussrand()*R0*sqrt(3*kT*2/gam/m/dt);
      rany = gaussrand()*R0*sqrt(3*kT*2/gam/m/dt);
      Xvel[i]  += ranx;
      Yvel[i]  += rany;
      Omega[i] = gaussrand()*R0*sqrt(3*kT*2/gam/m/dt);    
    }

    // INIT INTERACTION ANGLES
    for (int i=0; i<N; i++){
      avgtheta[i] = 0.0;
      thetacount[i] = 0;
    }

    for (int cell1=0; cell1<Nlist; cell1++){
      // CYCLE ON ALL ADIACENT CELLS
      for (int j=-1;j<2;j++){
	for (int k=-1;k<2;k++){
	  // SCALAR NUMBER OF CELL
	  int cell2 = neighcell(cell1, j, k, Nlistx, Nlist); 
	  // uno WILL CYCLE ON ATOMS IN FIRST CELL
	  int uno = head[cell1];
	  // due WILL CYCLE ON ATOMS IN SECOND CELL
	  int due;
	  double x, y, rr;
	  while (uno > -1){
	    // DEBUG
	    due = head[cell2];
	    while (due > -1){
	      if (uno < due){
		pbc(Xpos[uno], Xpos[due], Ypos[uno], Ypos[due], x, y);
		rr = x*x+y*y;
		if (rr < Rflok2){
		  avgtheta[uno]  +=Theta[due];
		  thetacount[uno]+=1;
		  avgtheta[due]  +=Theta[uno];
		  thetacount[due]+=1;
		}
		
		if (rr < R2){
		  F = fint(x, y);
		  Xvel[due] += F.fx;
		  Yvel[due] += F.fy;
		  Xvel[uno] -= F.fx;
		  Yvel[uno] -= F.fy;
		}
	      }
	      due = lscl[due];
	    }
	    uno = lscl[uno];
	  }
	}	
      }
    }

    for (int i=0; i<N;i++){
      if (thetacount[i]) Omega[i] += (avgtheta[i]/thetacount[i]);
      Xpos[i]  += (Xvel[i]+cos(Theta[i])*vel0)*dt/m/gam;
      Ypos[i]  += (Yvel[i]+sin(Theta[i])*vel0)*dt/m/gam;
      Theta[i] += (Omega[i] + omega0)*dt/m/gam;
    }


    // PBC enforce.
    pbc(Xpos, Ypos);  
    list_update(Xpos, Ypos, head, lscl);

    if (it%nprint == 0){print_data(it, N, Xpos, Ypos, Xvel, Yvel, Theta);}
    if (it%nprint == 0){print_config(it, N, Xpos, Ypos, Xvel, Yvel, Theta);} 
  }


  delete[] Xpos, Ypos, Xvel, Yvel, Theta, Omega;
  delete[] thetacount;
  delete[] avgtheta;
  fout.close();
  Xpos = Ypos = Xvel = Yvel = Theta = Omega = nullptr;
  thetacount = nullptr;
  avgtheta = nullptr;
  return 0;
}
