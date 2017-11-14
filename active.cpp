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
double R2{};
double Rflok{3};
double Rflok2{};
int N{1};
int Nmax{};
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
  // particles too dense -> problem with random init
  if( (lx > 0) && (ly > 0) ){
    Nmax = floor(lx/R0)*floor(ly/R0); 
  } else { Nmax = N+2;} 
  if (N>=Nmax-1) {return 1;}

return 0;
}
// PRINT PARAMETER FOR CHECK
void print_parameters(){
  cout << "#N  " << N     << "\n";
  cout << "#n  " << nstep << "\n"; 
  cout << "#v  " << vel0  << '\n'; 
  cout << "#g  " << gam   << '\n'; 
  cout << "#m  " << m     << '\n'; 
  cout << "#T  " << kT    << '\n'; 
  cout << "#U  " << U0    << '\n'; 
  cout << "#o  " << omega0<< '\n'; 
  cout << "#dt " << dt    << '\n'; 
  cout << "#q  " << qq    << '\n';
  cout << "#s  " << seed_number << '\n';
  cout << "#lx " << lx    << '\n'; 
  cout << "#ly " << ly    << '\n'; 
  cout << "#R0 " << R0    << '\n';
  cout << "Rflo" << Rflok << '\n'; 
  cout << "#p  " << nprint<< '\n';
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

// PARTICLE CLASS DEFINITION
class Particle{
  double *pos;         // position
  double *vel;         // force 
  double theta;  // angle of direction
  double omega;        // torque
public:
  // CONSTRUCTOR
  Particle(){
    pos = new double[2]{0.0,0.0};
    vel = new double[2]{0.0,0.0};
  }
  // DESTRUCTOR
  ~Particle(){
    delete[] pos;
    delete[] vel;
    pos = nullptr;
    vel = nullptr;
  }

  // BASIC INITIALIZATION FUNCTIONS
  void put_pos(double x, double y){
    pos[0] = x;
    pos[1] = y;
  }
  void put_vel(double x, double y){
    vel[0] = x;
    vel[1] = y;
  }
  void put_theta(double dir){theta = dir;}
  void put_force(ForcesXY F){vel[0]+=F.fx; vel[1]+=F.fy;}
  void add_omega(double f){omega += f;}
 //to be added when neighbour list 
 // void add_neigh(int neigh){NN.push_back(neigh);}
  void random_init(int, int);

  // OUTPUT OF DATA. QUITE SHITTY FOR NOW.
  string print_pos(){ return to_string(pos[0]) + " " + to_string(pos[1]);} 
  double give_posx(){return pos[0];}
  double give_posy(){return pos[1];}
  double give_theta(){return theta;}
  string print_vel(){return to_string(vel[0]) + " " + to_string(vel[1]);}
  string print_theta(double scale){return to_string(cos(theta)*scale) + " " + to_string(sin(theta)*scale);}

  // TIME EVOLUTION
  void step();
  void force();
  void force_random();
  void force_substrate();
  void force_interaction();
  ForcesXY fsub(double, double);
  // PBC
  void pbc();

}; // END OF CLASS SCOPE

// random init

void initialize_positions(Particle* P){
  vector<int> ran;
  int trial;
  int nx, ny;
  int nxmax;

  if (lx ==0 && ly ==0) {
    nxmax = int (sqrt(N*1.0));
  } else if (lx ==0 && ly >0) {
    nxmax = Nmax / floor (ly/R0) + 1; 
  } else {nxmax = floor(lx/R0);}

  for (int i=0;i<Nmax;i++){ran.push_back(i);}

  for (int i=0; i<N; i++){
    trial = rand()%(Nmax-i-1);
    nx = ran[trial]%nxmax;
    ny = ran[trial]/nxmax;
    P[i].random_init(nx,ny);
    ran.erase(ran.begin()+trial);
  }
}

void Particle::random_init(int nx, int ny){
  pos[0] = ( nx + 0.5) * R0;
  pos[1] = ( ny + 0.5) * R0;
  theta  = (double)rand() / RAND_MAX * M_PI;
}

// functions for keeping particle in BOX.
void Particle::pbc(){
  if (lx) {pos[0] -= floor(pos[0]/lx)*lx;}
  if (ly) {pos[1] -= floor(pos[1]/ly)*ly;}
}

// PBC TO GET REAL VECTOR BETWEEN TWO MOTHERFUCKING PARTICLES.
void pbc(Particle* P1, Particle* P2, double& x, double& y){
  x = (*P2).give_posx()-(*P1).give_posx();
  if (lx){x-=floor(x/lx+0.5)*lx;}
  y = (*P2).give_posy()-(*P1).give_posy();
  if (ly){y-=floor(y/ly+0.5)*ly;}
}

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

// function to make an evolution step
void Particle::step(){
  pos[0] += (vel[0]+cos(theta)*vel0)*dt/m/gam;
  pos[1] += (vel[1]+sin(theta)*vel0)*dt/m/gam;
  // PBC refolding
  if (lx || ly){
    pbc();
  }
  theta  += (omega + omega0)*dt/m/gam;
}

// function to calculate NON-INTERACTING forces
void Particle::force(){
  vel[0] = 0.0;
  vel[1] = 0.0;
  omega = 0.0; 
  force_random(); 
  force_substrate();
}

void Particle::force_random(){
  double fx, fy;
  fx =  gaussrand()*sqrt(kT*2/gam/m/dt);
  vel[0] += fx;
  fy = gaussrand()*sqrt(kT*2/gam/m/dt);
  vel[1] += fy;
  omega  += gaussrand()*sqrt(kT*2/gam/m/dt);
}

void Particle::force_substrate(){
  ForcesXY forces;
  forces = fsub(pos[0],pos[1]);
  vel[0] += forces.fx;
  vel[1] += forces.fy;
}

ForcesXY Particle::fsub(double x, double y){
  ForcesXY sub;
  sub.fx = U0*2*qq*sin(qq*x)*cos(qq*y/sqrt(3));
  sub.fy = U0*2*qq*cos(qq*x)*sin(qq*y/sqrt(3))+U0*sin(2*y/qq);
  return sub;
}

void print_data(int it, int N, Particle* P){
  cout << it*dt << " " ;
  for (int i=0; i < N; i++){
    cout << P[i].print_pos()<<" "<<P[i].print_theta(R0/2.)<<" ";
  }
  cout << "\n";
}

void print_config(int it, int N, Particle* P){
  fout << N << "\n# t=" << it*dt << "\nempty line ";
  for (int i=0; i < N; i++){
    fout << "P"<< i << " " << P[i].print_pos() << " "<< P[i].print_theta(R0/2.) << " " << P[i].print_vel() <<"\n";
  }
  fout << "\n\n";
}

// LIST ASSIGN
void list_update(const Particle* P, vector<int>& clist){
  //does nothing
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

  // SLIGHTLY LESS STUPID VERSION. 
  Particle* P = new Particle[N];
  ForcesXY F{0};
  double* avgtheta = new double[N];
  int* thetacount = new int[N];
  double x;
  double y;
  double rr;

  // LIST PARAMETERS
  vector<int>** clist;
  clist = new vector<int>*(N);
  for(int i=0; i<N;i++) clist[i] = new vector<int>;

  // initialization of the N particles
  initialize_positions(P);
  print_config(0, N, P);

  // first list call
  list_update(P, clist);

  for (int it=1; it<=nstep;it++){
    // SUBSTRATE AND RANDOM
    for (int i=0; i<N; i++){
      P[i].force();
      avgtheta[i] = 0.0;
      thetacount[i] = 0;
    }
    // INTERACTION BETWEEN PARTICLES
    for (int i=0; i<N; i++){
      for (int j=i+1; j<N; j++){
        //LET's TRY WITHOUT INTERACTIONS.
	pbc(&P[i],&P[j],x,y);
	rr = x*x+y*y;
	if (rr < Rflok2){

	  	  avgtheta[i]+=P[j].give_theta();
	  	  thetacount[i]+=1;
	  	  avgtheta[j]+=P[i].give_theta();
	  	  thetacount[j]+=1;

	  if (rr < R2){
	    F = fint(x, y);
	    P[j].put_force(F);
	    P[i].put_force(-F);
	  }
	}	
      }
    }
    
    for (int i=0; i<N;i++){
      if (thetacount[i]) P[i].add_omega(avgtheta[i]/thetacount[i]);
      P[i].step();
    }


    //TO ADD LATER
    if (it%nprint == 0){print_data(it, N, P);}
    if (it%nprint == 0){print_config(it, N, P);}
    
  }
  
  delete[] P;
  delete[] thetacount;
  delete[] avgtheta;
  fout.close();
  P = nullptr;
  thetacount = nullptr;
  avgtheta = nullptr;
  return 0;
}
