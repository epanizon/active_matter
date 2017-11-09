#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>

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
// double lx{10};
// double ly{10};
double qq{2*M_PI/1}; 
double U0{0.1};
double dt{0.01};
double m{1.0};
int N{1};
int nstep{1000};

// READ PARAMETERS FUNCTION FROM ROBERTO
int read_parameters(int argc, char *argv[]){
  char parameter,*endp;
  for(int i=1;i<=argc-1;i++){
    parameter=*argv[i];
    switch(parameter){
    case 'n':nstep =strtod(argv[i]+1,&endp);break;
    case 'v':vel0  =strtod(argv[i]+1,&endp);break;
    case 'g':gam   =strtod(argv[i]+1,&endp);break;
    case 'm':m     =strtod(argv[i]+1,&endp);break;
    case 'T':kT    =strtod(argv[i]+1,&endp);break;
    case 'U':U0    =strtod(argv[i]+1,&endp);break;
    case 'o':omega0=strtod(argv[i]+1,&endp);break;
    case 'd':dt    =strtod(argv[i]+1,&endp);break;
    case 'q':qq    =strtod(argv[i]+1,&endp);break;
    }
  }
return 0;
}


// struct of pair of double for force calculations
struct ForcesXY{
  double fx;
  double fy;
};

// PARTICLE CLASS DEFINITION
class Particle{
  double *pos;         // position
  double *vel;         // force 
  double theta {0.0};  // angle of direction
  double omega;        // torque
  // TO BE IMPLEMENTED
  //  vector<Particle*> NN{} // neighbors
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

  // OUTPUT OF DATA. QUITE SHITTY FOR NOW.
  string print_pos(){ return to_string(pos[0]) + " " + to_string(pos[1]);} 
  string print_vel(){ return to_string(vel[0]) + " " + to_string(vel[1]);}
  string print_theta(){return to_string(theta);}

  // TIME EVOLUTION
  void step();
  void force();
  void force_random();
  void force_substrate();
  void force_interaction();
  ForcesXY fsub(double, double);

};

// function to make an evolution step
void Particle::step(){
  pos[0] += (vel[0]+cos(theta)*vel0)*dt/m/gam;
  pos[1] += (vel[1]+sin(theta)*vel0)*dt/m/gam;
  theta  += (omega + omega0)*dt/m/gam;
}

// function to calculate forces
void Particle::force(){
  vel[0] = 0.0;
  vel[1] = 0.0;
  omega = 0.0;
  force_random();
  force_substrate();
  force_interaction();
}

void Particle::force_random(){
  vel[0] += gaussrand()*sqrt(kT*4*gam*m/dt);
  vel[1] += gaussrand()*sqrt(kT*4*gam*m/dt);
  omega  += gaussrand()*sqrt(kT*2*gam*m/dt);
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

void Particle::force_interaction(){
  // to be created
  //  if (NN == nullptr) return;  
  return;
}



  // MAIN 
int main(int argc, char** argv){
  // one day maybe parameters form file
  if (argc){
    cout << "# INPUT FROM LINE.\n";
    read_parameters(argc, argv);
  }
  Particle P1;
  P1.put_pos(1.1,3.2);
  P1.put_theta(0.2);  
  // start data
  cout<<"#initial pos is \n";
  cout<<"0.0 "<<P1.print_pos()<<" "<<P1.print_vel()<<" "<<P1.print_theta()<<"\n";

  for (int i=1; i<=nstep;i++){
    // evolution 
    P1.force();
    P1.step();
    // print
    cout<<i*dt<<" "<<P1.print_pos()<<" "<<P1.print_vel()<<" "<<P1.print_theta()<<"\n";
  }

  return 0;
}
