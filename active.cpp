#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>

using namespace std;

// gaussian random generator from Franco
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
double kT;
double gamma;
double lx;
double ly;

class Particle{
  double *pos;        // position
  double *vel;        // force 
  double theta {0.0}; // angle of direction
  double omega;       // torque
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
  string print_pos(){
    return to_string(pos[0]) + " " + to_string(pos[1]);
  } 
  string print_vel(){
    return to_string(vel[0]) + " " + to_string(vel[1]);
  }
  string print_theta(){return to_string(theta);}
  // TIME EVOLUTION
  void step(double);
  void force();
  void force_random();
  void force_substrate();
  void force_interaction();
};

void force(){
  vel[0] = 0.0;
  vel[1] = 0.0;
  omega = 0.0;
  force_random();
  force_substrate();
  force_interaction();
}

void Particle::step(double dt){
  pos[0] += (vel[0]+cos(theta)*vel0)*dt;
  pos[1] += (vel[1]+sin(theta)*vel0)*dt;
  theta  += (omega + omega0)*dt;
}

void Particle::force_random(){
  vel[0] += gaussrand()*kT*4*gamma;
  vel[1] += gaussrand()*kT*4*gamma;
  omega  += gaussrand()*kT*2*gamma;
}

  // MAIN 
int main(int argc, char** argv){
    if (argc){cout << "argc is "<< argc <<"\n";}
    Particle P1;
    P1.put_pos(1.1,3.2);
    P1.put_theta(0.2);  
    // check
    cout << "pos is " << P1.print_pos() << ' ' << P1.print_vel() << ' ' << P1.print_theta() << "\n";
    return 0;
}
