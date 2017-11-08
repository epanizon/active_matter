#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>

using namespace std;

class Particle{
  double *pos;
  double *vel;
  double theta {0.0};
public:
  Particle():{
    pos = new double[2]{0.0,0.0};
    vel = new double[2]{0.0,0.0};
  }
  ~Particle():{
    delete[] pos;
    delete[] vel;
    pos = nullptr;
    vel = nullptr;
  }
  void put_pos(double x, double y){
    pos[0] = x;
    pos[1] = y;
  }
  void put_vel(double x, double y){
    vel[0] = x;
    vel[1] = y;
  }
  void put_theta(double dir){theta = dir;}
  // shitty, but for now...
  ostream& print_pos(){
    ostream& os;
    os << pos[0] << " " << pos[1];
  } 
  void print_vel(){cout << vel[0] << " " << vel[1] << "\n";}
  void print_theta(){cout << theta << "\n";}

  

};


int main(int argc, char** argv){
    if (argc){cout << "argc is "<< argc <<"\n";}
    Particle P1;
    P1.put_pos(1.1,3.2);
    P1.put_theta(0.2);
    
    // mah
    cout << "pos is " << print_pos() << "\n";
    return 0;
}
