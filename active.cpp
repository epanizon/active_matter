#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>

using namespace std;

class Particle{
  double pos[2];
  double theta;
public:
  void put_pos(double x, double y){
    pos[0] = x;
    pos[1] = y;
  }
  void put_theta(double dir){
    theta = dir;
  }
  double* get_pos(){
    return pos;
  }
  double get_theta(){
    return theta;
  }
};


  int main(int argc, char** argv){
    if (argc){cout << "argc is "<< argc <<"\n";}
    Particle P1;
    P1.put_pos(1.1,3.2);
    P1.put_theta(0.2);
    
    cout << "pos is " << P1.get_pos() << "\n";
    cout << "dir is " << P1.get_theta() << "\n";
    return 0;
  }
