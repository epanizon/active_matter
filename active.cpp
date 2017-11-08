#include "Random.h"
#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

class Particle{
  vector<double> pos[2];
  double theta;
public:
  void put_pos(double x, double y){
    pos[1] = x;
    pos[2] = y;
  }
  void put_dir(double dir){
    theta = dir;
  }
  vector<double> get_pos(){
    return pos;
  }
  vector<double> get_dir(){
    return dir;
  }
}


  int main(int argc, char** argv){
    if (argc){cout << "argc is "<< argc <<"\n";}
    Particle P1;
    P1.put_pos(1.1,3.2);
    P1.put_dir(0.2);
    
    cout << "pos is " << P1.get_pos() << "\n";
    cout << "dir is " << P1.get_dir() << "\n";
    return 0;
  }
