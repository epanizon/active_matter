#include <iostream>
#include <omp.h>

using namespace std;
#define NUM_THREADS 1
int main(){  
  omp_set_num_threads(NUM_THREADS);
  # pragma omp parallel
  {
    int ID = omp_get_thread_num();
    cout << "hello(" << ID <<")";
    cout << " world(" << ID<<")\n";
  }

  static long num_steps = 100000000;
  double step;
  double t;
  double pi{};
  double sum{};
  step=1.0/(double)num_steps;
  t = omp_get_wtime();
  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel 
  {
    double x;
  #pragma omp for reduction (+:sum)
    for (int i = 0; i< num_steps; i++){
      x = (i+0.5)*step;
      sum = sum + 4.0/(1.0+x*x);
    }
  }
  pi = sum * step;
 t = omp_get_wtime()-t;
  cout << "pi="<< pi<< " t="<<t<<"\n"; 
}
