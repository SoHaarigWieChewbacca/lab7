#include <iostream>
#include <cmath>

using namespace std;

void function(double y1, double y2, double* f, const double eta);


int main() {
    const double eta = 0.5;
    double y[2];
    y[0] = 1E-5;
    y[1] = sqrt(eta)*y[0];
    
    const double dx = 0.01;
    const int x_end = 100;
    const int N = x_end/dx + 1;
    
    double k1[2];
    double k2[2];
    double k3[2];
    
    for(int i = 0; i < N; i++) {
	function(y[0], y[1], k1, eta);
	function(y[0] + 0.5*dx*k1[0], y[1] + 0.5*dx*k1[1], k2, eta);
	function(y[0] - dx*k1[0] + 2*dx*k2[0], y[1] - dx*k1[1] + 2*dx*k2[1], k3, eta);
		
	for(int j = 0; j < 2; j++) {
	    y[j] += dx/6.0 * (k1[j] + 4*k2[j] + k3[j]);
	}
	
	cout << i*dx << "\t" << y[0] << "\t" << y[1] << endl;
    }
    
    
    return 0;
}


void function(double y1, double y2, double* f, const double eta) {
    double temp = y1;
    f[0] = y2;
    f[1] = (eta - abs(temp)*abs(temp))*temp;
}