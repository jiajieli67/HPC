#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;


MatrixXd product_multi(const MatrixXd &A, const MatrixXd &B, int n) {
	MatrixXd C = MatrixXd::Zero(n,n);
	int nthreads,p;
	#pragma omp parallel shared(nthreads,p)
	{
		nthreads = omp_get_num_threads();
		//nthreads = 4;
		p = n/nthreads;
		//printf("p = %d\n", p);
		//#pragma omp for
		for(int i=0; i<=n-p; i+=p) {
			//printf("i = %d\n", i);
		  	for(int j=0; j<=n-p; j+=p) {
				MatrixXd sum = MatrixXd::Zero(p,p);
				for(int x=0; x<nthreads; x++) {
					  sum = sum + A.block(i,x*p,p,p) * B.block(x*p,j,p,p);
				}
				C.block(i,j,p,p) << sum;
		 	}
			//cout<<"block = "<<BlockC<<endl;
			//printf("i = %d, j = %d\n", i,j);
	 	}

 	}
	return C;
}



int main()
{

  int n = 8;

  /// nxn Matrix filled with random numbers between (-1,1)
  MatrixXd A = MatrixXd::Random(n,n);
  MatrixXd B = MatrixXd::Random(n,n);
  MatrixXd C = MatrixXd::Zero(n,n);

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      double sum = 0;
      for(int k=0; k<n; k++) {
        sum = A(i,k) * B(k,j);
      }
      C(i,j) = sum;
    }
  }
	cout<<"C1 = "<<A*B<<endl;
	//cout<<"A = "<<A<<endl;
	//cout<<"B = "<<B<<endl;
	//cout<<"blockA"<<A.block(0,0,2,2)<<endl;
	cout<<"C2 = "<<product_multi(A,B,n)<<endl;
}
