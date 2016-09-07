#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include "solver.cpp"
using namespace std;

//int N;

//void solve_matrix_eq(gsl_vector * solution,
//                     gsl_matrix * system,
//                     gsl_vector * RightPart)
////Solve SLE Ax=b, where A = system, b = RightPart, x = solution
//{
//    int i;
//    gsl_permutation * p = gsl_permutation_alloc (N*N);
//    gsl_linalg_LU_decomp (system, p, &i);
//    gsl_linalg_LU_solve (system, p, RightPart, solution);
//}


double (*phi)(double, double, int);

int main(/*int argc, char *argv[]*/)
{
    //double intStep = 4.;
    //initGaussInt();
    rect_area area;
    area.x0 = 0.;
    area.x1 = 1.;
    area.y0 = 0.;
    area.y1 = 1.;

//    basis test(1,3,0.,1.,0.,1.);
    basis test(1,3,area);

//    phi = (&basis::value_temp);
    cout << "Hello World!" <<test.value_temp(0.,0.1,0) <<endl;
    return 0;
}
