#include <iostream>
#include <stdio.h>
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



int main(/*int argc, char *argv[]*/)
{
    rect_area area;
    area.x0 = -M_PI;
    area.x1 =  M_PI;
    area.y0 = -M_PI;
    area.y1 =  M_PI;
    solver task1(1,8,area,1);
    task1.form_system();
    //task1.solve();
    task1.print_solution();
//    task1.plot();

    cout<<"\n\n"<<task1.gauss_integral_left(2,{0,0,0,0})<<"\n\n";

    //cout << "Hello World!" <<test.value_temp(0.,0.1,0) <<endl;
    return 0;
}
