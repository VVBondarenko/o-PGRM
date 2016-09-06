#include <iostream>
#include <gsl/gsl_linalg.h>

using namespace std;

int N;

void solve_matrix_eq(gsl_vector * solution,
                     gsl_matrix * system,
                     gsl_vector * RightPart)
//Solve SLE Ax=b, where A = system, b = RightPart, x = solution
{
    int i;
    gsl_permutation * p = gsl_permutation_alloc (N*N);
    gsl_linalg_LU_decomp (system, p, &i);
    gsl_linalg_LU_solve (system, p, RightPart, solution);
}

typedef struct rect_area {
    double x0, x1;
    double y0, y1;
} rect_area;


class basis{
public:
    basis();
    ~basis();
    int N;
    rect_area area;
    double value(double x, double y);
private:
    double simplePoly   (double x, double y);
    double chebyshevT   (double x, double y);
    double chebyshevU   (double x, double y);
    double Bsplines     (double x, double y);
    double fup3_poly    (double x, double y);
};

class solver
{
private:
    gsl_matrix *sys;
    gsl_vector *rightpart, *solution;

    double rightpart_f  (double x, double y);
    double boundary_f   (double x, double y);
    double omega        (double x, double y);
    double omega2       (double x, double y);
    double structure    (double x, double y, int n);

public:
    solver();
    ~solver();

    void form_matrix();
    void form_rightpart();
    void form_system();
    void solve();
    double value_at     (double x, double y);
};
int main(/*int argc, char *argv[]*/)
{
    cout << "Hello World!" << endl;
    return 0;
}
