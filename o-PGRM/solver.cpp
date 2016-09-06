
class solver
{
private:
    gsl_matrix *sys;
    gsl_vector *rightpart, *solution;

    basis basis_of_system;

    double (*rightpart_f)  (double x, double y);        //допустимо
    double (*boundary_phi) (double x, double y);        //может быть
    double (*omega)        (double x, double y);        //пока что самостоятельная
    double (*omega2)       (double x, double y);        //аналогично
    double (*structure)    (double x, double y, int n); //по-сути enum-овый выбор

public:
    solver();
    ~solver();

    void form_matrix();
    void form_rightpart();
    void form_system();
    void solve();
    double value_at         (double x, double y);
};

