#include <iostream>
#include <cmath>
#include <array>
#include <time.h>
#include <omp.h>
#include <fstream>

using namespace std;

struct Point{
    Point(double i_x = 0.0, double i_y = 0.0): x(i_x), y(i_y)
    {
    }

    double x;
    double y;
};

Point A(-2, 0), B(-1, 1), C(1, 1), D(2, 0), E(0, -2);
Point P_A(-3, -3), P_B(-3, 2), P_C(3, 2), P_D(3, -3);

//sides and inverse sides
double AB(double x)
{
    return x + 2;
}

double BC(double x)
{
    return 1; 
}

double CD(double x)
{
    return -x + 2; 
}

double ED(double x)
{
    return x - 2; 
}

double AE(double x)
{
    return -x - 2; 
}

double getX_CD(double y){
    return -y + 2;
}

double getX_ED(double y){
    return y + 2;
}

double getX_AB(double y){
    return y - 2;
}

double getX_AE(double y){
    return -y - 2;
}


// const int M = 160, N = 160; 
const int M = 80, N = 80;
const double h1 = (P_D.x - P_A.x) / M, h2 = (P_B.y - P_A.y) / N;
const double epsilon = h1 * h2, delta = 0.000001;

struct GridFunction
{
    GridFunction()
    {
	for(int i = 0; i <= M; ++i)
		for(int j = 0; j <= N; ++j)
			values[i][j] = 0.0;
    }

    double get(int i, int j) const { return values[i][j]; }

    double w_x(int i, int j) const { return (get(i + 1, j) - get(i, j)) / h1; }
    double w_y(int i, int j) const { return (get(i, j + 1) - get(i, j)) / h2; }

    void print(ostream& out) const 
    {
        for(const auto& value: values)
        {
            for(const auto& w: value)
                out << w << " ";
            out << endl;                
        }
    }

    array<array<double, N + 1>, M + 1> values;


};

GridFunction operator-(const GridFunction& i_l, const GridFunction& i_r)
{
    GridFunction result;
    int i;
    #pragma omp parallel default(shared) private(i)
    {
        int id = omp_get_thread_num();
        int numt = omp_get_num_threads();
        for(i = id + 1; i <= M - 1; i = i + numt)
            for(int j = 1; j <= N - 1; ++j)
                result.values[i][j] = i_l.values[i][j] - i_r.values[i][j];
    }
    return result;
}

GridFunction operator*(double i_l, const GridFunction& i_r)
{
    GridFunction result;
    int i;
    #pragma omp parallel default(shared) private(i) 
    {
        int id = omp_get_thread_num();
        int numt = omp_get_num_threads();
        for(i = id + 1; i <= M - 1; i = i + numt)
            for(int j = 1; j <= N - 1; ++j)
                result.values[i][j] = i_l * i_r.values[i][j];
    }
    return result;
}

double scalar(const GridFunction& u, const GridFunction& v)
{
    double sum = 0.0;
    int i;
    #pragma omp parallel default(shared) private(i) reduction(+: sum)
    {
        int id = omp_get_thread_num();
        int numt = omp_get_num_threads();
        for(i = id + 1; i <= M - 1; i = i + numt)
            for(int j = 1; j <= N - 1; ++j)
                sum += u.get(i, j) * v.get(i, j);
    
    }
    return h1 * h2 * sum;
}

double E_norm(const GridFunction& u){ return std::sqrt(scalar(u, u)); }



bool isBelongsToD(const Point& p)
{
    bool result;
    if(p.x <= A.x || p.x >= D.x || p.y >= B.y || p.y <= E.y)
        return false;

    if(p.y > A.y){
        if(p.x >= B.x && p.x <= C.x)
            return true;
        if(p.x < B.x)
            result = AB(p.x) > p.y;
        if(p.x > C.x)
            result = CD(p.x) > p.y;
    }else{
        if(p.x < E.x)
            result = AE(p.x) < p.y;
        if(p.x >= E.x)
            result = ED(p.x) < p.y;
    }
    return result;
}

double x_(double i) { return P_A.x + i * h1; }
double y_(double i) { return P_A.y + i * h2; }
double a(int i, int j)
{   

    double l = 0;
    Point P(x_(i - 0.5), y_(j - 0.5)), P_next(x_(i - 0.5), y_(j + 0.5));

    if(isBelongsToD(P) && isBelongsToD(P_next))
        return 1.;

    if(!isBelongsToD(P) && !isBelongsToD(P_next))
        return 1 / epsilon;
    
    //isBelongsToD(P) && !isBelongsToD(P_next)
    if(isBelongsToD(P))
    {
        if(P.x >= B.x && P.x <= C.x){
            l = B.y - P.y;
        }
        else if(P.x < B.x){
            l = AB(P.x) - P.y;
        }
        else if(P.x > C.x){
            l = CD(P.x) - P.y;
        }
        
    }

    //!isBelongsToD(P) && isBelongsToD(P_next)
    if(isBelongsToD(P_next)){
        if(P_next.x < E.x){
            l = P_next.y - AE(P_next.x);
        } else{
            l = P_next.y - ED(P_next.x);
        }
    }

    return l / h2 + (1 - l / h2) / epsilon;
}

double b(int i, int j)
{   
    double l = 0;
    Point P(x_(i - 0.5), y_(j - 0.5)), P_next(x_(i + 0.5), y_(j - 0.5));
    
    if(isBelongsToD(P) && isBelongsToD(P_next))
        return 1.;
    if(!isBelongsToD(P) && !isBelongsToD(P_next))
        return 1 / epsilon;

    //isBelongsToD(P) && !isBelongsToD(P_next)
    if(isBelongsToD(P)){
        if (P.y > A.y){
            l = getX_CD(P.y) - P.x;
        }else{
            l = getX_ED(P.y) - P.x;
        }
    }

    //!isBelongsToD(P) && isBelongsToD(P_next)
    if(isBelongsToD(P_next)){
        if (P_next.y > A.y){
            l = P_next.x - getX_AB(P_next.y);
        }else{
            l = P_next.x - getX_AE(P_next.y);
        }
    }
    return l / h1 + (1 - l / h1) / epsilon;
}

double S(int i, int j)
{   

    Point P1(x_(i - 0.5), y_(j - 0.5)), P2(x_(i + 0.5), y_(j - 0.5)), P3(x_(i + 0.5), y_(j + 0.5)), P4(x_(i - 0.5), y_(j + 0.5));
    int numberOfPointsInD = isBelongsToD(P1) + isBelongsToD(P2) + 
        isBelongsToD(P3) + isBelongsToD(P4);
   
    if(numberOfPointsInD == 4)
        return 1;
    if(numberOfPointsInD == 0)
        return 0;

    //partitioning within a semi-integer node area
    double parts_count = 5;

    double length_X = P2.x - P1.x;
    double length_Y = P3.y - P1.y;
    double h1_small = length_X / parts_count;
    double h2_small = length_Y / parts_count;
    double n = 0;
    double N = parts_count * parts_count;

    //counts of inside figure
    for(int i = 0; i < parts_count; ++i){
        for(int j = 0; j < parts_count; ++j){
            Point P_in(P1.x + i*h1_small, P1.y + j*h2_small);
            if(isBelongsToD(P_in)){
                n = n + 1;
            }
        }

    }

    return (length_X * length_Y) * (n / N) / (h1 * h2);    
}

double F(int i, int j)
{
    return S(i, j);
}

GridFunction A_operator(const GridFunction& i_w)
{
    GridFunction result;
    int i;
    #pragma omp parallel default(shared) private(i) 
    {
        int id = omp_get_thread_num();
        int numt = omp_get_num_threads();
        for(i = id + 1; i <= M - 1; i = i + numt)
            for(int j = 1; j <= N - 1; ++j)
            {
                result.values[i][j] = (a(i + 1, j) * i_w.w_x(i, j) - a(i, j) * i_w.w_x(i - 1, j)) / (-h1) 
                - (b(i, j + 1) * i_w.w_y(i, j) - b(i, j) * i_w.w_y(i, j - 1)) / h2;
            }
    }
    return result;
}

GridFunction B_()
{
    GridFunction B;
    int i;
    #pragma omp parallel default(shared) private(i)
    {
        int id = omp_get_thread_num();
        int numt = omp_get_num_threads();
        for(i = id + 1; i <= M - 1; i = i + numt)
            for(int j = 1; j <= N - 1; ++j)
                B.values[i][j] = F(i, j);
    }
    return B;
}

GridFunction r(const GridFunction& w)
{
    static GridFunction B = B_();
    return A_operator(w) - B;
}

double tau(const GridFunction& w)
{
    auto r_ = r(w);
    auto A_r = A_operator(r_);
    double E_A_r = E_norm(A_r);
    return scalar(A_r, r_) / (E_A_r * E_A_r);
}

int main()
{   

    GridFunction w, w_prev;


    double z;
    auto start = omp_get_wtime();

    do
    {   
        w_prev = w;
        w = w_prev - tau(w_prev) * r(w_prev);
        z = E_norm(w - w_prev);
        cout << z << endl;
    }
    while(z >= delta);

    cout << endl;
    cout << "Time = " << (double)(omp_get_wtime() - start) << " sec" <<  endl;
    
    //record u to txt 
    ofstream out;
    out.open("out.txt");
    if(out.is_open())
        w.print(out);
    out.close();
    return 0;
}
