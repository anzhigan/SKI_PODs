#include <iostream>
#include <cmath>
#include <array>
#include <time.h>
#include <omp.h>
#include <fstream>

#include "mpi.h"

struct Point{
    Point(double i_x = 0.0, double i_y = 0.0): x(i_x), y(i_y)
    {
    }

    double x;
    double y;
};

Point A(-2, 0), B(-1, 1), C(1, 1), D(2, 0), E(0, -2);
Point P_A(-3, -3), P_B(-3, 2), P_C(3, 2), P_D(3, -3);
double tooSmall = 1e-10;
int tooMuchElements = 1000;

double AB(double x) // Use only when D.x != A.x
{
    return x + 2;
}

double BC(double x) // Use only when C.x != B.x
{
    return 1; 
}

double CD(double x) // Use only when C.x != B.x
{
    return -x + 2; 
}

double ED(double x) // Use only when C.x != B.x
{
    return x - 2; 
}

double AE(double x) // Use only when C.x != B.x
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


const int M = 80, N = 80; // Grid' sizes
const double h1 = (P_D.x - P_A.x) / M, h2 = (P_B.y - P_A.y) / N;
const double epsilon = h1, delta = 0.000001;

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

    void print(std::ostream& out) const 
    {
        for(const auto& value: values)
        {
            for(const auto& w: value)
                out << w << " ";
            out << std::endl;                
        }
    }

    std::array<std::array<double, N + 1>, M + 1> values;


};

GridFunction operator+(const GridFunction& i_l, const GridFunction& i_r)
{
    GridFunction result;
    int i;
    #pragma omp parallel default(shared) private(i)
    {
        int id = omp_get_thread_num();
        int numt = omp_get_num_threads();
        for(i = id + 1; i <= M - 1; i = i + numt)
            for(int j = 1; j <= N - 1; ++j)
                result.values[i][j] = i_l.values[i][j] + i_r.values[i][j];
    }
    return result;
}


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
    double x_m_h = x_(i - 0.5), x_p_h = x_(i + 0.5);
    double y_m_h = y_(j - 0.5), y_p_h = y_(j + 0.5);
    Point P1(x_m_h, y_m_h), P2(x_p_h, y_m_h), P3(x_p_h, y_p_h), P4(x_m_h, y_p_h);
    int numberOfPointsInD = isBelongsToD(P1) + isBelongsToD(P2) + 
        isBelongsToD(P3) + isBelongsToD(P4);
   
    if(numberOfPointsInD == 4)
        return 1;
    if(numberOfPointsInD == 0)
        return 0;

    double parts_count = 5;

    double length_X = P2.x - P1.x;
    double length_Y = P3.y - P1.y;
    double h1_small = length_X / parts_count;
    double h2_small = length_Y / parts_count;
    double n = 0;
    double N = parts_count * parts_count;

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

GridFunction A_operator(const GridFunction& i_w, int rank, int size)
{
    int quadsPerProcess = 4 / size;
    int sendSize = ((M / 2) * (N / 2)) * quadsPerProcess;
    if(!rank && size > 1)
    {
        for(int i_num = 1; i_num < size; ++i_num)
        {
            int count = 0;
            double* sbuf = new double [sendSize + M + N];
            for(int iter = 0; iter < quadsPerProcess; ++iter)
            {
                int numOfQuad = quadsPerProcess * i_num + iter;
                for(int i = std::max(1, (numOfQuad/2)*M/2-1); i < std::min((1+numOfQuad/2)*M/2+1, M); ++i)
                    for(int j = std::max(1, (numOfQuad%2)*N/2-1); j < std::min((1+numOfQuad%2)*N/2+1, N); ++j)
                        sbuf[count++] = i_w.values[i][j];
            }
            if(sendSize + M + N < tooMuchElements)
                MPI_Send(sbuf, sendSize + M + N, MPI_DOUBLE, i_num, 1, MPI_COMM_WORLD);
            else
                for(int sendIter = 0; sendIter < (sendSize+M+N) / (M/2);++sendIter)
                    MPI_Send(sbuf+sendIter*(M/2), M/2, MPI_DOUBLE, i_num, 1, MPI_COMM_WORLD);
            delete[] sbuf;
        }
    }
    GridFunction result;
    double* send = nullptr;
    if(size > 1)
        send = new double [sendSize];
    int count = 0;
    for(int iter = 0; iter < quadsPerProcess; ++iter)
    {
        int numOfQuad = quadsPerProcess * rank + iter;
        for(int i = std::max(1, (numOfQuad/2)*M/2); i < std::min((1+numOfQuad/2)*M/2, M); ++i)
            for(int j = std::max(1, (numOfQuad%2)*N/2); j < std::min((1+numOfQuad%2)*N/2, N); ++j)
            {
                double value = (a(i + 1, j) * i_w.w_x(i, j) - a(i, j) * i_w.w_x(i - 1, j)) / (-h1) 
                - (b(i, j + 1) * i_w.w_y(i, j) - b(i, j) * i_w.w_y(i, j - 1)) / h2;
                if(!rank)
                {
                    result.values[i][j] = value;
                }
                else if(send)
                {
                    send[count] = value;
                    ++count;
                }
            }
    }
    if(!rank)
    {
        if(size > 1)
        {
            MPI_Status status;
            int resultsCount = 0; 
            auto recvSize = sendSize;
            double* recv = nullptr;
            while(true)
            {
                recv = new double[recvSize];
                if(recvSize < tooMuchElements)
                    MPI_Recv(recv, recvSize, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
                else
                {
                    MPI_Recv(recv, M/2, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
                    int source = status.MPI_SOURCE;
                    for(int recvIter = 1; recvIter < recvSize / (M/2); ++recvIter)
                        MPI_Recv(recv+recvIter*(M/2), M/2, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status);
                }
                ++resultsCount;
                count = 0;
                for(int iter = 0; iter < quadsPerProcess; ++iter)
                {
                    int numOfQuad = quadsPerProcess * status.MPI_SOURCE + iter;
                    for(int i = std::max(1, (numOfQuad/2)*M/2); i < std::min((1+numOfQuad/2)*M/2, M); ++i)
                        for(int j = std::max(1, (numOfQuad%2)*N/2); j < std::min((1+numOfQuad%2)*N/2, N); ++j)
                        {
                            result.values[i][j] = std::abs(recv[count]) > tooSmall ? recv[count]
                                                                        : 0.0;
                            ++count;
                        }
                }
                delete[] recv;
                if(resultsCount == size - 1)
                    break;
            }
        }
    }
    else
    {
        if(sendSize < tooMuchElements)
            MPI_Send(send, sendSize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        else
            for(int sendIter = 0; sendIter < sendSize / (M/2); ++sendIter)
                MPI_Send(send+sendIter*(M/2), M / 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    return result;
}

GridFunction B_right(int rank, int size)
{
    if(!rank && size > 1)
    {
        double buf;
        for(int i = 1; i < size; ++i)
            MPI_Send(&buf, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
    }
    GridFunction B;
    int quadsPerProcess = 4 / size;
    int sendSize = ((M / 2) * (N / 2)) * quadsPerProcess;
    double* send;
    if(rank)
        send = new double [sendSize];
    int count = 0;
    for(int iter = 0; iter < quadsPerProcess; ++iter)
    {
        int numOfQuad = quadsPerProcess * rank + iter;
        for(int i = std::max(1, (numOfQuad/2)*M/2); i < std::min((1+numOfQuad/2)*M/2, M); ++i)
            for(int j = std::max(1, (numOfQuad%2)*N/2); j < std::min((1+numOfQuad%2)*N/2, N); ++j)
            {
                double value = F(i, j);
                if(!rank)
                {
                    B.values[i][j] = value;
                }
                else if(send)
                {
                    send[count] = value;
                    ++count;
                }
            }
    }
    if(!rank)
    {
        if(size > 1)
        {
            MPI_Status status;
            int resultsCount = 0; 
            auto recvSize = sendSize;
            double* recv = nullptr;
            while(true)
            {
                recv = new double[recvSize];
                if(recvSize < tooMuchElements)
                    MPI_Recv(recv, recvSize, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
                else
                {
                    MPI_Recv(recv, M/2, MPI_DOUBLE, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
                    int source = status.MPI_SOURCE;
                    for(int recvIter = 1; recvIter < recvSize/(M/2); ++recvIter)
                        MPI_Recv(recv+recvIter*(M/2),M/2, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status);
                }
                ++resultsCount;
                count = 0;
                for(int iter = 0; iter < quadsPerProcess; ++iter)
                {
                    int numOfQuad = quadsPerProcess * status.MPI_SOURCE + iter;
                    for(int i = std::max(1, (numOfQuad/2)*M/2); i < std::min((1+numOfQuad/2)*M/2, M); ++i)
                        for(int j = std::max(1, (numOfQuad%2)*N/2); j < std::min((1+numOfQuad%2)*N/2, N); ++j)
                        {
                            B.values[i][j] = recv[count];
                            ++count;
                        }
                }
                delete[] recv;
                if(resultsCount == size - 1)
                    break;
            }
        }
    }
    else
    {
        if(sendSize < tooMuchElements)
            MPI_Send(send, sendSize, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        else
            for(int sendIter = 0; sendIter < sendSize/(M/2); ++sendIter)
                MPI_Send(send+sendIter*(M/2), M/2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
    return B;
}

GridFunction r(const GridFunction& i_w, int rank, int size)
{
    static GridFunction B = B_right(rank, size);
    return A_operator(i_w, rank, size) - B;
}

double tau(const GridFunction& i_r_k, int rank, int size)
{
    auto A_r_k = A_operator(i_r_k, rank, size);
    double E_A_r_k = E_norm(A_r_k);
    return scalar(A_r_k, i_r_k) / (E_A_r_k * E_A_r_k);
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    int quadsPerProcess = 4 / size;
    int recvSize = ((M / 2) * (N / 2)) * quadsPerProcess;
    double* recv = nullptr;
    if(size > 1)
        recv = new double [recvSize+M+N];
    GridFunction w, w_prev;
    double diff;
    if(!rank)
    { // main process
        auto start = MPI_Wtime();
        do
        {
            w_prev = w;
            auto r_k = r(w_prev, rank, size);
            w = w_prev - tau(r_k, rank, size) * r_k;
            diff = E_norm(w - w_prev);
            std::cout << "Difference = " << diff << std::endl;
        }
        while(diff >= delta);
        std::cout << std::endl;
        std::cout << "Computation time = " << (double)(MPI_Wtime() - start) << " sec" <<  std::endl;
        for(int i = 1; i < size; ++i)
            MPI_Send(&diff, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        std::ofstream out;
        out.open("out.txt");
        if(out.is_open())
            w.print(out);
        out.close();
    }
    else
    {
        while(true)
        {
            if(recvSize+M+N < tooMuchElements)
                MPI_Recv(recv, recvSize+M+N, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            else
                MPI_Recv(recv, M/2, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if(status.MPI_TAG == 0)
                break;
            if(status.MPI_TAG == 1)
            {
                if(recvSize+M+N >= tooMuchElements)
                    for(int recvIter = 1; recvIter < (recvSize+M+N)/(M/2); ++recvIter)
                        MPI_Recv(recv+recvIter*(M/2), M / 2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
                int quadsPerProcess = 4 / size;
                int count = 0;
                for(int iter = 0; iter < quadsPerProcess; ++iter)
                {
                    int numOfQuad = quadsPerProcess * rank + iter;
                    for(int i = std::max(1, (numOfQuad/2)*M/2-1); i < std::min((1+numOfQuad/2)*M/2+1, M); ++i)
                        for(int j = std::max(1, (numOfQuad%2)*N/2-1); j < std::min((1+numOfQuad%2)*N/2+1, N); ++j)
                        {    
                            w_prev.values[i][j] = std::abs(recv[count]) > tooSmall ? recv[count]
                                                                        : 0.0;
                            count++;
                        }
                }
                A_operator(w_prev, rank, size);
            }
            if(status.MPI_TAG == 2)
            {
                B_right(rank, size);
            }
        }
    }
    MPI_Finalize();
}