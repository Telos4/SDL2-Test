#ifndef CVECTOR_H
#define CVECTOR_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
typedef unsigned int number;
template <const number DIM> class CVector;
template <const number DIM>
std::ostream & operator<<(std::ostream & os, const CVector<DIM> & obj);

template <const number DIM>
std::istream & operator>>(std::istream & is, CVector<DIM> & obj);

template <const number DIM>
class CVector
{
    public:
        CVector();
        CVector(double * data, number dimension);
        CVector(const std::vector<double> & vec);
        CVector(const CVector& other);
        CVector(const std::string & filename);

        number getDimension() const;
        CVector & operator=(const CVector & other);

        CVector & operator+=(const CVector & other);
        CVector & operator-=(const CVector & other);
        CVector & operator*=(double scalar);

        CVector operator+(const CVector & other);
        CVector operator-(const CVector & other);
        CVector operator*(double scalar);

        double operator[](int index) const;
        double & operator[](int index);

        friend std::ostream & operator<< <>(std::ostream & os, const CVector & obj);
        friend std::istream & operator>> <>(std::istream & is, CVector & obj);

        double norm_1();    // 1-Norm
        double norm_2();    // 2-Norm
        double norm_inf();  // Unendlich-Norm

        void read(const std::string & filename);

        static CVector UnitVector(int k);    // Liefert den k-ten Einheitsvektor

    private:
        double m_data[DIM];    // Vektoreinträge
        static const number m_n = DIM;    // Dimension

        static int setw_size;
        static int setprecision_size;
};template <const number DIM>
int CVector<DIM>::setw_size = 10;

template <const number DIM>
int CVector<DIM>::setprecision_size = 5;
template <const number DIM>
CVector<DIM>::CVector()

{    for (int i = 0; i < m_n; i++)
        m_data[i] = 0.0;
}
template <const number DIM>
CVector<DIM>::CVector(double * data, number dimension)
{
    if (dimension >= DIM)
        exit(1);

    for (int i = 0; i < dimension; i++)
        m_data[i] = data[i];

    for (int i = dimension; i < m_n; i++)

        m_data[i] = 0.0;
}

template <const number DIM>
CVector<DIM>::CVector(const std::vector<double> & vec)

{

    int dim_vec = vec.size();

    if (dim_vec >= DIM)
        exit(1);

    for (int i = 0; i < dim_vec; i++)
        m_data[i] = vec[i];

    for (int i = dim_vec; i < m_n; i++)

        m_data[i] = 0.0;
}
template <const number DIM>
CVector<DIM>::CVector(const CVector<DIM> & other)
{
    for (int i = 0; i < m_n; i++)
        m_data[i] = other.m_data[i];
}
template <const number DIM>number CVector<DIM>::getDimension() const
{
    return m_n;
}
template <const number DIM>
CVector<DIM> & CVector<DIM>::operator=(const CVector<DIM> & rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    for (int i = 0; i < m_n; i++)
        m_data[i] = rhs.m_data[i];
    return *this;
}
template <const number DIM>
CVector<DIM> &CVector<DIM>::operator+=(const CVector<DIM> & other)
{
    if (m_n == other.m_n)
    {
        for (int i = 0; i < m_n; i++)
            m_data[i] += other.m_data[i];
    }
    else
        std::cerr << "Unterschiedliche Dimensionen!" << std::endl;
    return *this;
}

template <const number DIM>
CVector<DIM> & CVector<DIM>::operator-=(const CVector<DIM> & other)
{
    if (m_n == other.m_n)
    {
        for (int i = 0; i < m_n; i++)
            m_data[i] -= other.m_data[i];
    }
    else
        std::cerr << "Unterschiedliche Dimensionen!" << std::endl;
    return *this;
}
template <const number DIM>
CVector<DIM> & CVector<DIM>::operator*=(double scalar)
{
    for (int i = 0; i < m_n; i++)
        m_data[i] *= scalar;
    return *this;
}
template <const number DIM>
CVector<DIM> CVector<DIM>::operator+(const CVector & other)
{
    CVector<DIM> newVector(*this);
    newVector += other;
    return newVector;
}
template <const number DIM>
CVector<DIM> CVector<DIM>::operator-(const CVector<DIM> & other)
{
    CVector<DIM> newVector(*this);
    newVector -= other;
    return newVector;
}
template <const number DIM>
CVector<DIM> CVector<DIM>::operator*(double scalar)
{
    CVector<DIM> newVector(*this);
    newVector *= scalar;
    return newVector;
}
template <const number DIM>
CVector<DIM> operator*(double scalar, const CVector<DIM> & other)
{
    CVector<DIM> newVector(other);
    newVector *= scalar;
    return newVector;
}

template <const number DIM>
double CVector<DIM>::operator[](int index) const
{
    if (index < 0 || index >= m_n)
    {
        std::cerr << "Indexfehler!" << std::endl;
        return m_data[0];
    }
    return m_data[index];
}
template <const number DIM>
double & CVector<DIM>::operator[](int index)
{
    if (index < 0 || index >= m_n)
    {
        std::cerr << "Indexfehler!" << std::endl;
        return m_data[0];
    }
    return m_data[index];
}
template <const number DIM>
std::ostream & operator<<(std::ostream & os, const CVector<DIM> & obj)
{
    os << std::endl;
    os << std::setprecision(CVector<DIM>::setprecision_size);
    for (int i = 0; i < obj.m_n; i++)
        os << std::setw(CVector<DIM>::setw_size) << obj.m_data[i] << std::endl << std::endl;
    return os;
}
template <const number DIM>
std::istream & operator>>(std::istream & is, CVector<DIM> & obj)
{    for (int i = 0; i < obj.m_n; i++)
    {
        std::cin >> obj.m_data[i];
    }

    return is;
}
template <const number DIM>
double CVector<DIM>::norm_1()
{
    double res = 0;
    for (int i = 0; i < m_n; i++)
        res += fabs(m_data[i]);
    return res;
}
template <const number DIM>
double CVector<DIM>::norm_2()
{
    double res = 0;
    for (int i = 0; i < m_n; i++)
        res += m_data[i]*m_data[i];
    return sqrt(res);
}
template <const number DIM>
double CVector<DIM>::norm_inf()
{
    double res = 0.0;
    for (int i = 0; i < m_n; i++)
        if (fabs(m_data[i] > res))
            res = fabs(m_data[i]);
    return res;
}
template <const number DIM>
CVector<DIM> CVector<DIM>::UnitVector(int k)
{
    CVector<DIM> unitvector;
    if (k < DIM)
        unitvector.m_data[k] = 1.0;
    return unitvector;
}


#endif // CVECTOR_H
