#include <iostream>
#include "CVector.h"

using namespace std;

int main()
{
    CVector<3> v1, v2;

    cout << v1;
    v1[0] = 1.0;
    v1[1] = -1.0;
    v1[2] = 3.0;

    v2[2] = -7.0;

    cout << 3 * v1 + v2;
    cout << v1.getDimension();

    CVector<3> v3 = v1 - v2;
    cout << v3;
    v3 -= v1;
    cout << v3;

    v3 += v2 + 4 * v1;
    cout << v3;

    int k;
    CVector<5> v4;
    cin >> v4;

    cout << v4;

    return 0;
}
