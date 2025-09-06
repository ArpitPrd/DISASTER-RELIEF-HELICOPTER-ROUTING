#include <iostream>
#include <vector>
using namespace std;



int main() {
    vector<int> a {1, 3, 4};
    vector<int> b=a;
    
    b[2]=3;

    for (int i=0; i<3; i++) cout << a[i] << " ";
    cout << endl;

    for (int i=0; i<3; i++) cout << b[i] << " ";
    cout << endl;
}