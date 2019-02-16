// function template
#include <iostream>
using namespace std;

template <class T, int b>
T GetMax (T a) {
  return a*b;
}

int main () {
  int a=2, b=3;
  
  cout << GetMax<int,5>(a);

  return 0;
}

