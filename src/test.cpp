#include <iostream>
#include <cmath>
#include <stdio.h>
using namespace std;
int main()
{
  double velocity = 0.6;
  double yaw_rate = 0.2;
  double dist_theta = 0.3;
  double delta_t = 2;
  double dist_x = 1;
  double result =0.0;
  
 result = dist_x + (velocity/yaw_rate) * (sin(dist_theta) + yaw_rate * delta_t) - sin(dist_theta);
 cout<<result;
}

