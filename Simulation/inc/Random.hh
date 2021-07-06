
#ifndef Random_HH
#define Random_HH

#include <cstdlib>
#include <random>


class Random{
public:
  static double UniformRand(){
    return dis(generator);
  }

public:
   static std::mt19937 generator;
   static std::uniform_real_distribution<double> dis;



};

#endif // Random_HH
