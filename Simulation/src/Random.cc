#include "Random.hh"

std::mt19937 Random::generator(123);
std::uniform_real_distribution<double> Random::dis(0.0,1.0);
