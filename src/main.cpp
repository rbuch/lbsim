#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>
#include <array>
#include <iomanip>
#include <iostream>
#include <cassert>

#define CkAssert assert

constexpr void CkAssertMsg(bool cond, const char* msg)
{
  if (!cond)
  {
    std::cout << msg << std::endl;
    assert(cond);
  }
}

#include "vt_parse.h"
#include "TreeStrategyBase.h"
#include "greedy.h"
#include "kdlb.h"
#include "metislb.h"
#include "orb.h"
#include "fuzzyorb.h"

constexpr auto dimension = 2;
constexpr auto posDimension = 3;
using O = TreeStrategy::ObjPos<dimension>;
using P = TreeStrategy::Proc<dimension, false>;

template <typename O, typename P>
class Solution
{
  static constexpr auto dimension = O::dimension;
public:
  std::vector<std::array<double, dimension>> loads;

  Solution(size_t numProcs) { loads.resize(numProcs); }

  inline void assign(const O* o, P* p)
  {
    p->assign(o);
    auto& procLoad = loads[p->id];
    for (int i = 0; i < dimension; i++)
    {
      procLoad[i] += o->getLoad(i);
    }
  }

  inline void assign(const O& o, P& p) { assign(&o, &p); }
};

std::vector<std::array<double, dimension + 1>> generate(const size_t numObjs,
                                                        const size_t numProcs,
                                                        const int seed)
{
  std::random_device rd;
  std::default_random_engine engine{rd()};
  if (seed >= 0)
    engine.seed(seed);
  std::exponential_distribution<double> expo(0.15);
  std::normal_distribution<double> normal(10, 3);

  std::vector<std::array<double, dimension + 1>> loads;
  loads.reserve(numObjs);
  for (size_t i = 0; i < numObjs; i++)
  {
    std::array<double, dimension + 1> load = {0};
    for (int j = 1; j < load.size(); j++)
    {
      if (j & 1)
        load[j] = expo(engine);
      else
        load[j] = std::max(0.0, normal(engine));
      load[0] += load[j];
    }
    loads.push_back(load);
  }
  std::array<double, dimension> min = {std::numeric_limits<double>::max()}, max = {0},
                                sum = {0};
  for (const auto& load : loads)
  {
    for (int i = 0; i < dimension; i++)
    {
      min[i] = std::min(min[i], load[i+1]);
      max[i] = std::max(max[i], load[i+1]);
      sum[i] += load[i+1];
    }
  }

  std::cout << "d  min    mean    max    sum    sum per pe" << std::endl;
  for (int i = 0; i < dimension; i++)
  {
    std::cout << i << ": " <<  std::setprecision(4) << std::fixed << min[i] << " " << sum[i] / loads.size() << " "  << max[i] << " " << sum[i] << " " << sum[i] / numProcs << std::endl;
  }

  return loads;
}

std::vector<std::vector<float>> generatePositions(const size_t numObjs, const int seed)
{
  std::random_device rd;
  std::default_random_engine engine{rd()};
  if (seed >= 0)
    engine.seed(seed + 1);
  std::uniform_real_distribution<float> uniform(0, 100);


  std::vector<std::vector<float>> positions;
  positions.reserve(numObjs);
  for (size_t i = 0; i < numObjs; i++)
  {
    std::vector<float> pos;
    for (int j = 0; j < posDimension; j++)
    {
      pos.push_back(uniform(engine));
    }
    positions.push_back(pos);
  }
  std::array<float, posDimension> min = {std::numeric_limits<float>::max()}, max = {0},
                                sum = {0};
  for (const auto& p : positions)
  {
    for (int i = 0; i < posDimension; i++)
    {
      min[i] = std::min(min[i], p[i]);
      max[i] = std::max(max[i], p[i]);
      sum[i] += p[i];
    }
  }

  std::cout << "Positions" << std::endl;
  std::cout << "d  min    mean    max" << std::endl;
  for (int i = 0; i < posDimension; i++)
  {
    std::cout << i << ": " <<  std::setprecision(4) << std::fixed << min[i] << " " << sum[i] / positions.size() << " "  << max[i] << std::endl;
  }

  return positions;
}

void populate(std::vector<O>& objs, std::vector<P>& procs, const int seed)
{
  const auto objsPerProc = objs.size() / (double)procs.size();

  std::vector<std::array<double, dimension + 1>> loads =
      generate(objs.size(), procs.size(), seed);

  auto positions = generatePositions(objs.size(), seed);

  for (int i = 0; i < objs.size(); i++)
  {
    ptr(objs[i])->populate(i, loads[i].data(), i / objsPerProc);
    ptr(objs[i])->setPosition(positions[i]);
  }
  std::array<double, dimension + 1> empty = {0};
  for (int i = 0; i < procs.size(); i++)
  {
    ptr(procs[i])->populate(i, empty.data(), nullptr);
  }
}

template<template<typename, typename, typename, typename...> class T, typename O, typename P>
void testLB(std::vector<O> objs, std::vector<P> procs, std::string lb_name)
{
  constexpr int dimension = O::dimension;
  T<O, P, Solution<O, P>> strat;
  Solution<O, P> sol(procs.size());

  auto start = std::chrono::steady_clock::now();
  strat.solve(objs, procs, sol, false);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  std::cout << "Elapsed time for " << lb_name << ": "  << elapsed_seconds.count() << std::endl;
  std::array<double, dimension> maxloads = {0};
  std::array<double, dimension> totalloads = {0};
  for (const auto& loadVec : sol.loads)
  {
    for (int i = 0; i < maxloads.size(); i++)
    {
      maxloads[i] = std::max(maxloads[i], loadVec[i]);
      totalloads[i] += loadVec[i];
    }
  }
  double maxSum = 0;
  std::cout << "Maxloads: ";
  for (const auto& load : maxloads)
  {
    std::cout << load << " ";
    maxSum += load;
  }
  std::cout << "(∑=" << maxSum << ", max=" << *std::max_element(maxloads.begin(), maxloads.end()) << ")";
  std::cout << std::endl;

  double loadSum = 0;
  std::cout << "Ratio: ";
  for (int i = 0; i < dimension; i++)
  {
    std::cout << maxloads[i] / (totalloads[i] / procs.size()) << " ";
    loadSum += totalloads[i];
  }
  std::cout << "(∑=" << maxSum / (loadSum / procs.size())  << ")";
  std::cout << std::endl << std::endl;
}

template <int N>
void testLBHelper(size_t dim, const std::vector<std::vector<double>>& objLoads,
                  const std::vector<double>& bgLoads, const bool testScalar = false)
{
  if constexpr (N == 0)
    assert(false);
  else if (dim == N)
  {
    using ObjType = TreeStrategy::Obj<N>;
    using ProcType = TreeStrategy::Proc<N, false>;

    std::vector<ObjType> objs(objLoads.size());
    for (int i = 0; i < objs.size(); i++)
    {
      // TODO: Set the correct oldPe
      objs[i].populate(i, objLoads[i].data(), i / bgLoads.size());
    }
    std::vector<ProcType> procs(bgLoads.size());
    for (int i = 0; i < procs.size(); i++)
    {
      double curBgLoad = 0;//bgLoads[i] / dim;
      procs[i].populate(i, &curBgLoad, nullptr);
    }

    if (!testScalar)
    {
      testLB<TreeStrategy::Dummy, ObjType, ProcType>(objs, procs, "dummy");
      testLB<TreeStrategy::Random, ObjType, ProcType>(objs, procs, "random");
      testLB<TreeStrategy::Greedy, ObjType, ProcType>(objs, procs, "greedy");
      // testLB<TreeStrategy::GreedyNorm>(objs, procs, "greedynorm");
      // testLB<TreeStrategy::KdLB>(objs, procs, "kd");
      testLB<TreeStrategy::RKdLB, ObjType, ProcType>(objs, procs, "rkd");
      testLB<TreeStrategy::MetisLB, ObjType, ProcType>(objs, procs, "metis");
      // testLB<TreeStrategy::GreedySample>(objs, procs, "greedysample");
      // testLB<TreeStrategy::RandomScore>(objs, procs, "randomScore");
    }
    else
    {
      testLB<TreeStrategy::ScalarGreedy, ObjType, ProcType>(objs, procs, "scalargreedy");
    }
  }
  else
  {
    testLBHelper<N - 1>(dim, objLoads, bgLoads, testScalar);
  }
}

void testVtLogs(std::vector<std::vector<double>>& objLoads, std::vector<double>& bgLoads)
{
  assert(!objLoads.empty());

  std::vector<double> totalLoad(objLoads[0].size());

  const auto objDimension = objLoads[0].size() - 1;
  for (const auto& loadVec : objLoads)
  {
    assert(loadVec.size() - 1 == objDimension);
    for (int i = 0; i < totalLoad.size(); i++)
    {
      totalLoad[i] += loadVec[i];
    }
  }

  std::cout << "\n\nTotalLoad:\n";
  for (int i = 1; i < totalLoad.size(); i++)
  {
    std::cout << totalLoad[i] << " ";
  }
  std::cout << "\n\n";

  testLBHelper<15>(objDimension, objLoads, bgLoads, true);

  testLBHelper<15>(objDimension, objLoads, bgLoads);
}


int main(int argc, char* argv[])
{
  std::vector<O> objs;
  std::vector<P> procs;

  char** it = std::find(argv, argv + argc, (std::string)"-vt");
  if (it != argv + argc)
  {
    std::vector<std::vector<double>> loads;
    std::vector<double> bgloads;
    // Read vt log files
    while (++it != argv + argc)
    {
      std::cout << "Loading " << *it << std::endl;
      double bgload;
      vt::loadFile(*it, loads, bgload);
      bgloads.push_back(bgload);
    }

    testVtLogs(loads, bgloads);
    return 0;
  }

  const auto numProcs = (argc > 1) ? std::stoi(argv[1]) : 8192;
  const auto numObjs = (argc > 2) ? std::stoi(argv[2]) : 65536;
  const int seed = (argc > 3) ? std::stoi(argv[3]) : -1;

  objs.resize(numObjs);
  procs.resize(numProcs);
  populate(objs, procs, seed);

  std::cout << "Testing with " << numObjs << " objects and " << numProcs << " processors." << std::endl;

  testLB<TreeStrategy::Dummy>(objs, procs, "dummy");
  testLB<TreeStrategy::Random>(objs, procs, "random");
  testLB<TreeStrategy::Greedy>(objs, procs, "greedy");
  // testLB<TreeStrategy::GreedyNorm>(objs, procs, "greedynorm");
  // testLB<TreeStrategy::KdLB>(objs, procs, "kd");
  // testLB<TreeStrategy::RKdExpLB<1>::RKdLB>(objs, procs, "rkd1");
  testLB<TreeStrategy::RKdExpLB<2>::RKdLB>(objs, procs, "rkd2");
  testLB<TreeStrategy::RKdExpLB<4>::RKdLB>(objs, procs, "rkd4");
  testLB<TreeStrategy::RKdExpLB<8>::RKdLB>(objs, procs, "rkd8");
  testLB<TreeStrategy::RKdExpLB<16>::RKdLB>(objs, procs, "rkd16");
  testLB<TreeStrategy::MetisLB>(objs, procs, "metis");
  //testLB<TreeStrategy::GreedySample>(objs, procs, "greedysample");
  //testLB<TreeStrategy::RandomScore>(objs, procs, "randomScore");
  testLB<TreeStrategy::ORBScalar>(objs, procs, "orbScalar");
  testLB<TreeStrategy::ORBVector>(objs, procs, "orbVector");
  testLB<TreeStrategy::FuzzyORBScalar>(objs, procs, "fuzzyOrbScalar");

  return 0;
}
