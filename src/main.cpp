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
#include "load_config.h"

#include "TreeStrategyBase.h"
#include "greedy.h"
#include "kdlb.h"
#include "metislb.h"
#include "orb.h"
#include "fuzzyorb.h"

#include "hierarchical.h"

constexpr auto dimension = 2;
constexpr auto posDimension = 3;
using O = TreeStrategy::ObjPos<dimension>;
using P = TreeStrategy::Proc<dimension, false>;

bool onlyHierarch = false;

template <typename O, typename P>
class Solution
{
  static constexpr auto dimension = O::dimension;
public:
  std::vector<std::array<LoadFloatType, dimension>> loads;

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

std::vector<std::array<LoadFloatType, dimension + 1>> generate(const size_t numObjs,
                                                        const size_t numProcs,
                                                        const int seed)
{
  std::random_device rd;
  std::default_random_engine engine{rd()};
  if (seed >= 0)
    engine.seed(seed);
  std::exponential_distribution<LoadFloatType> expo(0.15);
  std::normal_distribution<LoadFloatType> normal(10, 3);

  std::vector<std::array<LoadFloatType, dimension + 1>> loads;
  loads.reserve(numObjs);
  for (size_t i = 0; i < numObjs; i++)
  {
    std::array<LoadFloatType, dimension + 1> load = {0};
    for (int j = 1; j < load.size(); j++)
    {
      if (j & 1)
        load[j] = expo(engine);
      else
        load[j] = std::max((LoadFloatType)0.0, normal(engine));
      load[0] += load[j];
    }
    loads.push_back(load);
  }
  std::array<LoadFloatType, dimension> min = {std::numeric_limits<LoadFloatType>::max()}, max = {0},
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
  const LoadFloatType objsPerProc = objs.size() / (double)procs.size();

  std::vector<std::array<LoadFloatType, dimension + 1>> loads =
      generate(objs.size(), procs.size(), seed);

  auto positions = generatePositions(objs.size(), seed);

  for (int i = 0; i < objs.size(); i++)
  {
    ptr(objs[i])->populate(i, loads[i].data(), i / objsPerProc);
    ptr(objs[i])->setPosition(positions[i]);
  }
  std::array<LoadFloatType, dimension + 1> empty = {0};
  for (int i = 0; i < procs.size(); i++)
  {
    ptr(procs[i])->populate(i, empty.data(), nullptr);
  }
}

template<template<typename, typename, typename, typename...> class T, typename O, typename P>
void testLB(std::vector<O> objs, std::vector<P> procs, std::string lb_name, const std::vector<LoadFloatType>& knownLoadSum, const bool warmup = false)
{
  if (!warmup)
      std::cout << "Testing " << lb_name << ": " << std::endl;
  constexpr int dimension = O::dimension;
  T<O, P, Solution<O, P>> strat;
  Solution<O, P> sol(procs.size());

  auto start = std::chrono::steady_clock::now();
  strat.solve(objs, procs, sol, false);
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  if (!warmup)
      std::cout << "Elapsed time for " << lb_name << ": "  << elapsed_seconds.count() << std::endl;
  std::array<LoadFloatType, dimension> maxloads = {0};
  std::array<LoadFloatType, dimension> totalloads = {0};

  LoadFloatType maxDimensionSum = 0;
  for (const auto& loadVec : sol.loads)
  {
    for (int i = 0; i < maxloads.size(); i++)
    {
      maxloads[i] = std::max(maxloads[i], loadVec[i]);
      totalloads[i] += loadVec[i];
    }

    maxDimensionSum += *std::max_element(loadVec.begin(), loadVec.end());
  }
  LoadFloatType maxSum = 0;
  if (!warmup)
      std::cout << "Maxloads: ";
  for (const auto& load : maxloads)
  {
    if (!warmup)
	std::cout << load << " ";
    maxSum += load;
  }
  if (!warmup)
  {
    std::cout << "(∑=" << maxSum
              << ", max=" << *std::max_element(maxloads.begin(), maxloads.end()) << ")";
    std::cout << std::endl;
  }

  LoadFloatType loadSum = 0;
  if (!warmup)
    std::cout << "Ratio: ";
  for (int i = 0; i < dimension; i++)
  {
    if (!warmup)
      std::cout << maxloads[i] / (totalloads[i] / procs.size()) << " ";
    loadSum += totalloads[i];
  }
  if (!warmup)
  {
    std::cout << "(∑=" << maxSum / (loadSum / procs.size()) << ", max="
              << *std::max_element(maxloads.begin(), maxloads.end()) /
                     (maxDimensionSum / procs.size())
              << ")";
    std::cout << std::endl << std::endl;
  }

  // Validate results
  bool validated = true;
  for (int i = 0; i < totalloads.size(); i++)
  {
    //CkAssertMsg(totalloads[i] != knownLoadSum[i], "Load does not validate!");
    constexpr auto EPSILON = 1.0e-5;
    if (std::fabs(totalloads[i] - knownLoadSum[i + 1]) >= EPSILON) {
      std::printf("%d: known: %f, calc: %f\n", i, knownLoadSum[i + 1], totalloads[i]);
      validated = false;
    }
  }
  if (!validated)
  {
    std::abort();
  }
}

template <int N>
void testLBHelper(size_t dim, const std::vector<std::vector<LoadFloatType>>& objLoads,
                  const std::vector<LoadFloatType>& bgLoads,
                  const std::vector<LoadFloatType>& knownLoadSum,
                  const bool testScalar = false)
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
      objs[i].populate(i, objLoads[i].data(), (int)(((float)i / objs.size())  * bgLoads.size()));
    }
    std::vector<ProcType> procs(bgLoads.size());
    for (int i = 0; i < procs.size(); i++)
    {
      LoadFloatType curBgLoad = 0;//bgLoads[i] / dim;
      procs[i].populate(i, &curBgLoad, nullptr);
    }

    if (!testScalar)
    {
      //testLB<TreeStrategy::Dummy, ObjType, ProcType>(objs, procs, "dummy");
      //testLB<TreeStrategy::Random, ObjType, ProcType>(objs, procs, "random");
      // testLB<TreeStrategy::Greedy, ObjType, ProcType>(objs, procs, "greedy", knownLoadSum);
      // testLB<TreeStrategy::GreedyNorm>(objs, procs, "greedynorm");
      // testLB<TreeStrategy::KdLB>(objs, procs, "kd");
      if (!onlyHierarch) {
        // testLB<TreeStrategy::RKdExpLB<2>::RKdLB, ObjType, ProcType>(objs, procs, "rkd2",
        //                                                             knownLoadSum);
        // testLB<TreeStrategy::RKdExpLBObjNorm<2>::RKdLB, ObjType, ProcType>(
        //     objs, procs, "rkd2ObjNorm", knownLoadSum);
        // testLB<TreeStrategy::RKdExpLBObjNormEarly<2>::RKdLB, ObjType, ProcType>(
        //     objs, procs, "rkd2ObjNormEarly", knownLoadSum);
        // testLB<TreeStrategy::RKdExpLBPareto<2>::RKdLB, ObjType, ProcType>(
        //     objs, procs, "rkd2Pareto", knownLoadSum);
        // testLB<TreeStrategy::RKdExpLB<4>::RKdLB, ObjType, ProcType>(objs, procs, "rkd4",
        //                                                             knownLoadSum);
        // testLB<TreeStrategy::RKdExpLBObjNorm<4>::RKdLB, ObjType, ProcType>(
        //     objs, procs, "rkd4ObjNorm", knownLoadSum);
	  testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 1>::RKdLB, ObjType, ProcType>(
										     objs, procs, "rkd4ObjNormEarly-1", knownLoadSum, true);
	  testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 1>::RKdLB, ObjType, ProcType>(
										     objs, procs, "rkd4ObjNormEarly-1", knownLoadSum, true);

        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 1>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-1", knownLoadSum);
        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 2>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-2", knownLoadSum);
        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 3>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-3", knownLoadSum);
        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 4>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-4", knownLoadSum);
        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 5>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-5", knownLoadSum);
	testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 6>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-6", knownLoadSum);
        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 7>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-7", knownLoadSum);
        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 8>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-8", knownLoadSum);
        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 9>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-9", knownLoadSum);
        testLB<TreeStrategy::RKdExpLBObjNormEarly<4, 10>::RKdLB, ObjType, ProcType>(
            objs, procs, "rkd4ObjNormEarly-10", knownLoadSum);

        // testLB<TreeStrategy::RKdExpLBPareto<4>::RKdLB, ObjType, ProcType>(
        //     objs, procs, "rkd4Pareto", knownLoadSum);
        // testLB<TreeStrategy::RKdExpLB<8>::RKdLB, ObjType, ProcType>(objs, procs,
        // "rkd8"); testLB<TreeStrategy::RKdExpLB<16>::RKdLB, ObjType, ProcType>(objs,
        // procs, "rkd16"); testLB<TreeStrategy::RKdExpLB<100>::RKdLB, ObjType,
        // ProcType>(objs, procs, "rkdInf");
        // testLB<TreeStrategy::MetisLB, ObjType, ProcType>(objs, procs, "metis",
        //                                                  knownLoadSum);
        // testLB<TreeStrategy::GreedySample>(objs, procs, "greedysample");
        // testLB<TreeStrategy::RandomScore>(objs, procs, "randomScore");
      }
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      //                                     8>::LB,
      //        ObjType, ProcType>(objs, procs, "hierarch<8>", knownLoadSum);
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      //                                     16>::LB,
      //        ObjType, ProcType>(objs, procs, "hierarch<16>", knownLoadSum);
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      //                                     32>::LB,
      //        ObjType, ProcType>(objs, procs, "hierarch<32>", knownLoadSum);
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      //                                     64>::LB,
      //        ObjType, ProcType>(objs, procs, "hierarch<64>", knownLoadSum);
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      // 					  128>::LB,
      //        ObjType, ProcType>(objs, procs, "hierarch<128>", knownLoadSum);
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      // 					  256>::LB,
      //        ObjType, ProcType>(objs, procs, "hierarch<256>", knownLoadSum);
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      // 					  512>::LB,
      //        ObjType, ProcType>(objs, procs, "hierarch<512>", knownLoadSum);
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      // 					  1024>::LB,
      //        ObjType, ProcType>(objs, procs, "hierarch<1024>", knownLoadSum);
      // if (procs.size() >= 2048)
      // 	  testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      // 					      2048>::LB,
      // 		 ObjType, ProcType>(objs, procs, "hierarch<2048>", knownLoadSum);
      // if (procs.size() >= 4096)
      // 	  testLB<TreeStrategy::HierarchicalLB<TreeStrategy::RKdExpLB<4>::RKdLB,
      // 					      4096>::LB,
      // 		 ObjType, ProcType>(objs, procs, "hierarch<4096>", knownLoadSum);

      //testLB<TreeStrategy::HierarchicalLB<TreeStrategy::Greedy>::LB, ObjType, ProcType>(objs, procs, "hierarch");
    }
    else
    {
      // testLB<TreeStrategy::ScalarGreedy, ObjType, ProcType>(objs, procs, "scalargreedy",
      //                                                       knownLoadSum, true);
      // testLB<TreeStrategy::ScalarGreedy, ObjType, ProcType>(objs, procs, "scalargreedy",
      //                                                       knownLoadSum);

      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::ScalarGreedy, 128>::LB, ObjType,
      //        ProcType>(objs, procs, "sgreedy<128>", knownLoadSum, true);
      // testLB<TreeStrategy::HierarchicalLB<TreeStrategy::ScalarGreedy, 128>::LB, ObjType,
      //        ProcType>(objs, procs, "sgreedy<128>", knownLoadSum);
    }
  }
  else
  {
    testLBHelper<N - 1>(dim, objLoads, bgLoads, knownLoadSum, testScalar);
  }
}

void testFromLogs(std::vector<std::vector<LoadFloatType>>& objLoads, std::vector<LoadFloatType>& bgLoads)
{
  assert(!objLoads.empty());

  std::cout << std::endl
            << "Testing with " << objLoads.size() << " objects and " << bgLoads.size()
            << " processors." << std::endl;

  std::vector<LoadFloatType> totalLoad(objLoads[0].size());

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

  testLBHelper<15>(objDimension, objLoads, bgLoads, totalLoad, true);

  testLBHelper<15>(objDimension, objLoads, bgLoads, totalLoad);
}


int main(int argc, char* argv[])
{
  std::vector<O> objs;
  std::vector<P> procs;

  char** hierarchIt = std::find(argv, argv + argc, (std::string)"-onlyhierarch");
  if (hierarchIt != argv + argc)
      onlyHierarch = true;

  char** it = std::find(argv, argv + argc, (std::string)"-vt");
  if (it != argv + argc)
  {
    char** phaseIt = std::find(argv, argv + argc, (std::string)"-phase");
    const int phase = (phaseIt != argv + argc) ? std::atoi(*(++phaseIt)) : -1;

    std::cout << "Using phase " << phase << " from VT log files:"  << std::endl;

    std::vector<std::vector<LoadFloatType>> loads;
    std::vector<LoadFloatType> bgloads;
    // Read vt log files
    while (++it != argv + argc)
    {
      std::cout << "Loading " << *it << std::endl;
      LoadFloatType bgload;
      vt::loadFile(*it, loads, bgload, phase);
      bgloads.push_back(bgload);
    }

    testFromLogs(loads, bgloads);
    return 0;
  }

  const auto numProcs = (argc > 1) ? std::stoi(argv[1]) : 8192;
  const auto numObjs = (argc > 2) ? std::stoi(argv[2]) : 65536;
  const int seed = (argc > 3) ? std::stoi(argv[3]) : -1;

  objs.resize(numObjs);
  procs.resize(numProcs);

  it = std::find(argv, argv + argc, (std::string)"-json");
  if (it != argv + argc)
  {
    std::vector<std::vector<LoadFloatType>> objs;
    conf::loadFile(*(++it), objs, numObjs, seed);
    std::vector<LoadFloatType> bgLoads(numProcs, 0);
    testFromLogs(objs, bgLoads);
    return 0;
  }
  else
  {
    populate(objs, procs, seed);
  }
/*
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
  testLB<TreeStrategy::RKdExpLB<100>::RKdLB>(objs, procs, "rkdInf");
  testLB<TreeStrategy::MetisLB>(objs, procs, "metis");
  //testLB<TreeStrategy::GreedySample>(objs, procs, "greedysample");
  //testLB<TreeStrategy::RandomScore>(objs, procs, "randomScore");
  testLB<TreeStrategy::ORBScalar>(objs, procs, "orbScalar");
  testLB<TreeStrategy::ORBVector>(objs, procs, "orbVector");
  testLB<TreeStrategy::FuzzyORBScalar>(objs, procs, "fuzzyOrbScalar");
*/
  return 0;
}
