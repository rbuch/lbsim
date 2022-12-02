#ifndef KDLB_H
#define KDLB_H

#include "TreeStrategyBase.h"
#include "kd.h"
#include <array>
#include <iostream>

#include <omp.h>

namespace TreeStrategy
{
template <typename O, typename P, typename S, typename T>
class BaseKdLB : public Strategy<O, P, S>
{
public:
  BaseKdLB() = default;
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    auto objsIter = objs.begin();
    T* tree = nullptr;
    for (int i = 0; i < procs.size() && objsIter != objs.end(); i++, objsIter++)
    {
      solution.assign(*objsIter, procs[i]);
      tree = T::insert(tree, procs[i]);
    }

    for (; objsIter != objs.end(); objsIter++)
    {
      auto proc = *(T::findMinNorm(tree, *objsIter));
      tree = T::remove(tree, proc);
      solution.assign(*objsIter, proc);
      tree = T::insert(tree, proc);
    }
    T::freeAll();
  }
};

template <typename O, typename P, typename S, typename T>
class BaseKdLBObjNorm : public Strategy<O, P, S>
{
public:
  BaseKdLBObjNorm() = default;
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    auto objsIter = objs.begin();
    T* tree = nullptr;
    for (int i = 0; i < procs.size() && objsIter != objs.end(); i++, objsIter++)
    {
      solution.assign(*objsIter, procs[i]);
      tree = T::insert(tree, procs[i]);
    }

    for (; objsIter != objs.end(); objsIter++)
    {
      auto proc = *(T::findMinNormObjNorm(tree, *objsIter));
      tree = T::remove(tree, proc);
      solution.assign(*objsIter, proc);
      tree = T::insert(tree, proc);
    }

    T::freeAll();
  }
};

  template <typename O, typename P, typename S, typename T>
class BaseKdLBObjNormEarly : public Strategy<O, P, S>
{
  private:
    void updateMax(std::array<KDFloatType, O::dimension>& maxLoads, const P& proc)
    {
      for (int i = 0; i < O::dimension; i++)
      {
        maxLoads[i] = std::max(maxLoads[i], proc[i]);
      }
    }

public:
  BaseKdLBObjNormEarly() = default;
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    std::array<KDFloatType, O::dimension> maxLoads = {0};
    auto objsIter = objs.begin();
    T* tree = nullptr;
    for (int i = 0; i < procs.size() && objsIter != objs.end(); i++, objsIter++)
    {
      solution.assign(*objsIter, procs[i]);
      updateMax(maxLoads, procs[i]);
      tree = T::insert(tree, procs[i]);
    }

    for (; objsIter != objs.end(); objsIter++)
    {
      auto proc = *(T::findMinNormObjNormEarly(tree, *objsIter, maxLoads));
      tree = T::remove(tree, proc);
      solution.assign(*objsIter, proc);
      updateMax(maxLoads, proc);
      tree = T::insert(tree, proc);
    }

    T::freeAll();
  }
};

template <typename O, typename P, typename S, typename T>
class BaseKdLBPareto : public Strategy<O, P, S>
{
public:
  BaseKdLBPareto() = default;
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    auto objsIter = objs.begin();
    T* tree = nullptr;
    for (int i = 0; i < procs.size() && objsIter != objs.end(); i++, objsIter++)
    {
      solution.assign(*objsIter, procs[i]);
      tree = T::insert(tree, procs[i]);
    }

    T* paretoFrontier = T::getParetoFrontier(tree);
    std::array<KDFloatType, O::dimension> lastRemovedProc = {0};
    for (; objsIter != objs.end(); objsIter++)
    {
      auto proc = *(T::findMinNormObjNorm(paretoFrontier, *objsIter));

      for (int i = 0; i < O::dimension; i++)
        lastRemovedProc[i] = proc[i];

      tree = T::remove(tree, proc);
      paretoFrontier = T::remove(paretoFrontier, proc);

      solution.assign(*objsIter, proc);
      tree = T::insert(tree, proc);
      if (paretoFrontier != nullptr)
      {
        const auto nn = T::getNN(paretoFrontier, lastRemovedProc);
        paretoFrontier = T::updateParetoFrontier(tree, lastRemovedProc, paretoFrontier, nn);
      }
      else
        paretoFrontier = T::getParetoFrontier(tree);
    }

    T::freeAll();
  }
};

template <typename O, typename P, typename S, typename T>
class BaseKdLBParallel : public Strategy<O, P, S>
{
public:
  BaseKdLBParallel() = default;
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    auto objsIter = objs.begin();
    T* tree = nullptr;
    for (int i = 0; i < procs.size() && objsIter != objs.end(); i++, objsIter++)
    {
      solution.assign(*objsIter, procs[i]);
      tree = T::insert(tree, procs[i]);
    }

    const auto numProcs = 2;//omp_get_num_procs();
    omp_set_num_threads(numProcs);
    std::cout << "Running with " << numProcs << " threads" << std::endl;
    for (;objsIter != objs.end();)
    {
      std::vector<P> candidates(numProcs);

      // Look for candidate PEs in parallel
#pragma omp parallel num_threads(numProcs) shared(candidates)
      {
        const auto threadId = omp_get_thread_num();
        if (std::distance(objs.begin(), threadId + objsIter) < objs.size())
        {
          candidates[threadId] = *(T::findMinNorm(tree, *(objsIter + threadId)));
        }
      }

      const auto start = std::distance(objs.begin(), objsIter);
      const auto end = std::min((size_t)objs.size(), (size_t)std::distance(objs.begin(), objsIter + candidates.size()));
      std::vector<int> usedIds;
      usedIds.reserve(candidates.size());
      for (int i = start; i < end; i++)
      {
        const auto offset = i - start;
        auto candidate = candidates[offset];

        // If we've assigned to this PE already, then look for the best candidate again
        if (std::find(usedIds.begin(), usedIds.end(), candidate.id) != usedIds.end())
        {
          candidate = *(T::findMinNorm(tree, *objsIter));
        }

        tree = T::remove(tree, candidate);
        solution.assign(*objsIter, candidate);
        tree = T::insert(tree, candidate);

        usedIds.push_back(candidate.id);

        objsIter++;
      }
    }
    T::freeAll();
  }
};

template <typename O, typename P, typename S>
class KdLB : public BaseKdLB<O, P, S, KDNode<P>>
{
};

template <int Exp>
class KdExpLB
{
public:
  template <typename O, typename P, typename S>
  class KdLB : public BaseKdLB<O, P, S, KDNode<P, Exp>>
  {
  };
};

template <typename O, typename P, typename S>
class RKdLB : public BaseKdLB<O, P, S, RKDNode<P>>
{
};

template <int Exp>
class RKdExpLB
{
public:
  template <typename O, typename P, typename S>
  class RKdLB : public BaseKdLB<O, P, S, RKDNode<P, Exp>>
  {
  };
};

template <int Exp>
class RKdExpLBObjNorm
{
public:
  template <typename O, typename P, typename S>
  class RKdLB : public BaseKdLBObjNorm<O, P, S, RKDNode<P, Exp>>
  {
  };
};

template <int Exp>
class RKdExpLBObjNormEarly
{
public:
  template <typename O, typename P, typename S>
  class RKdLB : public BaseKdLBObjNormEarly<O, P, S, RKDNode<P, Exp>>
  {
  };
};

template <int Exp>
class RKdExpLBPareto
{
public:
  template <typename O, typename P, typename S>
  class RKdLB : public BaseKdLBPareto<O, P, S, RKDNode<P, Exp>>
  {
  };
};

template <int Exp>
class RKdExpLBParallel
{
public:
  template <typename O, typename P, typename S>
  class RKdLB : public BaseKdLBParallel<O, P, S, RKDNode<P, Exp>>
  {
  };
};


template <typename O, typename P, typename S, typename T>
class BaseKdConstraintLB : public Strategy<O, P, S>
{
  private:
    const std::vector<KDFloatType>& constraints;

public:
  BaseKdConstraintLB(const std::vector<KDFloatType>& constraints) : constraints(constraints) {}

  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    const int numConstraints = constraints.size();

    // TODO: Sort only by non-constraint dimensions
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    auto objsIter = objs.begin();
    T* tree = nullptr;
    for (int i = 0; i < procs.size() && objsIter != objs.end(); i++, objsIter++)
    {
      solution.assign(*objsIter, procs[i]);
      tree = T::insert(tree, procs[i], numConstraints);
    }

    for (; objsIter != objs.end(); objsIter++)
    {
      auto proc = *(T::findMinNormConstraints(tree, *objsIter, constraints));

      tree = T::remove(tree, proc);
      //CkAssert(T::checkTree(tree);
      solution.assign(*objsIter, proc);
      tree = T::insert(tree, proc, numConstraints);
      //CkAssert(T::checkTree(tree);
    }
  }
};

template <typename O, typename P, typename S>
class KdConstraintLB : public BaseKdConstraintLB<O, P, S, KDNode<P>>
{
public:
  using BaseKdConstraintLB<O, P, S, KDNode<P>>::BaseKdConstraintLB;
};

template <typename O, typename P, typename S>
class RKdConstraintLB : public BaseKdConstraintLB<O, P, S, RKDNode<P>>
{
public:
  using BaseKdConstraintLB<O, P, S, RKDNode<P>>::BaseKdConstraintLB;
};

}  // namespace TreeStrategy

#endif /* KDLB_H */
