#ifndef KDLB_H
#define KDLB_H

#include "TreeStrategyBase.h"
#include "kd.h"
#include <array>
#include <iostream>
#include <numeric>

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

template <typename O, typename P, typename S, typename T, int THRESH=1>
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
      bool didEarly = false;
      auto proc =
          *(T::findMinNormObjNormEarly(tree, *objsIter, maxLoads, THRESH, didEarly));
      if (didEarly)
	  std::cout << std::distance(objs.begin(), objsIter) << std::endl;
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

template <int Exp, int THRESH=1>
class RKdExpLBObjNormEarly
{
public:
  template <typename O, typename P, typename S>
  class RKdLB : public BaseKdLBObjNormEarly<O, P, S, RKDNode<P, Exp>, THRESH>
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


template <typename O, typename P, typename S, typename T>
class BaseKdConstraintLB : public Strategy<O, P, S>
{
  private:
    const std::vector<KDFloatType>& constraints;
    std::vector<bool> isRelative;

public:
  BaseKdConstraintLB(const std::vector<KDFloatType>& constraints)
      : constraints(constraints)
  {
    isRelative.resize(constraints.size(), false);
  }
  BaseKdConstraintLB(const std::vector<KDFloatType>& constraints,
                     const std::vector<bool>& isRelative)
      : constraints(constraints), isRelative(isRelative)
  {
  }

  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    auto currentConstraints = constraints;
    const int numConstraints = currentConstraints.size();

    for (int c_i = 0, o_i = O::dimension - numConstraints; c_i < numConstraints;
         c_i++, o_i++)
    {
      // If the constraint is relative, then set it to the product of the user
      // specified constraint value and the mean load (sum / num PEs) of the
      // corresponding dimension of the object load vectors.
      //
      // e.g. a user specified constraint of 1.5 for 100 objects and 10 PEs
      // where each object's load vector has 1.0 in the corresponding dimension
      // means the LB will keep the resulting load on every PE <= 15 = (1.0 *
      // 100 / 10 * 1.5)
      if (isRelative[c_i])
      {
        KDFloatType sum = 0.0;
        for (const auto& obj : objs)
        {
          sum += obj.getLoad(o_i);
        }
        currentConstraints[c_i] *= sum / procs.size();
      }
    }

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
      auto proc = *(T::findMinNormConstraints(tree, *objsIter, currentConstraints));

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
