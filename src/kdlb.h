#ifndef KDLB_H
#define KDLB_H

#include "TreeStrategyBase.h"
#include "kd.h"
#include <iostream>

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
