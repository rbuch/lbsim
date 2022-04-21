#ifndef BISECT_H
#define BISECT_H

#include "TreeStrategyBase.h"
#include "pheap.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <random>
#include <vector>

namespace TreeStrategy
{
template <typename O, typename P, typename S>
class Bisect : public Strategy<O, P, S>
{
 public:
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    std::vector<ProcHeap<P>> heaps;
    for (int i = 0; i < O::dimension; i++)
    {
      heaps.push_back(ProcHeap<P>(procs, i));
    }

    std::vector<std::priority_queue<O>> objHeaps;
    for (int i = 0; i < O::dimension; i++)
    {
      auto dimComp = [=](const O& o1, const O& o2) {
	return o1.load[i] < o2.load[i];
      };
      objHeaps.push_back(std::priority_queue<O>(dimComp));
    }


    std::array<float, O::dimension> halfLoad;

    for (const auto& o : objs)
    {
      for (int i = 0; i < O::dimension; i++)
      {
        halfLoad[i] += o.load[i];
      }
    }

    for (int i = 0; i < O::dimension; i++)
    {
      halfLoad[i] /= 2;
    }

    for (const auto& o : objs)
    {
      int maxdimension = 0;
      float maxfactor = 0;
      for (int i = 0; i < O::dimension; i++)
      {
        if (o.load[i] / halfLoad[i] > maxfactor)
        {
          maxfactor = o.load[i] / halfLoad[i];
          maxdimension = i;
        }
      }
      P p = heaps[maxdimension].top();
      solution.assign(o, p);  // update solution (assumes solution updates processor load)
      for (int i = 0; i < O::dimension; i++)
      {
        heaps[i].remove(p);
        heaps[i].push(p);
      }
    }
  }
};

}  // namespace TreeStrategy

#endif /* BISECT_H */
