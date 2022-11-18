#ifndef HIERARCHICAL_H
#define HIERARCHICAL_H

#include "TreeStrategyBase.h"
#include <array>
#include <chrono>
#include <iostream>
#include <limits>

namespace TreeStrategy {
template <template <typename, typename, typename, typename...> class T>
class HierarchicalLB {
public:
  template <typename O, typename P, typename S>
  class LB : public Strategy<O, P, S> {
  private:
    const int pesPerNode = 100;

    template <typename _O, typename _P> class HierSolution {
      static constexpr auto dimension = _O::dimension;

    public:
      std::vector<std::vector<_O>> objs;

      HierSolution(size_t numProcs) { objs.resize(numProcs); }

      inline void assign(const _O *o, _P *p) {
        p->assign(*o);
        objs[p->id].push_back(*o);
      }

      inline void assign(const _O &o, _P &p) { assign(&o, &p); }
    };

  public:
    LB() = default;

    void solve(std::vector<O> &objs, std::vector<P> &procs, S &solution,
               bool objsSorted) {

      // First, LB across "nodes"
      const auto rootStart = std::chrono::steady_clock::now();
      std::vector<P> nodes(std::ceil(procs.size() / (float)pesPerNode));

      std::array<LoadFloatType, O::dimension + 1> empty = {0};
      for (int i = 0; i < nodes.size(); i++) {
        ptr(nodes[i])->populate(i, empty.data(), nullptr);
      }

      HierSolution<O, P> rootSol(nodes.size());
      T<O, P, HierSolution<O, P>> rootStrategy;
      rootStrategy.solve(objs, nodes, rootSol, false);
      const auto rootEnd = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed_seconds = rootEnd - rootStart;
      std::cout << "\tTime for root LB: " << elapsed_seconds.count() << std::endl;

      // Now, LB within each "node"
      T<O, P, S> leafStrategy;
      double leafMax = 0;
      for (int i = 0; i < rootSol.objs.size(); i++) {
        const auto leafStart = std::chrono::steady_clock::now();
        std::vector<P> procSubset(((i + 1) * pesPerNode <= procs.size())
                                      ? pesPerNode
                                      : procs.size() % pesPerNode);
        for (int j = 0; j < procSubset.size(); j++) {
          ptr(procSubset[j])
              ->populate(j + i * pesPerNode, empty.data(), nullptr);
        }
        leafStrategy.solve(rootSol.objs[i], procSubset, solution, false);
        const auto leafEnd = std::chrono::steady_clock::now();
        const auto leafElapsed = std::chrono::duration<double>(leafEnd - leafStart).count();
        leafMax = std::max(leafElapsed, leafMax);
      }
      std::cout << "\tMax time for leaf LB: " << leafMax << std::endl;
      std::cout << "\tAdj. time for hierarchical LB: " << leafMax + elapsed_seconds.count() << std::endl;
    }
  };
};
} // namespace TreeStrategy

#endif /* HIERARCHICAL_H */
