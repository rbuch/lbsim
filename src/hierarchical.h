#ifndef HIERARCHICAL_H
#define HIERARCHICAL_H

#include "TreeStrategyBase.h"
#include <array>
#include <iostream>

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
      std::vector<P> nodes(std::ceil(procs.size() / (float)pesPerNode));

      std::array<LoadFloatType, O::dimension + 1> empty = {0};
      for (int i = 0; i < nodes.size(); i++) {
        ptr(nodes[i])->populate(i, empty.data(), nullptr);
      }

      HierSolution<O, P> rootSol(nodes.size());
      T<O, P, HierSolution<O, P>> rootStrategy;
      rootStrategy.solve(objs, nodes, rootSol, false);

      T<O, P, S> leafStrategy;
      for (int i = 0; i < rootSol.objs.size(); i++) {
        std::vector<P> procSubset(((i + 1) * pesPerNode <= procs.size())
                                      ? pesPerNode
                                      : procs.size() % pesPerNode);
        for (int j = 0; j < procSubset.size(); j++) {
          ptr(procSubset[j])
              ->populate(j + i * pesPerNode, empty.data(), nullptr);
        }
        leafStrategy.solve(rootSol.objs[i], procSubset, solution, false);
      }
    }
  };
};
} // namespace TreeStrategy

#endif /* HIERARCHICAL_H */
