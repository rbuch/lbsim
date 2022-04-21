#ifndef GREEDY_H
#define GREEDY_H

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
class Greedy : public Strategy<O, P, S>
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

    std::array<float, O::dimension> averageLoad;

    for (const auto& o : objs)
    {
      for (int i = 0; i < O::dimension; i++)
      {
        averageLoad[i] += o.load[i];
      }
    }

    for (int i = 0; i < O::dimension; i++)
    {
      averageLoad[i] /= objs.size();
    }

    for (const auto& o : objs)
    {
      int maxdimension = 0;
      float maxfactor = 0;
      for (int i = 0; i < O::dimension; i++)
      {
        if (o.load[i] / averageLoad[i] > maxfactor)
        {
          maxfactor = o.load[i] / averageLoad[i];
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

template <typename P, typename S>
class Greedy<Obj<1>, P, S> : public Strategy<Obj<1>, P, S>
{
public:
  void solve(std::vector<Obj<1>>& objs, std::vector<P>& procs, S& solution,
             bool objsSorted)
  {
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<Obj<1>>());
    std::priority_queue<P, std::vector<P>, CmpLoadGreater<P>> procHeap(
        CmpLoadGreater<P>(), procs);

    for (const auto& o : objs)
    {
      P p = procHeap.top();
      procHeap.pop();
      solution.assign(o, p);  // update solution (assumes solution updates processor load)
      procHeap.push(p);
    }
  }
};

template <typename O, typename P, typename S>
class ScalarGreedy : public Strategy<O, P, S>
{
public:
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution,
             bool objsSorted)
  {
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());
    std::priority_queue<P, std::vector<P>, CmpLoadGreater<P>> procHeap(
        CmpLoadGreater<P>(), procs);

    for (const auto& o : objs)
    {
      P p = procHeap.top();
      procHeap.pop();
      solution.assign(o, p);  // update solution (assumes solution updates processor load)
      procHeap.push(p);
    }
  }
};

template <typename O, typename P>
struct GreedySolution
{
  std::array<float, O::dimension> maxload;

  GreedySolution()
  {
    maxload.fill(0);
  }

  inline void assign(const O& o, P& p)
  {
    ptr(p)->assign(o);
    for (int i = 0; i < O::dimension; i++)
    {
      maxload[i] = std::max(maxload[i], ptr(p)->getLoad(i));
    }
  }
};

template <typename P>
struct GreedySolution<Obj<1>, P>
{
  inline void assign(const Obj<1>& o, P& p)
  {
    ptr(p)->assign(o);
    maxload = std::max(maxload, ptr(p)->getLoad());
  }
  float maxload = 0;
};

// NOTE: this will modify order of objects in 'objs' if objsSorted is false
template <typename O, typename P>
typename std::conditional<(O::dimension > 1), std::array<float, O::dimension>,
                          float>::type
calcGreedyMaxload(std::vector<O>& objs, std::vector<P>& procs, bool objsSorted)
{
  GreedySolution<O, P> greedy_sol;
  Greedy<O, P, GreedySolution<O, P>> greedy;
  greedy.solve(objs, procs, greedy_sol, objsSorted);
  // greedy will modify my copy of processors only if they are passed as pointers
  if (std::is_pointer<P>::value)
    for (auto& p : procs) ptr(p)->resetLoad();
  return greedy_sol.maxload;
}

template <typename O, typename P, typename S>
class GreedyRefine : public Strategy<O, P, S>
{
private:
  float tolerance = 1;  // tolerance above greedy maxload (not average load!)

public:
  GreedyRefine() = default;

  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    auto M = calcGreedyMaxload(objs, procs, objsSorted);

    for (auto& value : M)
    {
      value *= tolerance;
    }

    // need custom heap that allows removal of elements from any position
    ProcHeap<P> procHeap(procs);
    P p;
    for (const auto& o : objs)
    {
      // TODO improve the case where the proc is not in my list of processors (because
      // it belongs to a foreing domain). procHeap API should return an error?
      P& oldPe = procHeap.getProc(ptr(o)->oldPe);
      if ((oldPe.id >= 0) && (oldPe.getLoad() + ptr(o)->getLoad() <= M[0]))
        p = oldPe;
      else
        p = procHeap.top();
      procHeap.remove(p);
      solution.assign(o, p);
      procHeap.push(p);
      M[0] = std::max(M[0], ptr(p)->getLoad());
    }
  }
};

template <typename P, typename S>
class GreedyRefine<Obj<1>, P, S> : public Strategy<Obj<1>, P, S>
{
private:
  float tolerance = 1;  // tolerance above greedy maxload (not average load!)

public:
  GreedyRefine() = default;
  void solve(std::vector<Obj<1>>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    auto M = calcGreedyMaxload(objs, procs, objsSorted);

    M *= tolerance;

    // need custom heap that allows removal of elements from any position
    ProcHeap<P> procHeap(procs);
    P p;
    for (const auto& o : objs)
    {
      // TODO improve the case where the proc is not in my list of processors (because
      // it belongs to a foreing domain). procHeap API should return an error?
      P& oldPe = procHeap.getProc(ptr(o)->oldPe);
      if ((oldPe.id >= 0) && (oldPe.getLoad() + ptr(o)->getLoad() <= M))
        p = oldPe;
      else
        p = procHeap.top();
      procHeap.remove(p);
      solution.assign(o, p);
      procHeap.push(p);
      M = std::max(M, ptr(p)->getLoad());
    }
  }
};

template <typename O, typename P, typename S>
class GreedyNorm : public Strategy<O, P, S>
{
private:
  // Calculates p-norm of vector x
  static constexpr LoadFloatType norm(const std::array<LoadFloatType, O::dimension>& x)
  {
    LoadFloatType result = 0;
    for (const auto& element : x)
    {
      const auto elementSq = element * element;
      result += elementSq * elementSq;
    }
    return result;
  }

public:
  GreedyNorm() = default;
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    for (const auto& o : objs)
    {
      auto minNorm = std::numeric_limits<float>::max();
      P* minProc = nullptr;
      for (auto& p : procs)
      {
        auto temp = o.load;
        for (int i = 0; i < O::dimension; i++)
        {
          temp[i] += p.load[i];
        }
        const auto tempNorm = norm(temp);
        if (tempNorm < minNorm)
        {
          minNorm = tempNorm;
          minProc = &p;
        }
      }
      solution.assign(&o, minProc);  // update solution (assumes solution updates processor load)
    }
  }
};

template <typename P, typename S>
class GreedyNorm<Obj<1>, P, S> : public Strategy<Obj<1>, P, S>
{
private:
  Greedy<Obj<1>, P, S> greedy;
public:
  GreedyNorm() = default;

  void solve(std::vector<Obj<1>>& objs, std::vector<P>& procs, S& solution,
             bool objsSorted)
  {
    greedy.solve(objs, procs, solution, objsSorted);
  }
};

template <typename O, typename P, typename S, size_t NUM_SAMPLES = 25>
class GreedySample : public Strategy<O, P, S>
{
private:
  // Calculates p-norm of vector x
  static constexpr LoadFloatType norm(const std::array<LoadFloatType, O::dimension>& x)
  {
    LoadFloatType result = 0;
    for (const auto& element : x)
    {
      const auto elementSq = element * element;
      result += elementSq * elementSq;
    }
    return result;
  }

public:
  GreedySample() = default;
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    if (procs.size() <= 256)
    {
      GreedyNorm<O, P, S> gn;
      gn.solve(objs, procs, solution, objsSorted);
      return;
    }
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    std::default_random_engine rng;

    const int numSamples = std::max(NUM_SAMPLES, (size_t)std::log2(procs.size()));
    for (const auto& o : objs)
    {
      auto minNorm = std::numeric_limits<float>::max();
      P* minProc = nullptr;
      for (size_t j = 0; j < numSamples; j++)
      {
	auto& p = procs[std::uniform_int_distribution<size_t>(0, procs.size() - 1)(rng)];
        auto temp = o.load;
        for (int i = 0; i < O::dimension; i++)
        {
          temp[i] += p.load[i];
        }
        const auto tempNorm = norm(temp);
        if (tempNorm < minNorm)
        {
          minNorm = tempNorm;
          minProc = &p;
        }
      }
      solution.assign(&o, minProc);  // update solution (assumes solution updates processor load)
    }
  }
};


template <typename O, typename P, typename S, size_t NUM_SAMPLES = 25>
class RandomScore : public Strategy<O, P, S>
{
private:
  // Calculates p-norm of vector x
  static constexpr LoadFloatType norm(const std::array<LoadFloatType, O::dimension>& x)
  {
    LoadFloatType result = 0;
    for (const auto& element : x)
    {
      const auto elementSq = element * element;
      result += elementSq * elementSq;
    }
    return result;
  }

public:
  RandomScore() = default;
  void solve(std::vector<O>& objs, std::vector<P>& procs, S& solution, bool objsSorted)
  {
    if (procs.size() <= 256)
    {
      GreedyNorm<O, P, S> gn;
      gn.solve(objs, procs, solution, objsSorted);
      return;
    }
    // Sorts by maxload in vector
    if (!objsSorted) std::sort(objs.begin(), objs.end(), CmpLoadGreater<O>());

    std::default_random_engine rng;

    std::array<float, O::dimension> maxLoad;
    maxLoad.fill(0);

    const int numSamples = std::max(NUM_SAMPLES, (size_t)std::log2(procs.size()));
    for (const auto& o : objs)
    {
      auto maxScore = std::numeric_limits<float>::lowest();
      P* minProc = nullptr;
      for (size_t j = 0; j < numSamples; j++)
      {
	auto& p = procs[std::uniform_int_distribution<size_t>(0, procs.size() - 1)(rng)];
	float score = 0;
        for (int i = 0; i < O::dimension; i++)
        {
	  score += maxLoad[i] - (p.load[i] + o.load[i]);
        }
	if (score > maxScore)
	{
	  maxScore = score;
	  minProc = &p;
	}
      }
      solution.assign(&o, minProc);  // update solution (assumes solution updates processor load)
      for (int i = 0; i < O::dimension; i++)
      {
	maxLoad[i] = std::max(maxLoad[i], (float)minProc->load[i]);
      }
    }
  }
};

}  // namespace TreeStrategy

#endif /* GREEDY_H */
