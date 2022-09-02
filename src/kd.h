#ifndef KD_H
#define KD_H

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
#define DEBUG
#ifdef DEBUG
#  include <iostream>
#endif
#include <limits>
#include <random>
#include <stack>
#include <utility>

using KDFloatType = double;

template <typename TreeType, typename Elem, int Exp, int N = Elem::dimension>
class BaseKDNode
{
private:
  template <size_t chunksPerBlock>
  class PoolAllocator
  {
  private:
    struct Chunk
    {
      Chunk* next;
    };

    Chunk* currentChunk = nullptr;
    std::vector<Chunk*> blocks;

    Chunk* allocateBlock(size_t chunkSize)
    {
      const size_t blockSize = chunksPerBlock * chunkSize;
      Chunk* blockBegin = reinterpret_cast<Chunk*>(malloc(blockSize));
      Chunk* chunk = blockBegin;

      for (int i = 0; i < chunksPerBlock - 1; ++i)
      {
        chunk->next =
            reinterpret_cast<Chunk*>(reinterpret_cast<char*>(chunk) + chunkSize);
        chunk = chunk->next;
      }

      chunk->next = nullptr;

      return blockBegin;
    }

  public:
    PoolAllocator() {}

    void* allocate(size_t size)
    {
      if (currentChunk == nullptr)
      {
        currentChunk = allocateBlock(size);
        blocks.push_back(currentChunk);
      }

      Chunk* freeChunk = currentChunk;
      currentChunk = currentChunk->next;
      return freeChunk;
    }

    void deallocate(void* chunk, size_t size)
    {
      reinterpret_cast<Chunk*>(chunk)->next = currentChunk;
      currentChunk = reinterpret_cast<Chunk*>(chunk);
    }

    void reset()
    {
      for (Chunk* block : blocks) free(block);
      blocks.clear();
      currentChunk = nullptr;
    }
  };

public:
  static void* operator new(std::size_t size) { return alloc.allocate(size); }
  static void operator delete(void* p, std::size_t size)
  {
    alloc.deallocate((TreeType*)p, size);
  }

  static void freeAll() { alloc.reset(); }

protected:
  Elem data;
  TreeType* left = nullptr;
  TreeType* right = nullptr;
  KDFloatType norm;
  unsigned int size = 1;

  static inline PoolAllocator<4096> alloc;

  BaseKDNode(const Elem& key, const int numConstraints = 0) : data(key), norm(calcNorm(data, numConstraints)) {}
  ~BaseKDNode()
  {
    if (left != nullptr) delete left;
    if (right != nullptr) delete right;
  }

  template <typename A>
  static KDFloatType calcNorm(const A& x)
  {
    KDFloatType norm = 0;
    if constexpr (Exp >= 32)
    {
      norm = *std::max_element(x.begin(), x.end());
    }
    else
    {
      for (int i = 0; i < N; i++)
      {
        norm += std::pow(x[i], Exp);
      }
    }
    return norm;
  }

  template <typename A, typename B,
            typename std::enable_if<!std::is_integral<B>::value, bool>::type = true>
  static KDFloatType calcNorm(const A& a, const B& b)
  {
    KDFloatType norm = 0;
    if constexpr (Exp >= 32)
    {
      std::array<KDFloatType, N> sum;
      for (int i = 0; i < N; i++)
      {
        sum[i] = a[i] + b[i];
      }
      norm = *std::max_element(sum.begin(), sum.end());
    }
    else
    {
      for (int i = 0; i < N; i++)
      {
        norm += std::pow(a[i] + b[i], Exp);
      }
    }
    return norm;
  }

  template <typename A, typename B = size_t,
            typename std::enable_if<std::is_integral<B>::value, bool>::type = true>
  static KDFloatType calcNorm(const A& x, const B numConstraints)
  {
    KDFloatType norm = 0;
    if constexpr (Exp >= 32)
    {
      norm = *std::max_element(x.begin(), x.end() - numConstraints);
    }
    else
    {
      for (int i = 0; i < N - numConstraints; i++)
      {
        norm += std::pow(x[i], Exp);
      }
    }
    return norm;
  }

  template <typename A, typename B, typename C = size_t,
            typename std::enable_if<!std::is_integral<B>::value, bool>::type = true>
  static KDFloatType calcNorm(const A& a, const B& b, const C numConstraints)
  {
    KDFloatType norm = 0;
    if constexpr (Exp >= 32)
    {
      std::array<KDFloatType, N> sum;
      for (int i = 0; i < N; i++)
      {
        sum[i] = a[i] + b[i];
      }
      norm = *std::max_element(sum.begin(), sum.end() - numConstraints);
    }
    else
    {
      for (int i = 0; i < N - numConstraints; i++)
      {
        norm += std::pow(a[i] + b[i], Exp);
      }
    }
    return norm;
  }

  template <typename A, typename B>
  static std::array<KDFloatType, N> addVecs(const A& a, const B& b)
  {
    std::array<KDFloatType, N> sum = {0};
    for (int i = 0; i < N; i++)
    {
      sum[i] = a[i] + b[i];
    }

    return sum;
  }

  template <typename A, typename B>
  static KDFloatType distance2(const A& a, const B& b)
  {
    KDFloatType distance2 = 0;
    for (int i = 0; i < N; i++)
    {
      distance2 += std::pow(a[i] - b[i], 2);
    }

    return distance2;
  }

  // true if a dominates b
  template <typename A, typename B>
  static bool isDominated(const A& a, const B& b)
  {
    for (int i = 0; i < N; i++)
    {
      if (b[i] < a[i])
        return false;
    }
    return true;
  }

  virtual int getSplitDim(int depth = 0) = 0;

#ifdef DEBUG
public:
  static void printTree(TreeType* tree, int depth = 0, std::string prefix = "")
  {
    if (tree == nullptr) return;
    std::cout << prefix;
    std::cout << tree->data.id << " (" << tree->size << ", " << tree->getSplitDim(depth)
              << "): ";
    for (int i = 0; i < N; i++) std::cout << tree->data[i] << " ";
    std::cout << std::endl;
    if (tree->left != nullptr) printTree(tree->left, depth + 1, prefix + "L  ");
    if (tree->right != nullptr) printTree(tree->right, depth + 1, prefix + "R  ");
  }

  static size_t countNodes(TreeType* tree)
  {
    size_t count = 0;
    if (tree == nullptr)
      return count;
    else
      count++;
    if (tree->left != nullptr) count += countNodes(tree->left);
    if (tree->right != nullptr) count += countNodes(tree->right);

    return count;
  }

  static bool checkTree(TreeType* tree)
  {
    std::vector<std::pair<TreeType*, bool>> accum;
    return checkTreeHelper(tree, accum);
  }

private:
  static bool checkTreeHelper(TreeType* tree,
                              std::vector<std::pair<TreeType*, bool>>& accum)
  {
    bool valid = true;
    for (int depth = 0; depth < accum.size(); depth++)
    {
      TreeType* node = accum[depth].first;
      bool isLeft = accum[depth].second;
      int splitDim = node->getSplitDim(depth);

      valid &= (isLeft) ? tree->data[splitDim] < node->data[splitDim]
        : tree->data[splitDim] >= node->data[splitDim];

      if (!valid) return false;
    }

    if (tree->left != nullptr && valid)
    {
      accum.emplace_back(tree, true);
      valid &= checkTreeHelper(tree->left, accum);
      accum.pop_back();
    }

    if (tree->right != nullptr && valid)
    {
      accum.emplace_back(tree, false);
      valid &= checkTreeHelper(tree->right, accum);
      accum.pop_back();
    }

    return valid;
  }

#endif
};

template <typename Elem, int Exp = 4, int N = Elem::dimension>
class KDNode : public BaseKDNode<KDNode<Elem, Exp, N>, Elem, Exp>
{
  using base = BaseKDNode<KDNode<Elem, Exp, N>, Elem, Exp>;
  using kdt = KDNode*;
public:
  KDNode(const Elem& key, const int numConstraints = 0) : base(key, numConstraints) {}

  int getSplitDim(int depth = 0) override { return depth % N; }

  static kdt insert(kdt t, const Elem& x, const int numConstraints = 0, unsigned int depth = 0)
  {
    if (t == nullptr)
      t = new KDNode(x, numConstraints);
    else
    {
      t->size++;
      const auto dim = depth % N;
      if (x[dim] < t->data[dim])
	t->left = insert(t->left, x, numConstraints, depth + 1);
      else
	t->right = insert(t->right, x, numConstraints, depth + 1);
    }
    return t;
  }

  static Elem findMin(kdt t, unsigned int targetDim, unsigned int depth = 0)
  {
    assert(t != nullptr);
    const auto dim = depth % N;
    if (dim == targetDim)
    {
      if (t->left == nullptr)
	return t->data;
      else
	return findMin(t->left, targetDim, depth + 1);
    }
    else
    {
      auto obj = t->data;
      if (t->left != nullptr)
      {
	const auto left = findMin(t->left, targetDim, depth + 1);
	if (left[targetDim] < obj[targetDim])
	  obj = left;
      }
      if (t->right != nullptr)
      {
	const auto right = findMin(t->right, targetDim, depth + 1);
	if (right[targetDim] < obj[targetDim])
	  obj = right;
      }
      return obj;
    }
  }

  static kdt remove(kdt t, const Elem& x, unsigned int depth = 0)
  {
    const auto dim = depth % N;
    t->size--;

    if (t->data == x)
    {
      if (t->right != nullptr)
      {
	t->data = findMin(t->right, dim, depth + 1);
	t->right = remove(t->right, t->data, depth + 1);
      }
      else if (t->left != nullptr)
      {
	t->data = findMin(t->left, dim, depth + 1);
	t->right = remove(t->left, t->data, depth + 1);
	t->left = nullptr;
      }
      else
      {
	delete t;
	t = nullptr;
      }
    }
    else if (x[dim] < t->data[dim])
      t->left = remove(t->left, x, depth + 1);
    else
      t->right = remove(t->right, x, depth + 1);

    return t;
  }

  template <typename T>
  static Elem* findMinNorm(kdt t, const T& x)
  {
    std::array<KDFloatType, N> mins = {0};
    KDFloatType bestNorm = std::numeric_limits<KDFloatType>::max();
    return findMinNormHelper(t, x, 0, nullptr, bestNorm, mins);
  }

  template <typename T>
  static Elem* findMinNormConstraints(kdt t, const T& x, const std::vector<KDFloatType>& constraints)
  {
    std::array<KDFloatType, N> mins = {0};
    KDFloatType bestNorm = std::numeric_limits<KDFloatType>::max();
    Elem* result = findMinNormHelperC(t, x, 0, nullptr, bestNorm, mins, constraints);
    CkAssertMsg(result != nullptr, "Cannot satisfy constraint!");
    return result;
  }

private:
  template <typename T>
  static Elem* findMinNormHelper(kdt t, const T& x, unsigned int depth, Elem* bestObj,
                                 KDFloatType& bestNorm, std::array<KDFloatType, N>& minBounds)
  {
    const auto dim = depth % N;

    if (t->left != nullptr)
    {
      bestObj = findMinNormHelper(t->left, x, depth + 1, bestObj, bestNorm, minBounds);
    }
    if (t->norm < bestNorm)
    {
      const auto rootNorm = base::calcNorm(x, t->data);
      if (rootNorm < bestNorm)
      {
	bestObj = &(t->data);
	bestNorm = rootNorm;
      }
    }
    if (t->right != nullptr)
    {
      const auto oldMin = minBounds[dim];
      minBounds[dim] = t->data[dim];
      if (base::calcNorm(x, minBounds) < bestNorm)
      {
	bestObj = findMinNormHelper(t->right, x, depth + 1, bestObj, bestNorm, minBounds);
      }
      minBounds[dim] = oldMin;
    }

    return bestObj;
  }

  template <typename T>
  static Elem* findMinNormHelperC(kdt t, const T& x, unsigned int depth,
                                            Elem* bestObj, KDFloatType& bestNorm,
                                            std::array<KDFloatType, N>& minBounds,
                                            const std::vector<KDFloatType>& constraints)
  {
    const auto dim = depth % N;
    const int numConstraints = constraints.size();

    if (t->left != nullptr)
    {
      bestObj = findMinNormHelperC(t->left, x, depth + 1, bestObj, bestNorm, minBounds, constraints);
    }
    if (t->norm < bestNorm)
    {
      const auto sum = base::addVecs(x, t->data);
      const auto rootNorm = base::calcNorm(sum, numConstraints);
      const auto constraintsMet =
          std::equal(constraints.begin(), constraints.end(), sum.end() - numConstraints,
                     [&](const KDFloatType& a, const KDFloatType& b) { return a >= b; });
      if (rootNorm < bestNorm && constraintsMet)
      {
        bestObj = &(t->data);
        bestNorm = rootNorm;
      }
    }
    if (t->right != nullptr)
    {
      const auto oldMin = minBounds[dim];
      minBounds[dim] = t->data[dim];
      const auto constraintMet =
          (dim >= N - numConstraints)
              ? minBounds[dim] + x[dim] < constraints[dim - (N - numConstraints)]
              : true;
      if (constraintMet && base::calcNorm(x, minBounds, numConstraints) < bestNorm)
      {
        bestObj = findMinNormHelperC(t->right, x, depth + 1, bestObj, bestNorm, minBounds, constraints);
      }
      minBounds[dim] = oldMin;
    }

    return bestObj;
  }
};

template <typename Elem, int Exp = 4, int N = Elem::dimension,
          typename Engine = std::default_random_engine>
class RKDNode : public BaseKDNode<RKDNode<Elem, Exp, N, Engine>, Elem, Exp>
{
  using base = BaseKDNode<RKDNode<Elem, Exp, N, Engine>, Elem, Exp>;
  using rkdt = RKDNode*;
private:
  int discr;

  static int random(int min, int max)
  {
    static Engine rng;
    return std::uniform_int_distribution<int>(min, max)(rng);
  }

public:
  RKDNode(const Elem& key, const int numConstraints = 0) : base(key, numConstraints), discr(random(0, N - 1)) {}

  int getSplitDim(int depth = 0) override { return discr; }

  static rkdt insert(rkdt t, const Elem& x, const int numConstraints = 0)
  {
    std::stack<std::pair<rkdt, bool>> stack;
    while (t != nullptr && random(0, t->size) > 0)
    {
      t->size++;
      const int i = t->discr;
      const bool isLeftChild = x[i] < t->data[i];
      stack.emplace(t, isLeftChild);
      if (isLeftChild)
      {
        t = t->left;
      }
      else
      {
        t = t->right;
      }
    }

    t = insert_at_root(t, x, numConstraints);

    while (!stack.empty())
    {
      const std::pair<rkdt, bool>& entry = stack.top();
      const rkdt parent = entry.first;
      const bool isLeftChild = entry.second;
      if (isLeftChild)
        parent->left = t;
      else
        parent->right = t;
      t = parent;
      stack.pop();
    }
    return t;
  }

  static rkdt remove(rkdt t, const Elem& x)
  {
    if (t == nullptr) return nullptr;
    const int i = t->discr;
    if (t->data == x)
    {
      const auto newRoot = join(t->left, t->right, i);
      t->left = nullptr;
      t->right = nullptr;
      delete t;
      return newRoot;
    }
    t->size--;
    if (x[i] < t->data[i])
      t->left = remove(t->left, x);
    else
      t->right = remove(t->right, x);
    return t;
  }

  static rkdt insert_at_root(rkdt t, const Elem& x, const int numConstraints = 0)
  {
    rkdt r = new RKDNode(x, numConstraints);
    if (t != nullptr) r->size += t->size;
    auto p = split(t, r);
    r->left = p.first;
    r->right = p.second;
    return r;
  }

  // This splits t's subtrees at whatever r specifies
  // t might be nullptr, r cannot be
  static std::pair<rkdt, rkdt> split(rkdt t, rkdt r)
  {
    if (t == nullptr) return std::make_pair(nullptr, nullptr);
    const int i = r->discr;
    const int j = t->discr;

    // t and r discriminate in the same dimension, so just resplit the appropriate subtree
    // of t
    if (i == j)
    {
      if (t->data[i] < r->data[i])
      {
        const auto p = split(t->right, r);
        t->right = p.first;
        const int splitSize = (p.second == nullptr) ? 0 : p.second->size;
        t->size -= splitSize;
        return std::make_pair(t, p.second);
      }
      else
      {
        const auto p = split(t->left, r);
        t->left = p.second;
        const int splitSize = (p.first == nullptr) ? 0 : p.first->size;
        t->size -= splitSize;
        return std::make_pair(p.first, t);
      }
    }
    // t and r discriminate on different dimensions, so recursively split both subtrees of
    // t
    else
    {
      const auto L = split(t->left, r);
      const auto R = split(t->right, r);

      if (t->data[i] < r->data[i])
      {
        t->left = L.first;
        t->right = R.first;
        const int splitSize = ((L.second == nullptr) ? 0 : L.second->size) +
                              ((R.second == nullptr) ? 0 : R.second->size);
        t->size -= splitSize;

        return std::make_pair(t, join(L.second, R.second, j));
      }
      else
      {
        t->left = L.second;
        t->right = R.second;
        const int splitSize = ((L.first == nullptr) ? 0 : L.first->size) +
                              ((R.first == nullptr) ? 0 : R.first->size);
        t->size -= splitSize;

        return std::make_pair(join(L.first, R.first, j), t);
      }
    }
  }

  static rkdt join(rkdt l, rkdt r, int dim)
  {
    if (l == nullptr) return r;
    if (r == nullptr) return l;

    const int m = l->size;
    const int n = r->size;
    const int u = random(0, m + n - 1);
    if (u < m)
    {
      l->size += r->size;
      if (l->discr == dim)
      {
        l->right = join(l->right, r, dim);
        return l;
      }
      else
      {
        auto R = split(r, l);
        l->left = join(l->left, R.first, dim);
        l->right = join(l->right, R.second, dim);
        return l;
      }
    }
    else
    {
      r->size += l->size;
      if (r->discr == dim)
      {
        r->left = join(l, r->left, dim);
        return r;
      }
      else
      {
        auto L = split(l, r);
        r->left = join(L.first, r->left, dim);
        r->right = join(L.second, r->right, dim);
        return r;
      }
    }
  }

  template <typename T>
  static std::array<KDFloatType, N> getNN(rkdt t, const T& x)
  {
    KDFloatType distance2 = std::numeric_limits<KDFloatType>::max();
    const auto nn = getNNHelper(t, x, distance2);
    std::array<KDFloatType, N> result = {0};
    for (int i = 0; i < N; i++)
    {
      result[i] = nn->data[i];
    }

    return result;
  }

  template <typename T>
  static rkdt getNNHelper(rkdt t, const T& x, KDFloatType& bestDistance2)
  {
    const auto discr = t->discr;
    const bool leftFirst = x[discr] < t->data[discr];
    const auto first = (leftFirst) ? t->left : t->right;
    const auto second = (leftFirst) ? t->right : t->left;
    rkdt best = (first == nullptr) ? t : getNNHelper(first, x, bestDistance2);

    const auto tDistance2 = base::distance2(x, t->data);
    if (tDistance2 < bestDistance2)
    {
      best = t;
      bestDistance2 = tDistance2;
    }

    if (second != nullptr && std::pow(x[discr] - t->data[discr], 2) < bestDistance2)
    {
      auto candidateDistance2 = bestDistance2;
      const auto candidate = getNNHelper(second, x, candidateDistance2);
      if (candidateDistance2 < bestDistance2)
      {
        best = candidate;
        bestDistance2 = candidateDistance2;
      }
    }

    return best;
  }

  static rkdt getParetoFrontier(rkdt t)
  {
    std::array<KDFloatType, N> minBounds = {0};
    return getParetoFrontierHelper(t, minBounds, nullptr);
  }

  static rkdt getParetoFrontier(const rkdt t, rkdt paretoFrontier = nullptr)
  {
    std::array<KDFloatType, N> minBounds = {0};
    return getParetoFrontierHelper(t, minBounds, paretoFrontier);
  }

  static rkdt updateParetoFrontier(const rkdt t, std::array<KDFloatType, N>& minBounds,
                                  rkdt paretoFrontier, const std::array<KDFloatType, N>& nn)
  {
    if (t == nullptr)
      return paretoFrontier;

    const auto dim = t->discr;
    const auto tVal = t->data[dim];
    const auto mVal = minBounds[dim];

    if (mVal < tVal && t->left != nullptr)
    {
      paretoFrontier = updateParetoFrontier(t->left, minBounds, paretoFrontier, nn);
    }

    if (mVal <= tVal && !base::isDominated(nn, t->data) && !isDominated(paretoFrontier, t->data))
    {
      paretoFrontier = insert(paretoFrontier, t->data);
    }

    if (t->right != nullptr)
    {
      if (mVal < tVal)
      {
        minBounds[dim] = tVal;
        if (!base::isDominated(nn, minBounds) && !isDominated(paretoFrontier, minBounds))
        {
          paretoFrontier = updateParetoFrontier(t->right, minBounds, paretoFrontier, nn);
        }
        minBounds[dim] = mVal;
      }
      else
        paretoFrontier = updateParetoFrontier(t->right, minBounds, paretoFrontier, nn);
    }

    return paretoFrontier;
  }

  static rkdt getParetoFrontierHelper(const rkdt t, std::array<KDFloatType, N>& minBounds,
                                      rkdt paretoFrontier)
  {
    if (t == nullptr) return paretoFrontier;

    if (t->left != nullptr)
      paretoFrontier = getParetoFrontierHelper(t->left, minBounds, paretoFrontier);

    if (!isDominated(paretoFrontier, t->data))
      paretoFrontier = insert(paretoFrontier, t->data);

    if (t->right != nullptr)
    {
      const auto dim = t->discr;
      const auto temp = minBounds[dim];
      minBounds[dim] = t->data[dim];
      if (!isDominated(paretoFrontier, minBounds))
      {
        paretoFrontier = getParetoFrontierHelper(t->right, minBounds, paretoFrontier);
      }
      minBounds[dim] = temp;
    }

    return paretoFrontier;
  }

  // true if x is dominated by some node in t (i.e. some node is <= x in every dim), else
  // false
  template <typename T>
  static bool isDominated(rkdt t, const T& x)
  {
    std::vector<rkdt> stack;
    stack.reserve(128);

    if (t != nullptr)
      stack.push_back(t);

    while (!stack.empty())
    {
      const auto current = stack.back();
      stack.pop_back();

      const auto dim = current->discr;
      if (x[dim] < current->data[dim])
      {
        if (current->left != nullptr)
          stack.push_back(current->left);
        continue;
      }

      if (base::isDominated(current->data, x))
        return true;

      if (current->left != nullptr)
        stack.push_back(current->left);
      if (current->right != nullptr)
        stack.push_back(current->right);
    }

    return false;
  }

  template <typename T>
  static Elem* findMinNorm(rkdt t, const T& x)
  {
    std::array<KDFloatType, N> mins = {0};
    KDFloatType bestNorm = std::numeric_limits<KDFloatType>::max();
    return findMinNormHelper(t, x, nullptr, bestNorm, mins);
  }

  template <typename T>
  static Elem* findMinNormObjNorm(rkdt t, const T& x)
  {
    std::array<KDFloatType, N> mins = {0};
    KDFloatType bestNorm = std::numeric_limits<KDFloatType>::max();
    return findMinNormHelperObjNorm(t, x, nullptr, bestNorm, mins, base::calcNorm(x));
  }

  template <typename T>
  static Elem* findMinNormObjNormEarly(rkdt t, const T& x, const std::array<KDFloatType, N>& maxLoads)
  {
    std::array<KDFloatType, N> mins = {0};
    KDFloatType bestNorm = std::numeric_limits<KDFloatType>::max();
    int earlyExit = 5;
    return findMinNormHelperObjNormEarly(t, x, nullptr, bestNorm, mins, base::calcNorm(x), maxLoads, earlyExit);
  }

  template <typename T>
  static Elem* findMinNormConstraints(rkdt t, const T& x,
                                      const std::vector<KDFloatType>& constraints)
  {
    std::array<KDFloatType, N> mins = {0};
    KDFloatType bestNorm = std::numeric_limits<KDFloatType>::max();
    Elem* result = findMinNormHelperC(t, x, 0, nullptr, bestNorm, mins, constraints);
    CkAssertMsg(result != nullptr, "Cannot satisfy constraint!");
    return result;
  }

private:
  template <typename T>
  static Elem* findMinNormHelper(rkdt t, const T& x, Elem* bestObj, KDFloatType& bestNorm,
                                 std::array<KDFloatType, N>& minBounds)
  {
    if (t->left != nullptr)
    {
      bestObj = findMinNormHelper(t->left, x, bestObj, bestNorm, minBounds);
    }
    if (t->norm < bestNorm)
    {
      const auto rootNorm = base::calcNorm(x, t->data);
      if (rootNorm < bestNorm)
      {
        bestObj = &(t->data);
        bestNorm = rootNorm;
      }
    }
    if (t->right != nullptr)
    {
      const auto dim = t->discr;
      const auto oldMin = minBounds[dim];
      minBounds[dim] = t->data[dim];
      if (base::calcNorm(x, minBounds) < bestNorm)
      {
        bestObj = findMinNormHelper(t->right, x, bestObj, bestNorm, minBounds);
      }
      minBounds[dim] = oldMin;
    }

    return bestObj;
  }

  template <typename T>
  static Elem* findMinNormHelperObjNorm(rkdt t, const T& x, Elem* bestObj, KDFloatType& bestNorm,
                                 std::array<KDFloatType, N>& minBounds, const KDFloatType xNorm)
  {
    if (t->left != nullptr)
    {
      bestObj = findMinNormHelperObjNorm(t->left, x, bestObj, bestNorm, minBounds, xNorm);
    }
    if (t->norm + xNorm < bestNorm)
    {
      const auto rootNorm = base::calcNorm(x, t->data);
      if (rootNorm < bestNorm)
      {
        bestObj = &(t->data);
        bestNorm = rootNorm;
      }
    }
    if (t->right != nullptr)
    {
      const auto dim = t->discr;
      const auto oldMin = minBounds[dim];
      minBounds[dim] = t->data[dim];
      if (base::calcNorm(x, minBounds) < bestNorm)
      {
        bestObj = findMinNormHelperObjNorm(t->right, x, bestObj, bestNorm, minBounds, xNorm);
      }
      minBounds[dim] = oldMin;
    }

    return bestObj;
  }

  template <typename T>
  static Elem* findMinNormHelperObjNormEarly(rkdt t, const T& x, Elem* bestObj,
                                             KDFloatType& bestNorm,
                                             std::array<KDFloatType, N>& minBounds,
                                             const KDFloatType xNorm,
                                             const std::array<KDFloatType, N>& maxLoads,
                                             int& earlyExit)
  {
    if (t->left != nullptr)
    {
      bestObj = findMinNormHelperObjNormEarly(t->left, x, bestObj, bestNorm, minBounds, xNorm, maxLoads, earlyExit);
    }
    if (earlyExit > 0 && t->norm + xNorm < bestNorm)
    {
      const auto rootNorm = base::calcNorm(x, t->data);
      if (rootNorm < bestNorm)
      {
        bestObj = &(t->data);
        bestNorm = rootNorm;

        bool found = true;
        for (int i = 0; i < N; i++)
        {
          if (t->data[i] + x[i] > maxLoads[i])
          {
            found = false;
            break;
          }
        }
        if (found)
        {
          earlyExit--;
        }
      }
    }
    if (earlyExit > 0 && t->right != nullptr)
    {
      const auto dim = t->discr;
      const auto oldMin = minBounds[dim];
      minBounds[dim] = t->data[dim];
      if (base::calcNorm(x, minBounds) < bestNorm)
      {
        bestObj =
            findMinNormHelperObjNormEarly(t->right, x, bestObj, bestNorm, minBounds, xNorm, maxLoads, earlyExit);
      }
      minBounds[dim] = oldMin;
    }

    return bestObj;
  }

  template <typename T>
  static Elem* findMinNormHelperC(rkdt t, const T& x, unsigned int depth, Elem* bestObj,
                                  KDFloatType& bestNorm,
                                  std::array<KDFloatType, N>& minBounds,
                                  const std::vector<KDFloatType>& constraints)
  {
    const auto dim = t->discr;
    const int numConstraints = constraints.size();

    if (t->left != nullptr)
    {
      bestObj = findMinNormHelperC(t->left, x, depth + 1, bestObj, bestNorm, minBounds,
                                   constraints);
    }
    if (t->norm < bestNorm)
    {
      const auto sum = base::addVecs(x, t->data);
      const auto rootNorm = base::calcNorm(sum, numConstraints);
      const auto constraintsMet =
          std::equal(constraints.begin(), constraints.end(), sum.end() - numConstraints,
                     [&](const KDFloatType& a, const KDFloatType& b) { return a >= b; });
      if (rootNorm < bestNorm && constraintsMet)
      {
        bestObj = &(t->data);
        bestNorm = rootNorm;
      }
    }
    if (t->right != nullptr)
    {
      const auto oldMin = minBounds[dim];
      minBounds[dim] = t->data[dim];
      const auto constraintMet =
          (dim >= N - numConstraints)
              ? minBounds[dim] + x[dim] < constraints[dim - (N - numConstraints)]
              : true;
      if (constraintMet && base::calcNorm(x, minBounds, numConstraints) < bestNorm)
      {
        bestObj = findMinNormHelperC(t->right, x, depth + 1, bestObj, bestNorm, minBounds,
                                     constraints);
      }
      minBounds[dim] = oldMin;
    }

    return bestObj;
  }

  static Elem findMin(rkdt t, unsigned int targetDim)
  {
    assert(t != nullptr);
    const auto dim = t->discr;
    if (dim == targetDim)
    {
      if (t->left == nullptr)
        return t->data;
      else
        return findMin(t->left, targetDim);
    }
    else
    {
      auto obj = t->data;
      if (t->left != nullptr)
      {
        const auto left = findMin(t->left, targetDim);
        if (left[targetDim] < obj[targetDim]) obj = left;
      }
      if (t->right != nullptr)
      {
        const auto right = findMin(t->right, targetDim);
        if (right[targetDim] < obj[targetDim]) obj = right;
      }
      return obj;
    }
  }
};

#endif /* KD_H */
