#ifndef KD_H
#define KD_H

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>
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
protected:
  Elem data;
  TreeType* left = nullptr;
  TreeType* right = nullptr;
  KDFloatType norm;
  unsigned int size = 1;

  BaseKDNode(const Elem& key, const int numConstraints = 0) : data(key), norm(calcNorm(data, numConstraints)) {}

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

  virtual int getSplitDim(int depth = 0) = 0;

#ifdef DEBUG
public:
  static void printTree(TreeType* tree, int depth = 0, std::string prefix = "")
  {
    std::cout << prefix;
    std::cout << tree->data.id << " (" << tree->size << ", " << tree->getSplitDim(depth) << "): ";
    for (int i = 0; i < N; i++) std::cout << tree->data[i] << " ";
    std::cout << std::endl;
    if (tree->left != nullptr) printTree(tree->left, depth + 1, prefix + "L  ");
    if (tree->right != nullptr) printTree(tree->right, depth + 1, prefix + "R  ");
  }

  static size_t countNodes(TreeType* tree)
  {
    size_t count = 0;
    if (tree != nullptr) count++;
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
  RKDNode(const Elem& key, const int numConstraints) : base(key, numConstraints), discr(random(0, N - 1)) {}

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
  static Elem* findMinNorm(rkdt t, const T& x)
  {
    std::array<KDFloatType, N> mins = {0};
    KDFloatType bestNorm = std::numeric_limits<KDFloatType>::max();
    return findMinNormHelper(t, x, nullptr, bestNorm, mins);
  }

  template <typename T>
  static Elem* findMinNormNoRecurse(rkdt t, const T& x)
  {
    std::array<KDFloatType, N> mins = {0};
    KDFloatType bestNorm = std::numeric_limits<KDFloatType>::max();
    return findMinNormHelperNoRecurse(t, x, nullptr, bestNorm, mins);
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
  static Elem* findMinNormHelperNoRecurse(rkdt t, const T& x, Elem* bestObj,
                                          KDFloatType& bestNorm,
                                          std::array<KDFloatType, N>& minBounds)
  {
    using entry_t = std::pair<rkdt, std::array<KDFloatType, N>>;
    std::stack<entry_t, std::vector<entry_t>> stack;

    rkdt current = t;
    std::array<KDFloatType, N> currentBounds = minBounds;

    while (current != nullptr || !stack.empty())
    {
      while (current != nullptr)
      {
        stack.emplace(current, currentBounds);
        current = current->left;
      }

      if (!stack.empty())
      {
        auto entry = stack.top();
        rkdt node = entry.first;
        currentBounds = std::move(entry.second);
        stack.pop();

        if (node->norm < bestNorm)
        {
          const auto nodeNorm = base::calcNorm(x, node->data);
          if (nodeNorm < bestNorm)
          {
            bestObj = &(node->data);
            bestNorm = nodeNorm;
          }
        }

        if (node->right != nullptr)
        {
          const auto dim = node->discr;
          currentBounds[dim] = node->data[dim];

          // If this is true, then search the right branch
          if (base::calcNorm(x, currentBounds) < bestNorm)
          {
            current = node->right;
          }
        }
      }
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
