#ifndef LOAD_CONF_H_
#define LOAD_CONF_H_

#include <algorithm>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "json.hpp"

namespace conf
{
std::vector<double> configureHelper(const nlohmann::json& configJson,
                                           std::minstd_rand& eng, const int objIndex,
                                           int first /* = 0*/, int size /* = numChares */)
{
  std::vector<double> results;
  for (const auto& entry : configJson)
  {
    for (const auto& element : entry.items())
    {
      const auto key = element.key();
      const auto values = element.value();
      if (key == "constant")
      {
        const auto value = values["value"].get<double>();
        results.push_back(value);
      }
      else if (key == "linear")
      {
        const auto base = values["base"].get<double>();
        const auto increment = values["increment"].get<double>();
        const auto shift = values.contains("shift") ? values["shift"].get<int>() : 0;
        results.push_back(base + increment * ((size + objIndex - shift) % size));
      }
      else if (key == "normal")
      {
        const auto mean = values["mean"].get<double>();
        const auto stddev = values["stddev"].get<double>();
        std::normal_distribution<> dist(mean, stddev);
        const auto value = std::max(0.0, dist(eng));
        results.push_back(value);
      }
      else if (key == "exponential")
      {
        const auto lambda = values["lambda"].get<double>();
        std::exponential_distribution<> dist(lambda);
        results.push_back(dist(eng));
      }
      else if (key == "nestedProb")
      {
        auto ratio = values["ratio"].get<std::vector<double>>();
        std::partial_sum(ratio.begin(), ratio.end(), ratio.begin());

        std::vector<double> nestedResults =
            configureHelper(values["dists"], eng, objIndex, first, size);
        std::uniform_real_distribution<> dist(0, ratio.back());
        const auto sample = dist(eng);
        const auto index =
            std::distance(ratio.begin(), std::find_if(ratio.begin(), ratio.end(),
                                                      [=](double value)
                                                      { return sample <= value; }));
        results.push_back(nestedResults[index]);
      }
      else if (key == "nestedBlock")
      {
        auto ratio = values["ratio"].get<std::vector<double>>();
        std::partial_sum(ratio.begin(), ratio.end(), ratio.begin());
        std::transform(ratio.begin(), ratio.end(), ratio.begin(),
                       [=](double value) { return value / ratio.back() * size; });
        // At this point, ratio contains the exclusive upper bound of each block in order,
        // e.g. an original ratio of [2,1] and size of 60 would result in ratio=[40, 60]
        const auto index =
            std::distance(ratio.begin(), std::find_if(ratio.begin(), ratio.end(),
                                                      [=](double value)
                                                      { return objIndex < value; }));
        const auto blockStart = (index > 0) ? ratio[index - 1] : 0;
        const auto blockSize = ratio[index] - ((index > 0) ? ratio[index - 1] : 0);
        std::vector<double> nestedResults =
            configureHelper(nlohmann::json::array({values["dists"][index]}), eng,
                            objIndex, blockStart, blockSize);
        results.push_back(nestedResults[0]);
      }
      else
      {
        printf("No such config key \"%s\"", key.c_str());
        abort();
      }
    }
  }
  return results;
}

void loadFile(std::string filename, std::vector<std::vector<double>>& loads, const int size,
              const int seed)
{
  std::ifstream file(filename);
  nlohmann::json configJson;
  file >> configJson;

  for (int i = 0; i < size; i++)
  {
    std::minstd_rand eng(i + seed);
    auto results = configureHelper(configJson["phases"], eng, i, 0, size);
    results.insert(results.begin(), std::accumulate(results.begin(), results.end(),
                                                    decltype(results)::value_type(0)));
    loads.push_back(results);
  }
}
};  // namespace conf

#endif
