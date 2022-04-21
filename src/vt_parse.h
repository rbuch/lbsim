#ifndef VT_PARSE_H_
#define VT_PARSE_H_

#include "TreeStrategyBase.h"
#include "json.hpp"
#include <fstream>
#include <string>

namespace vt
{
using json = nlohmann::json;

class Obj
{
    std::vector<double> loadVec;
    int oldPe;
    int id;

    Obj(std::vector<double> loadVec, int oldPe, int id) : loadVec(loadVec), oldPe(oldPe), id(id) {}
};

void loadFile(std::string filename, std::vector<std::vector<double>>& loads,
              double& bgload)
{
  std::ifstream file(filename);
  json data;
  file >> data;

  auto phases = data["phases"].get<std::vector<json>>();
  auto tasks = phases[phases.size() - 1]["tasks"].get<std::vector<json>>();
  for (const auto& task : tasks)
  {
    if (task["entity"]["migratable"].get<bool>())
    {
      auto subphases = task["subphases"].get<std::vector<json>>();
      std::vector<double> load;
      load.push_back(task["time"].get<double>());
      for (const auto& subphase : subphases)
      {
        load.push_back(subphase["time"].get<double>());
      }
      loads.push_back(load);
    }
    else
    {
      bgload += task["time"].get<double>();
    }
  }
}
};      // namespace vt
#endif  // VT_PARSE_H_
