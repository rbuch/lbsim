#include "vt_parse.h"
#include <iostream>

int main(int argc, char* argv[])
{
    std::cout << argc << std::endl;

    for (int i = 1; i < argc; i++)
    {
        auto loads = parseFile(argv[i]);
        std::cout << "\n\n" << argv[i] << ":\n";
        for (const auto& loadVector : loads)
        {
            for (const auto& load : loadVector)
                std::cout << load << " ";
            std::cout << std::endl;
        }
    }

    return 0;
}
