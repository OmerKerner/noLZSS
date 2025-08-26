#pragma once
#include <vector>
#include <utility>
#include <string_view>

namespace noLZSS {

struct Factor { unsigned int start; unsigned int length; };
std::vector<Factor> factorize(std::string_view text);

}
