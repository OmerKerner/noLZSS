#pragma once
#include <vector>
#include <utility>
#include <string_view>

namespace noLZSS {

// Factor: start offset and length (use size_t for portability on large inputs)
struct Factor { size_t start; size_t length; };
std::vector<Factor> factorize(std::string_view text);

}
