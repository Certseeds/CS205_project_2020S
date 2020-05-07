/*  CS205_C_CPP 
    Copyright (C) 2020  nanoseeds

    CS205_C_CPP is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    CS205_C_CPP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    */
/**
 * @Github: https://github.com/Certseeds/CS205_C_CPP
 * @Organization: SUSTech
 * @Author: nanoseeds
 * @Date: 2020-05-07 11:30:27 
 * @LastEditors  : nanoseeds
 */
#include "./../../catch.hpp"
#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using Catch::Matchers::Equals;
auto function = [](const std::vector<int32_t> &t1,
                   const std::vector<int32_t> &t2) -> bool {
    assert(t1.size() == t2.size());
    for (uint32_t i = 0; i < t1.size(); i++) {
        if (t1[i] > t2[i]) {
            return true;
        } else if (t1[i] < t2[i]) {
            return false;
        }
    }
    return true;
};

TEST_CASE("test 0", "[test1]") {
    REQUIRE(function({CATCH_VERSION_MAJOR, CATCH_VERSION_MINOR, CATCH_VERSION_PATCH}, {2, 12, 1}));
    REQUIRE(function({2, 13, 0}, {2, 12, 1}));
    REQUIRE(function({3, 0, 12}, {2, 12, 1}));
}

TEST_CASE("test 1", "[test1]") {
    CHECK(1 * 2 == 1 + 1);
}