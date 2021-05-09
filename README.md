<!--
 * @Organization: SUSTech
 * @Author: nanoseeds
 * @Date: 2020-05-07 10:59:03
 * @LastEditTime: 2021-05-09 20:53:22
 * @License: CC-BY-NC-SA_V4_0 or any later version 
 -->
## CS205_C_CPP_2020S_project : matrix  

details look at `./material/project-matrix.md`  

1. All code(in this project,include *.c,*.cpp,*.h,*.hpp,CMakeLists.txt,etc) based on AGPL3.0(or any later version).

2. All *.md files are based on CC-BY-NC-SA-4.0(or any later version).

<!-- 3. 如果将此project作为子module使用.
  + 确保上层文件夹所在目录下存在
    + `catch.hpp` : Catch2 的头文件
    + `catch_main.hpp` : 
      只需要两行
      ``` cpp
      #define CATCH_CONFIG_MAIN  
      #include  "./catch.hpp"
      ```
  + 还不清楚可以参考 https://github.com/Certseeds/CS205_C_CPP 的结构
  + 然后添加到`CS205_project_2020S/src`到`Cmakelist.txt`中即可. -->

3. 由于opencv在github actions构建时间太离谱,默认不需要也没测试有opencv的情况下的构建与测试

4. 现在只在wsl,windows,ubuntu环境下进行了测试,没测macos下的情况.

5. 构建过程参考[github actions workflow](./.github/workflows/compile.yaml)内的命令行

6. 
  + 单元测试: test_matrix_1.cpp
  + 整体测试: test_matrix_2.cpp
  + O0,O1,O2,O3下均无异常.

6. group-member
  + [huyuhao](https://github.com/huyuhao412)
    + inverse
    + reshape and slice.  
  + [Wjia](https://github.com/Wjia0628)
    + eigenvalue
    + eigenvectors
    + max,min,sum,avg and their row/col's version.
  + me: other codes, code review, etc.

![C++](https://img.shields.io/badge/C%2B%2B-17-orange)

[![Catch2](https://img.shields.io/badge/Catch2-2.12.2-orange)][Catch2_2.12.2]

[![OpenCV](https://img.shields.io/badge/OpenCV-3.4.10-orange)][OpenCV_3.4.10]

[![AGPL3.0 Licence](https://img.shields.io/badge/License-AGPL_V3-orange)][agpl_3_0]

[![AGPL_V3](https://www.gnu.org/graphics/agplv3-with-text-162x68.png)][agpl_3_0]

[![CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-orange)][cc_by_nc_sa_4_0]

[![CC BY-SA 4.0][cc_by_nc_sa_4_0_image]][cc_by_nc_sa_4_0]

[cc_by_nc_sa_4_0]: https://creativecommons.org/licenses/by-nc-sa/4.0/

[cc_by_nc_sa_4_0_image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png

[agpl_3_0]: https://opensource.org/licenses/AGPL-3.0

[Catch2_2.12.2]: https://github.com/catchorg/Catch2/releases/tag/v2.12.2

[OpenCV_3.4.10]: https://github.com/opencv/opencv/releases/tag/3.4.10