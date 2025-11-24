cmake_minimum_required(VERSION 3.22.1)
project(Euler_DG_Solver C)

# 设置C标准
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# 编译选项
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O2")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 -DDEBUG")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -march=native")

# 添加数学库
find_library(MATH_LIBRARY m)

# 包含目录
include_directories(${CMAKE_SOURCE_DIR}/include)

# 源文件（只包含初始化相关的，暂时不编译时间和空间离散）
set(INIT_SOURCES
    src/Gauss_quadrature.c
    src/Legendre_basis.c
    src/initial_conditions.c
    src/io_utils.c
    src/TVD_Limiter.c
    src/spatial_discretization.c
    src/time_integration.c
    src/Euler_exact.c
    
)

# ============ 可执行文件: 测试初始化 ============
add_executable(test_init
    apps/main.c
    ${INIT_SOURCES}
)
target_link_libraries(test_init ${MATH_LIBRARY})

# ============ 可执行文件: 完整求解器（等时间空间离散写好后使用）============
# add_executable(euler_solver
#     apps/main.c
#     ${INIT_SOURCES}
#     src/time_integration.c
#     src/spatial_discretization.c
# )
# target_link_libraries(euler_solver ${MATH_LIBRARY})

# ============ 输出目录设置 ============
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# 创建输出目录
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)

# ============ 安装规则 ============
install(TARGETS test_init DESTINATION bin)

# ============ 打印配置信息 ============
message(STATUS "")
message(STATUS "==================================================")
message(STATUS "  1D Euler DG Solver - Initialization Test")
message(STATUS "==================================================")
message(STATUS "  C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C Flags: ${CMAKE_C_FLAGS}")
message(STATUS "  Source directory: ${CMAKE_SOURCE_DIR}")
message(STATUS "  Binary directory: ${CMAKE_BINARY_DIR}")
message(STATUS "  Executable: test_init")
message(STATUS "==================================================")
message(STATUS "")