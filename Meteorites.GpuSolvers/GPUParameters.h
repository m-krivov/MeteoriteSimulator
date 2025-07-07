#pragma once

constexpr unsigned int ITERS_PER_KERNEL = 200;
constexpr unsigned int CASES_PER_THREAD = 1;
// There is only examination of grid_idx in AdamsKernel. It means CASE_NUM
// must be a multiple of CASES_PER_THREAD. But we have Solve() method for one case.
constexpr unsigned int THREADS_PER_BLOCK = 32;
//constexpr unsigned int BLOCKS_PER_SM = 2; // 2 is the most effective (theoretically)
constexpr unsigned int CASE_NUM = 100000; // !!! (can't include Main.cpp)
// Ok, i made Solve() method for one case, but on stage2 calculations will be very very slow.
// (Maybe i should write method with vector of cases as parameter?) 