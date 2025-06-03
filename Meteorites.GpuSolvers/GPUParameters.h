#pragma once

constexpr unsigned int ITERS_PER_KERNEL = 200;
constexpr unsigned int CASES_PER_THREAD = 1;
constexpr unsigned int THREADS_PER_BLOCK = 32;
constexpr unsigned int BLOCKS_PER_SM = 1; // 2 is the most effective (theoretically)
constexpr unsigned int CASE_NUM = 1792;   // update for actual GPU