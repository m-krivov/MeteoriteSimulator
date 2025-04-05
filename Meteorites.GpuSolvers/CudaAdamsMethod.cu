#include "Meteorites.Core/Case.h"
#include "Meteorites.Core/Functionals.h"
#include "Meteorites.Core/ResultFormatters.h"

#include <cuda_runtime.h>

#include <vector>
#include <iostream>

constexpr unsigned int ITERS_PER_KERNEL = 10000;
constexpr unsigned int CASES_PER_THREAD = 1;

static void
HandleError(cudaError_t err, const char *file, int line)
{
    if (err != cudaSuccess) {
        printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
        exit(EXIT_FAILURE);
    }
}
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

__device__ inline double cudaExp(double x) { return exp(x); }
__device__ inline float cudaExp(float x) { return expf(x); }
__device__ inline double cudaPow(double x, double y) { return pow(x, y); }
__device__ inline float cudaPow(float x, float y) { return powf(x, y); }
__device__ inline double cudaSin(double x) { return sin(x); }
__device__ inline float cudaSin(float x) { return sinf(x); }
__device__ inline double cudaCos(double x) { return cos(x); }
__device__ inline float cudaCos(float x) { return cosf(x); }

class CudaConstants
{
  public:
    __device__ static constexpr real
    R()
    {
        return (real) 6371000;
    }

    __device__ static constexpr real
    g(real h)
    {
        return ((real) (6.67428 * 5.9726) / ((6371 + h / 1000) * (6371 + h / 1000)) * (real) 1e7);
    }

    __device__ static constexpr real
    g()
    {
        return g(0);
    }

    __device__ static real
    rho_a(real h)
    {
        constexpr real rho_0 = (real) 1.125;
        constexpr real alpha = (real) (0.029 * 9.8 / (8.31 * 273));
        return rho_0 * cudaExp(-alpha * h);
    }

    __device__ static real
    Midsection(real M, real rho)
    {
        return (real) (M_PI * pow(3.0 / (4.0 * M_PI), 2.0 / 3.0)) * pow(M / rho, (real) (2.0 / 3.0));
    }
};

class CudaCase
{
  public:
    real H_, ch_, rho_, cd_, cl_;
    real m0_, v0_, h0_, gamma0_;
};

class Layer
{
  public:
    real V_, gamma_, h_, l_, M_;
    real fV_, fgamma_, fh_, fl_, fM_;

    __device__ void
    Set(real new_V, real new_gamma, real new_h, real new_l, real new_M, const CudaCase &params)
    {
        V_ = new_V;
        gamma_ = new_gamma;
        h_ = new_h;
        l_ = new_l;
        M_ = new_M;

        if (new_M <= (real) 0.0) // probably, 'dt' is too large
        {
            fV_ = fh_ = fl_ = fM_ = (real) 0.0;
        } else {
            auto sin_gamma = cudaSin(gamma_);
            auto cos_gamma = cudaCos(gamma_);
            auto g = CudaConstants::g(h_);
            auto rho_a = CudaConstants::rho_a(h_);
            auto midsection = CudaConstants::Midsection(M_, params.rho_);

            fV_ = -params.cd_ * rho_a * V_ * V_ * midsection / (2 * M_) + g * sin_gamma;
            fgamma_ = +g * cos_gamma / V_ - V_ * cos_gamma / CudaConstants::R() -
                      params.cl_ * rho_a * V_ * midsection / (2 * M_);
            fh_ = -V_ * sin_gamma;
            fl_ = V_ * (CudaConstants::R() / (CudaConstants::R() + h_)) * cos_gamma;
            fM_ = -(params.ch_ * rho_a * V_ * V_ * V_ * midsection / 2) / params.H_;
        }
    }
};

class Record
{
  public:
    real t_, m_, v_, h_, l_, gamma_;
};

template <unsigned int STEPS, unsigned int ITERS> class ThreadContext
{
  public:
    Layer steps[STEPS + 1];
    CudaCase params;
    // Record records[ITERS + STEPS + 1];
    size_t nxt;
    real t;
    size_t timestamp;
    size_t curr_case_num;
    bool ended;
};

template <unsigned int STEPS, unsigned int ITERS>
__device__ void UniStepAdams(ThreadContext<STEPS, ITERS> &c, size_t nxt, real dt);

__device__ void
OneStepAdams(Layer &res, const Layer &f0, real dt, CudaCase &params)
{
    res.Set(f0.V_ + f0.fV_ * dt, f0.gamma_ + f0.fgamma_ * dt, f0.h_ + f0.fh_ * dt, f0.l_ + f0.fl_ * dt,
            f0.M_ + f0.fM_ * dt, params);
}

template <>
__device__ void
UniStepAdams(ThreadContext<1, ITERS_PER_KERNEL> &c, size_t nxt, real dt)
{
    OneStepAdams(c.steps[nxt], c.steps[(nxt + 1) & 1], dt, c.params);
}

__device__ void
TwoStepAdams(Layer &res, const Layer &f1, const Layer &f0, real dt, CudaCase &params)
{
    const auto c1 = (real) 1.5;
    const auto c0 = -(real) 0.5;

    res.Set(f1.V_ + (c1 * f1.fV_ + c0 * f0.fV_) * dt, f1.gamma_ + (c1 * f1.fgamma_ + c0 * f0.fgamma_) * dt,
            f1.h_ + (c1 * f1.fh_ + c0 * f0.fh_) * dt, f1.l_ + (c1 * f1.fl_ + c0 * f0.fl_) * dt,
            f1.M_ + (c1 * f1.fM_ + c0 * f0.fM_) * dt, params);
}

template <>
__device__ void
UniStepAdams(ThreadContext<2, ITERS_PER_KERNEL> &c, size_t nxt, real dt)
{
    TwoStepAdams(c.steps[nxt], c.steps[(nxt + 1) % 3], c.steps[(nxt + 2) % 3], dt, c.params);
}

__device__ void
ThreeStepAdams(Layer &res, const Layer &f2, const Layer &f1, const Layer &f0, real dt, CudaCase &params)
{
    const auto c2 = (real) 23 / 12;
    const auto c1 = -(real) 16 / 12;
    const auto c0 = (real) 5 / 12;

    res.Set(f2.V_ + (c2 * f2.fV_ + c1 * f1.fV_ + c0 * f0.fV_) * dt,
            f2.gamma_ + (c2 * f2.fgamma_ + c1 * f1.fgamma_ + c0 * f0.fgamma_) * dt,
            f2.h_ + (c2 * f2.fh_ + c1 * f1.fh_ + c0 * f0.fh_) * dt,
            f2.l_ + (c2 * f2.fl_ + c1 * f1.fl_ + c0 * f0.fl_) * dt,
            f2.M_ + (c2 * f2.fM_ + c1 * f1.fM_ + c0 * f0.fM_) * dt, params);
}

template <>
__device__ void
UniStepAdams(ThreadContext<3, ITERS_PER_KERNEL> &c, size_t nxt, real dt)
{
    ThreeStepAdams(c.steps[nxt], c.steps[(nxt + 1) & 3], c.steps[(nxt + 2) & 3], c.steps[(nxt + 3) & 3], dt, c.params);
}

template <unsigned int STEPS, unsigned int ITERS>
__device__ void
InitContext(ThreadContext<STEPS, ITERS> &c, CudaCase &curr_case, real dt, size_t grid_idx, Record *&curr_record)
{
    c.params = curr_case;
    c.steps[STEPS].Set(curr_case.v0_, curr_case.gamma0_, curr_case.h0_, 0.0, // Case::l0()
                       curr_case.m0_, c.params);
    *curr_record = {
        0.0, c.steps[STEPS].M_, c.steps[STEPS].V_, c.steps[STEPS].h_, c.steps[STEPS].l_, c.steps[STEPS].gamma_};
    curr_record++;
    OneStepAdams(c.steps[STEPS - 1], c.steps[STEPS], dt, c.params);
    *curr_record = {dt,
                    c.steps[STEPS - 1].M_,
                    c.steps[STEPS - 1].V_,
                    c.steps[STEPS - 1].h_,
                    c.steps[STEPS - 1].l_,
                    c.steps[STEPS - 1].gamma_};
    curr_record++;
    if (STEPS >= 2) {
        TwoStepAdams(c.steps[STEPS - 2], c.steps[STEPS - 1], c.steps[STEPS], dt, c.params);
        *curr_record = {dt,
                        c.steps[STEPS - 2].M_,
                        c.steps[STEPS - 2].V_,
                        c.steps[STEPS - 2].h_,
                        c.steps[STEPS - 2].l_,
                        c.steps[STEPS - 2].gamma_};
        curr_record++;
    }

    if (STEPS >= 3) {
        ThreeStepAdams(c.steps[STEPS - 3], c.steps[STEPS - 2], c.steps[STEPS - 1], c.steps[STEPS], dt, c.params);
        *curr_record = {dt,
                        c.steps[STEPS - 1].M_,
                        c.steps[STEPS - 1].V_,
                        c.steps[STEPS - 1].h_,
                        c.steps[STEPS - 1].l_,
                        c.steps[STEPS - 1].gamma_};
        curr_record++;
    }
}

template <unsigned int STEPS, unsigned int ITERS>
__global__ void
AdamsKernel(ThreadContext<STEPS, ITERS> *contexts, CudaCase *problems, real dt, real timeout, real *timestamps,
            size_t n_timestamps, real *functional_args, Record *records, unsigned int *active_threads)
{
    size_t grid_idx = blockIdx.x * blockDim.x + threadIdx.x;
    ThreadContext<STEPS, ITERS> &c = contexts[grid_idx];

    if (c.ended) {
        return;
    }
    if (c.curr_case_num == CASES_PER_THREAD) {
        c.ended = true;
        atomicSub(active_threads, 1);
    }

    CudaCase &current_case = problems[grid_idx * CASES_PER_THREAD + c.curr_case_num];
    Record *curr_record = records + (ITERS_PER_KERNEL + STEPS) * grid_idx;
    real t;
    size_t record_num; // can be calculated dynamically
    real *V_arg, *h_arg;
    size_t nxt, timestamp;

    // restore context
    if (c.t == 0.0) { // new case
        InitContext(c, current_case, dt, grid_idx, curr_record);
        t = dt * (real) STEPS;
        record_num = STEPS + 1;
        V_arg = functional_args + (grid_idx * n_timestamps * 2);
        h_arg = functional_args + (grid_idx * n_timestamps * 2) + n_timestamps;
        nxt = STEPS;
        timestamp = 0;
    } else { // old case
        t = c.t;
        record_num = c.t / dt;
        V_arg = functional_args + (grid_idx * n_timestamps * 2) + c.timestamp;
        h_arg = functional_args + (grid_idx * n_timestamps * 2) + n_timestamps + c.timestamp;
        nxt = c.nxt;
        timestamp = c.timestamp;
    }

    // main loop
    unsigned int iters_count = 0;
    while (t < timeout && iters_count < ITERS) {
        // If necessery, update the functional's arguments
        if (timestamp < n_timestamps && t >= timestamps[timestamp]) {
            const auto &step = c.steps[(nxt + 1) % (STEPS + 1)];
            V_arg[timestamp] = step.V_;
            h_arg[timestamp] = step.h_;
            timestamp++;
        }

        // Compute values for the next step, store them
        UniStepAdams<STEPS, ITERS>(c, nxt, dt);
        auto M = c.steps[nxt].M_;
        auto h = c.steps[nxt].h_;

        t += dt;
        *curr_record = {t, M, c.steps[nxt].V_, h, c.steps[nxt].l_, c.steps[nxt].gamma_};
        curr_record++;
        record_num++;
        nxt = (nxt + STEPS) % (STEPS + 1);

        // Check, should we stop the simulation?
        if (M <= (real) 0.01 || h <= (real) 0.0) {
            c.t = 0.0;
            c.curr_case_num++;
        }
        iters_count++;
    }
};

template <unsigned int STEPS, unsigned int ITERS>
void
CudaAdamsMethod(const std::vector<Case> &problem_vector, const IFunctional &functional, real dt, real timeout,
                IResultFormatter &results)
{
    size_t count_of_threads = problem_vector.size() / CASES_PER_THREAD;

    // preparing buffers

    // thread contexts (inner, frequently used)
    ThreadContext<STEPS, ITERS> *dev_thread_sandbox_arr;
    size_t thread_sandbox_arr_size = sizeof(ThreadContext<STEPS, ITERS>) * count_of_threads;
    HANDLE_ERROR(cudaMalloc(&dev_thread_sandbox_arr, thread_sandbox_arr_size));

    // problems (read-only)
    CudaCase *dev_problems;
    size_t problems_size = sizeof(CudaCase) * problem_vector.size();
    HANDLE_ERROR(cudaMalloc(&dev_problems, problems_size));
    HANDLE_ERROR(cudaMemcpy(dev_problems, problem_vector.data(), problems_size, cudaMemcpyHostToDevice)); // damn reinterpret cast

    // meteorite's timestamps (read-only)
    const real *timestamps = nullptr;
    size_t n_timestamps = 0;
    functional.GetTimeStamps(n_timestamps, timestamps);
    real *dev_timestamps;
    HANDLE_ERROR(cudaMalloc(&dev_timestamps, sizeof(real) * n_timestamps));
    HANDLE_ERROR(cudaMemcpy(dev_timestamps, timestamps, sizeof(real) * n_timestamps, cudaMemcpyHostToDevice));

    // vector of (V_arg[] and h_arg[] for functional)'s arrays
    // load to vector every iteration (if timestamp = n_timestamps ???)
    // damn fragmentation
    std::vector<real *> functional_args;
    size_t functional_args_size = n_timestamps * sizeof(real) * 2 * count_of_threads;
    real *dev_functional_args;
    HANDLE_ERROR(cudaMalloc(&dev_functional_args, functional_args_size));

    // vector of Record's arrays
    // load to vector every iteration
    std::vector<Record *> records;
    size_t records_size = sizeof(Record) * (ITERS_PER_KERNEL + STEPS) * count_of_threads;
    Record *dev_records;
    HANDLE_ERROR(cudaMalloc(&dev_records, records_size));

    // active threads atomic counter (make reduction)
    unsigned int active_threads = count_of_threads, *dev_active_threads;
    HANDLE_ERROR(cudaMalloc(&dev_active_threads, sizeof(active_threads)));
    HANDLE_ERROR(cudaMemcpy(dev_active_threads, &active_threads, sizeof(unsigned int), cudaMemcpyHostToDevice));

    std::cout << "all cases:" << problem_vector.size() << " threads:" << count_of_threads << "\n";
    getchar();
    int iter = 0;
    while (active_threads) {
        iter++;
        std::cout << "while iter:" << iter << " active_threads:" << active_threads << "\n";

        AdamsKernel<STEPS, ITERS>
        <<<problem_vector.size() / CASES_PER_THREAD / 32, 32>>>
        (dev_thread_sandbox_arr, dev_problems, dt, timeout, dev_timestamps, n_timestamps,
         dev_functional_args, dev_records, dev_active_threads);
  
        functional_args.push_back((real *)malloc(functional_args_size));
        HANDLE_ERROR(cudaMemcpy(functional_args[functional_args.size() - 1], dev_functional_args,
                                functional_args_size, cudaMemcpyDeviceToHost));

        records.push_back((Record *)malloc(records_size));
        HANDLE_ERROR(cudaMemcpy(records[records.size() - 1], dev_records,
                                records_size, cudaMemcpyDeviceToHost));

        HANDLE_ERROR(cudaMemcpy(&active_threads, dev_active_threads, sizeof(active_threads), cudaMemcpyDeviceToHost));
    }

    std::cout << "\nend\n";

    cudaFree(dev_thread_sandbox_arr);
    cudaFree(dev_problems);
    cudaFree(dev_timestamps);
    cudaFree(dev_functional_args);
    cudaFree(dev_records);
    cudaFree(dev_active_threads);
}

// template-specified functions compiles only this way
template void CudaAdamsMethod<1u, ITERS_PER_KERNEL>(const std::vector<Case> &, const IFunctional &, real, real,
                                                    IResultFormatter &);
template void CudaAdamsMethod<2u, ITERS_PER_KERNEL>(const std::vector<Case> &, const IFunctional &, real, real,
                                                    IResultFormatter &);
template void CudaAdamsMethod<3u, ITERS_PER_KERNEL>(const std::vector<Case> &, const IFunctional &, real, real,
                                                    IResultFormatter &);
