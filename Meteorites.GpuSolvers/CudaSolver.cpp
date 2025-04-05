#include "CudaSolver.h"
#include "CudaAdamsMethod.h"

constexpr unsigned int ITERS_PER_KERNEL = 10000;

void
CudaSolver::Configure(NumericalAlgorithm alg, real dt, real timeout)
{
    size_t steps = 0;
    switch (alg) {
    case NumericalAlgorithm::ONE_STEP_ADAMS:
        steps = 1;
        break;

    case NumericalAlgorithm::TWO_STEP_ADAMS:
        steps = 2;
        break;

    case NumericalAlgorithm::THREE_STEP_ADAMS:
        steps = 3;
        break;

    default:
        throw std::runtime_error("unsupported numerical algorithm");
    }

    if (dt <= (real) 0.0 || timeout < dt * (steps + 1)) {
        throw std::runtime_error("wrong time step and/or timeout");
    }

    adams_steps_ = steps;
    dt_ = dt;
    timeout_ = timeout;
}

void
CudaSolver::Solve(ICaseGenerator &generator, const IFunctional &functional, IResultFormatter &results)
{
    std::vector<Case> problem_vector;
    Case problem;
    //while (generator.Next(problem)) {
    for (int i = 0; i < 1024; i++) {
        generator.Next(problem);
        problem_vector.emplace_back(std::move(problem));
    }

    switch (AdamsSteps()) {
    case 1:
        CudaAdamsMethod<1u, ITERS_PER_KERNEL>(problem_vector, functional, Dt(), Timeout(), results);
        break;

    case 2:
        CudaAdamsMethod<2u, ITERS_PER_KERNEL>(problem_vector, functional, Dt(), Timeout(), results);
        break;

    case 3:
        CudaAdamsMethod<3u, ITERS_PER_KERNEL>(problem_vector, functional, Dt(), Timeout(), results);
        break;

    default:
        assert(false);
    }
}
