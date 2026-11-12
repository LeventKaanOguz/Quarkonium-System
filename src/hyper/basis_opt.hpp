#pragma once
// Basis hyperparameters (notes 7.12): they are not physics, so they are chosen
// variationally and their footprint on every result is reported as a numerical
// error sigma_num from nested doubling n -> 2n-1.
#include <functional>
#include <vector>

#include "data/dataset.hpp"
#include "models/model.hpp"
#include "solve/spectrum.hpp"

namespace qgem {

struct BasisOptResult {
  BasisSpec spec;
  double objective = 0;
  double start_objective = 0;
  int evaluations = 0;
};

/// Variational objective of a whole dataset at theta: for every channel the
/// data name, the sum of its lowest K levels (K = highest level the data use
/// + 1). By the Ky Fan trace principle this sum is an upper bound on its exact
/// value, so a lower value is a more converged basis. Infinity if a channel
/// collapses (a level below zero mass: H not bounded below), a level the data
/// need is missing, or the ranges leave [r1_min, rmax_max].
double dataset_trace(const Model& model, const Vec& theta, const Dataset& data,
                     const BasisSpec& spec);

/// Minimise dataset_trace over (ln r1, ln rmax) at fixed n (Nelder-Mead),
/// within r1 in [0.002, 1] and rmax in [3, 400] GeV^-1.
BasisOptResult optimize_dataset_basis(const Model& model, const Vec& theta,
                                      const Dataset& data, BasisSpec start);

/// Mass of every datum at theta in one basis (NaN where the level is missing
/// or the channel collapses).
Vec dataset_levels(const Model& model, const Vec& theta, const Dataset& data,
                   const BasisSpec& spec);

/// Choice of the basis size by nested convergence. The sizes n_0, 2n_0 - 1,
/// 4n_0 - 3, ... at fixed (r1, rmax) form nested bases (each contains the
/// previous one), so by Hylleraas-Undheim-MacDonald every level can only
/// decrease along the sequence. The chosen size is the smallest n_k whose
/// largest level change to n_{k+1} is below `target` (GeV).
struct BasisSizeChoice {
  std::vector<int> n;              ///< the nested sizes
  std::vector<double> max_change;  ///< max_i |M_i(n_k) - M_i(n_k+1)|, GeV
  std::vector<std::string> worst;  ///< label of the datum giving max_change
  double max_increase = 0;    ///< largest rise of a level along the sequence
                              ///< (0 for an exactly nested space; > 0 only
                              ///< through eigenvectors dropped below eta)
  int chosen = 0;             ///< 0: no size meets the target
  double gap_to_largest = 0;  ///< max_i |M_i(chosen) - M_i(largest size)|, GeV
};
BasisSizeChoice choose_basis_size(const Model& model, const Vec& theta,
                                  const Dataset& data, BasisSpec spec, int n0,
                                  int steps, double target);

struct NestedSequence {
  std::vector<int> n;         ///< basis sizes n, 2n-1, 4n-3, ...
  std::vector<double> value;  ///< the observable at each size
  double sigma_num = 0;       ///< |P(last) - P(previous)|
  double aitken =
      0;  ///< Aitken limit (valid if the last three differences are geometric)
  bool monotone = false;  ///< for energies: HUM says it must decrease
};

/// An observable evaluated on a solved channel (e.g. a mass, or |R(0)|^2
/// through obs/).
using Observable = std::function<double(const SolvedChannel&, MatrixCache&)>;

/// Evaluate `obs` on the nested sequence of bases starting from `start`
/// (`steps` doublings, >= 2).
NestedSequence nested_convergence(const Model& model, const Vec& theta,
                                  const Channel& ch, const BasisSpec& start,
                                  const Observable& obs, int steps = 3);

}  // namespace qgem
