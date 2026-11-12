#include "linalg/gep.hpp"

#include <stdexcept>

namespace qgem {

Mat orthonormalizer(const Mat& S, double eta) {
  Eigen::SelfAdjointEigenSolver<Mat> es(S);
  const Vec& s = es.eigenvalues();
  const int n = static_cast<int>(s.size());
  int first = 0;
  while (first < n && s(first) <= eta * s(n - 1)) ++first;
  const int m = n - first;
  if (m == 0)
    throw std::runtime_error(
        "orthonormalizer: no overlap eigenvalue above threshold");
  return es.eigenvectors().rightCols(m) *
         s.tail(m).cwiseSqrt().cwiseInverse().asDiagonal();
}

GepResult solve_gep(const Mat& H, const Mat& S, double eta) {
  if (H.rows() != H.cols() || H.rows() != S.rows() || S.rows() != S.cols())
    throw std::invalid_argument(
        "solve_gep: H and S must be square and equal in size");
  Eigen::SelfAdjointEigenSolver<Mat> es(S);
  const Vec& s = es.eigenvalues();
  const int n = static_cast<int>(s.size());
  int first = 0;
  while (first < n && s(first) <= eta * s(n - 1)) ++first;
  const int m = n - first;
  if (m == 0)
    throw std::runtime_error(
        "solve_gep: no overlap eigenvalue above threshold");
  Mat X = es.eigenvectors().rightCols(m) *
          s.tail(m).cwiseSqrt().cwiseInverse().asDiagonal();
  Mat Hp = X.transpose() * H * X;
  Hp = 0.5 * (Hp + Hp.transpose());
  Eigen::SelfAdjointEigenSolver<Mat> eh(Hp);
  return {eh.eigenvalues(), X * eh.eigenvectors(), X, m, s(first), s(n - 1)};
}

GepResult solve_gep_cholesky(const Mat& H, const Mat& S) {
  Eigen::LLT<Mat> llt(S);
  if (llt.info() != Eigen::Success)
    throw std::runtime_error("solve_gep_cholesky: S not positive definite");
  Mat L = llt.matrixL();
  // H' = L^{-1} H L^{-T}
  Mat Y = L.triangularView<Eigen::Lower>().solve(H);
  Mat Hp = L.triangularView<Eigen::Lower>().solve(Y.transpose()).transpose();
  Hp = 0.5 * (Hp + Hp.transpose());
  Eigen::SelfAdjointEigenSolver<Mat> eh(Hp);
  Mat C = L.transpose().triangularView<Eigen::Upper>().solve(eh.eigenvectors());
  Eigen::SelfAdjointEigenSolver<Mat> es(S);
  return {eh.eigenvalues(),
          C,
          Mat(),
          static_cast<int>(S.rows()),
          es.eigenvalues()(0),
          es.eigenvalues()(S.rows() - 1)};
}

}  // namespace qgem
