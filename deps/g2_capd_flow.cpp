#include "capd/mpcapdlib.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace capd;

static MpFloat widen_amount(const MpFloat& x) {
  return abs(x) * MpFloat("1e-60") + MpFloat("1e-200");
}

static MpFloat widen_lo(const MpFloat& x) { return x - widen_amount(x); }
static MpFloat widen_hi(const MpFloat& x) { return x + widen_amount(x); }

static void print_interval(const char* tag, int i, int j,
                           const MpInterval& x) {
  std::ostringstream out;
  out.precision(90);
  out << tag << " " << i;
  if (j >= 0) out << " " << j;
  out << " " << widen_lo(x.leftBound())
      << " " << widen_hi(x.rightBound()) << "\n";
  std::cout << out.str();
}

static void print_lower(const char* tag, int i, const MpFloat& x) {
  std::ostringstream out;
  out.precision(90);
  out << tag << " " << i << " " << widen_lo(x) << "\n";
  std::cout << out.str();
}

static std::string read_file(const char* path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("could not open CAPD field file");
  std::ostringstream text;
  text << input.rdbuf();
  return text.str();
}

static void update_chamber(int kind, const MpIVector& z,
                           MpFloat (&minimum)[3]) {
  MpInterval d[3];
  if (kind == 4) {
    // z=(q1,...,q8,T)
    d[0] = z[1] + z[6];
    d[1] = z[2] + z[5];
    d[2] = -(z[3] + z[4]);
  } else {
    // z=(q1,q2,q4,q5,q6,q8,T), with q3=q2 and q7=q6.
    d[0] = z[1] + z[4];
    d[1] = z[1] + z[4];
    d[2] = -(z[2] + z[3]);
  }
  for (int i = 0; i < 3; ++i) {
    MpFloat lower = d[i].leftBound();
    if (!(lower > MpFloat(0)))
      throw std::runtime_error("real G2 chamber was not certified on a CAPD step");
    if (lower < minimum[i]) minimum[i] = lower;
  }
}

int main(int argc, char* argv[]) {
  try {
    // argv: KIND FIELD ORDER TOL BITS LAM_LO LAM_HI DELTA,
    //       followed by lower/upper endpoints of every initial coordinate.
    if (argc < 10) throw std::runtime_error("insufficient arguments");
    const int kind = std::atoi(argv[1]);
    if (kind != 4 && kind != 23) throw std::runtime_error("unknown case");
    const int dimension = kind == 4 ? 9 : 7;
    if (argc != 9 + 2 * dimension)
      throw std::runtime_error("wrong number of interval endpoints");

    const int order = std::atoi(argv[3]);
    const double tolerance = std::strtod(argv[4], nullptr);
    MpFloat::setDefaultPrecision(std::atoi(argv[5]));

    const std::string field = read_file(argv[2]);
    MpIMap vector_field(field.c_str());
    const MpFloat lambda_lower(argv[6]), lambda_upper(argv[7]);
    vector_field.setParameter(
        "lam", MpInterval(widen_lo(lambda_lower), widen_hi(lambda_upper)));
    vector_field.setParameter("delta", MpInterval(MpFloat(argv[8])));

    MpIVector initial(dimension);
    for (int i = 0; i < dimension; ++i) {
      const MpFloat lower(argv[9 + 2 * i]);
      const MpFloat upper(argv[10 + 2 * i]);
      initial[i] = MpInterval(widen_lo(lower), widen_hi(upper));
    }

    MpIOdeSolver solver(vector_field, order);
    solver.setAbsoluteTolerance(tolerance);
    solver.setRelativeTolerance(tolerance);
    MpITimeMap time_map(solver);
    MpC1Rect2Set set(initial);
    MpIVector result(dimension);
    MpIMatrix derivative(dimension, dimension);
    MpFloat margins[3] = {MpFloat("1e300"), MpFloat("1e300"),
                          MpFloat("1e300")};
    update_chamber(kind, initial, margins);

    long steps = 0;
    time_map.stopAfterStep(true);
    do {
      result = time_map(MpInterval(1), set, derivative);
      ++steps;
      const MpIOdeSolver::SolutionCurve& curve = solver.getCurve();
      MpIVector tube = curve(MpInterval(curve.getLeftDomain(),
                                        curve.getRightDomain()));
      update_chamber(kind, tube, margins);
      if (steps % 40 == 0 || time_map.completed()) {
        MpFloat width(0);
        for (int i = 0; i < dimension; ++i) {
          MpFloat w = result[i].rightBound() - result[i].leftBound();
          if (w > width) width = w;
        }
        std::ostringstream progress;
        progress.precision(12);
        progress << "CAPD: s=" << set.getCurrentTime().rightBound();
        progress.precision(4);
        progress << " diameter<=" << width << " steps=" << steps << "\n";
        std::cerr << progress.str();
        std::cerr.flush();
      }
    } while (!time_map.completed());

    for (int i = 0; i < dimension; ++i)
      print_interval("STATE", i, -1, result[i]);
    for (int i = 0; i < dimension; ++i)
      for (int j = 0; j < dimension; ++j)
        print_interval("MAT", i, j, derivative[i][j]);
    for (int i = 0; i < 3; ++i) print_lower("MARGIN", i, margins[i]);
    std::cout << "STEPS " << steps << "\nOK\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "CAPD ERROR: " << error.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "CAPD ERROR: unknown exception\n";
    return 1;
  }
}
