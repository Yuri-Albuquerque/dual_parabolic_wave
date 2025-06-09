#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cstdint>
#include <eigen3/Eigen/Eigen>
#include <omp.h>

using namespace std;
namespace py = pybind11;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;

MatrixXd resizeSlice(py::array_t<double> mat, int index, int nz, int nx)
{
  auto u = mat.mutable_unchecked<3>();
  MatrixXd resize = MatrixXd::Zero(nz + 2, nx + 2);
  for (ssize_t i = 0; i < nz; i++)
  {
    for (ssize_t j = 0; j < nx; j++)
    {
      resize(i + 1, j + 1) = u(i, j, index);
    };
  };

  return resize;
};

py::array_t<double> solve_wv(py::array_t<double> parameters, py::array_t<double> w, py::array_t<double> vel, py::array_t<double> eta, py::array_t<double> S)
{
  auto u = w.mutable_unchecked<3>();
  auto c = vel.unchecked<2>();
  auto damp = eta.unchecked<2>();
  auto source = S.unchecked<3>();
  auto par = parameters.unchecked<1>();
  double xMin = par(0);
  double xMax = par(1);
  double zMin = par(2);
  double zMax = par(3);
  double tMin = par(4);
  double tMax = par(5);
  double hx = par(6);
  double hz = par(7);
  double ht = par(8);
  int nz = u.shape(0);
  int nx = u.shape(1);
  int nt = u.shape(2);
  VectorXd Cxx = VectorXd::Zero(3);
  VectorXd Czz = VectorXd::Zero(3);
  Cxx << 1.0, -2.0, 1.0;
  Czz = Cxx;
  double dzz;
  double dxx;
  double q0;
  double q1;
  double q2;
  double q3;
  double init_cond_u0 = 0.0;
#pragma omp parallel for
  for (ssize_t i = 0; i < nz; i++)
  {
    for (ssize_t j = 0; j < nx; j++)
    {
      u(i, j, 0) = init_cond_u0;
      u(i, j, 1) = u(i, j, 0) + ht * init_cond_u0;
    };
  };

#pragma omp parallel for ordered schedule(static)
  for (ssize_t t = 1; t < nt - 1; t++)
  {
    auto resized = resizeSlice(w, t, nz, nx);
    for (ssize_t i = 0; i < nz; i++)
    {
      for (ssize_t j = 0; j < nx; j++)
      {
        q0 = c(i, j) * c(i, j) * ht;
        q1 = c(i, j) * c(i, j) * ht * ht;
        q2 = ((c(i, j) * ht) / hx) * ((c(i, j) * ht) / hx);
        q3 = ((c(i, j) * ht) / hz) * ((c(i, j) * ht) / hz);
        dzz = Czz(0) * resized(i, j + 1) + Czz(1) * resized(i + 1, j + 1) + Czz(2) * resized(i + 2, j + 1);
        dxx = Cxx(0) * resized(i + 1, j) + Cxx(1) * resized(i + 1, j + 1) + Cxx(2) * resized(i + 1, j + 2);
        if (i == 0)
        {
          u(0, j, t + 1) = (-1.0 * (u(0, j, t - 1) - 2.0 * u(0, j, t)) + q0 * damp(0, j) * u(0, j, t - 1) + q1 * source(0, j, t + 1) + q2 * dxx + q3 * 2.0 * (u(1, j, t) - u(0, j, t))) * (1.0 / (1.0 + damp(0, j) * q0));
        }
        else
        {
          u(i, j, t + 1) = (-1.0 * (u(i, j, t - 1) - 2.0 * u(i, j, t)) + q0 * damp(i, j) * u(i, j, t - 1) + q1 * source(i, j, t + 1) + q2 * dxx + q3 * dzz) * (1.0 / (1.0 + damp(i, j) * q0));
        };
      };
    };
  };
  return w;
};

py::array_t<double> reverse_time(py::array_t<double> p, py::array_t<double> p_tilde)
{
  auto P = p.mutable_unchecked<3>();
  auto ptilde = p_tilde.mutable_unchecked<3>();
  int nz = ptilde.shape(0);
  int nx = ptilde.shape(1);
  int nt = ptilde.shape(2);
#pragma omp parallel for ordered schedule(static)
  for (ssize_t t = 0; t < nt; t++)
  {
    for (ssize_t i = 0; i < nz; i++)
    {
      for (ssize_t j = 0; j < nx; j++)
      {
        P(i, j, t) = ptilde(i, j, (nt - 1) - t);
      };
    };
  };

  return p;
};

// The level set method for 2d - domain
py::array_t<double> hj(py::array_t<double> v1, py::array_t<double> v2, py::array_t<double> phiraw, py::array_t<double> parameters)
{
  auto phi = phiraw.mutable_unchecked<2>();
  auto V1 = v1.unchecked<2>();
  auto V2 = v2.unchecked<2>();
  auto par = parameters.unchecked<1>();
  double beta = par(0);
  double itermax = par(1);
  int nz = phi.shape(0);
  int nx = phi.shape(1);
  double h = min(1.0 / nz, 1.0 / nx);
  double maxv;
  double dt;
  MatrixXd maxv_mat = MatrixXd::Zero(nz, nx);
  MatrixXd g = MatrixXd::Zero(nz, nx);
  MatrixXd Dxp = MatrixXd::Zero(nz, nx);
  MatrixXd Dyp = MatrixXd::Zero(nz, nx);
  // -------- x --------
  MatrixXd Dxm = MatrixXd::Zero(nz, nx);
  MatrixXd Dym = MatrixXd::Zero(nz, nx);
  for (ssize_t k = 0; k < itermax; k++)
  {
    // Compute first-order terms in the HJ equation
    Dxp.setZero();
    Dyp.setZero();
    Dxm.setZero();
    Dym.setZero();
    /* Be careful here,
       the x-axis is in column
       the y-axis is in rows
     */
#pragma omp parallel for
    for (ssize_t i = 0; i < nz; i++)
    {
      for (ssize_t j = 0; j < nx - 1; j++)
      {
        Dxp(i, j) = phi(i, j + 1) * nx - phi(i, j) * nx;
      };
    };
#pragma omp parallel for
    for (ssize_t i = 0; i < nz; i++)
    {
      for (ssize_t j = 1; j < nx; j++)
      {
        Dxm(i, j) = phi(i, j) * nx - phi(i, j - 1) * nx;
      };
    };
#pragma omp parallel for
    for (ssize_t i = 0; i < nz - 1; i++)
    {
      for (ssize_t j = 0; j < nx; j++)
      {
        Dyp(i, j) = phi(i + 1, j) * nz - phi(i, j) * nz;
      };
    };
#pragma omp parallel for
    for (ssize_t i = 1; i < nz - 1; i++)
    {
      for (ssize_t j = 0; j < nx; j++)
      {
        Dym(i, j) = phi(i, j) * nz - phi(i - 1, j) * nz;
      };
    };
#pragma omp parallel for
    for (ssize_t i = 0; i < nz; i++)
    {
      Dxp(i, nx - 1) = Dxp(i, nx - 2);
      Dxm(i, 0) = Dxm(i, 1);
    };
#pragma omp parallel for
    for (ssize_t j = 0; j < nx; j++)
    {
      Dyp(nz - 1, j) = Dxp(nz - 2, j);
      Dym(0, j) = Dym(1, j);
    };
#pragma omp parallel for
    for (ssize_t i = 0; i < nz; i++)
    {
      for (ssize_t j = 0; j < nx; j++)
      {
        g(i, j) = 0.5 * (V1(i, j) * (Dxp(i, j) + Dxm(i, j)) + V2(i, j) * (Dyp(i, j) + Dym(i, j))) - 0.5 * abs(V1(i, j)) * (Dxp(i, j) - Dxm(i, j)) - 0.5 * abs(V2(i, j)) * (Dyp(i, j) - Dym(i, j));
        maxv_mat(i, j) = abs(V1(i, j)) + abs(V2(i, j));
      };
    };
    maxv = maxv_mat.maxCoeff();
    dt = beta * h / (2 * sqrt(2) * maxv);
#pragma omp parallel for ordered schedule(static)
    for (ssize_t i = 0; i < nz; i++)
    {
      for (ssize_t j = 0; j < nx; j++)
      {
        phi(i, j) = phi(i, j) - dt * g(i, j);
      };
    };
  };
  return phiraw;
};

// The absorbing damping layer for the acoustic wave equation
py::array_t<double> damping_function(py::array_t<double> eta_raw, py::array_t<double> parameters)
{
  auto par = parameters.unchecked<1>();
  double xMin = par(0);
  double xMax = par(1);
  double zMin = par(2);
  double zMax = par(3);
  double tMin = par(4);
  double tMax = par(5);
  double hx = par(6);
  double hz = par(7);
  double ht = par(8);
  double dmp_xmin = par(10);
  double dmp_xmax = par(11);
  double dmp_zmax = par(12);
  //
  int nz = (zMax - zMin) / hz;
  int nx = (xMax - xMin) / hx;
  int nt = (tMax - tMin) / ht;
  VectorXd grid_x, grid_z;
  grid_x.setLinSpaced(nx, xMin, xMax);
  grid_z.setLinSpaced(nz, zMin, zMax);
  auto dmp_func = [](auto x, auto lim) { return 1e4 * pow(x - lim, 2); };
  //
  // MatrixXd eta = MatrixXd::Zero(nz, nx);
  auto eta = eta_raw.mutable_unchecked<2>();
  double a1 = (zMax - dmp_zmax) / (xMin - dmp_xmin);
  double a2 = (zMax - dmp_zmax) / (xMax - dmp_xmax);
  double b1 = (dmp_zmax * xMin - zMax * dmp_xmin) / (xMin - dmp_xmin);
  double b2 = (dmp_zmax * xMax - zMax * dmp_xmax) / (xMax - dmp_xmax);
  //
#pragma omp parallel for
  for (size_t i = 0; i < nz; i++)
  {
    for (size_t j = 0; j < nx; j++)
    {
      if ((grid_x(j) < dmp_xmin && grid_z(i) < dmp_zmax) || (grid_z(i) >= dmp_zmax && grid_z(i) <= a1 * grid_x(j) + b1))
      {
        eta(i, j) = dmp_func(grid_x(j), dmp_xmin);
      }else if ((grid_x(j) > dmp_xmax && grid_z(i) < dmp_zmax) || (grid_z(i) >= dmp_zmax && grid_z(i) <= a2 * grid_x(j) + b2))
      {
        eta(i, j) = dmp_func(grid_x(j), dmp_xmax);
      }else  if (grid_z(i) >= dmp_zmax && ((grid_x(j) >= 0 && grid_x(j) <= 1.0 / a1 * (grid_z(i) - b1)) || grid_x(j) <= 1.0 / a2 * (grid_z(i) - b2)))
      {
        eta(i, j) = dmp_func(grid_z(i), dmp_zmax);
      }
    }
  }
  return eta_raw;
}

// ---------------- X -------------------
PYBIND11_MODULE(wave_solver, m)
{
  m.doc() = R"pbdoc(
        Pybind11 wave solver plugin
        -----------------------
     
        .. currentmodule:: wave_solver
     
        .. autosummary::
           :toctree: _generate
         
           solver
    )pbdoc";

  m.def("solve_wv", &solve_wv, R"pbdoc(

        Gives the solution through finite difference of acoustic wave PDE rho*u_tt - eta*u_t + div(grad(u) = f. Given the velocity field c = 1/rho**2, some damping function eta and the source f. 
    )pbdoc");

  m.def("reverse_time", &reverse_time, R"pbdoc(

        Invert order of the third dimension of some tree-dimensional Matrix;
    )pbdoc");

  m.def(
      "subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
     
        Some other explanation about the subtract function.
    )pbdoc");

  m.def("hj", &hj, R"pbdoc(
          
        Compute level set function;
    )pbdoc");

  m.def("damping_function", &damping_function, R"pbdoc(
          
        Compute the absorbing sponge layer for the acoustic wave equation;
    )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = VERSION_INFO;
#else
  m.attr("__version__") = "dev";
#endif
}