#include "InterTemperatureAdvection.h"

template <>
InputParameters
validParams<InterTemperatureAdvection>()
{
  InputParameters params = validParams<InterfaceKernel>();
  params.addRequiredParam<RealVectorValue>("velocity", "Velocity vector");
  params.addClassDescription("DG upwinding for the convection");
  params.addParam<Real>("heat_source", "Should have units of power/length^2");
  return params;
}

InterTemperatureAdvection::InterTemperatureAdvection(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _velocity(getParam<RealVectorValue>("velocity")),
    _rho(getMaterialProperty<Real>("rho")),
    _cp(getMaterialProperty<Real>("cp")),
    _heat_source(getParam<Real>("heat_source"))
{
}

Real
InterTemperatureAdvection::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  Real vdotn = _rho[_qp] * _cp[_qp] * _velocity * _normals[_qp];

  switch (type)
  {
    case Moose::Element:
      if (vdotn >= 0)
        r += vdotn * _u[_qp] * _test[_i][_qp];
      else
        r += vdotn * _neighbor_value[_qp] * _test[_i][_qp] - _heat_source;
      break;

    case Moose::Neighbor:
      if (vdotn >= 0)
        r -= vdotn * _u[_qp] * _test_neighbor[_i][_qp] - _heat_source;
      else
        r -= vdotn * _neighbor_value[_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}

Real
InterTemperatureAdvection::computeQpJacobian(Moose::DGJacobianType type)
{
  Real r = 0;

  Real vdotn = _rho[_qp] * _cp[_qp] * _velocity * _normals[_qp];

  switch (type)
  {
    case Moose::ElementElement:
      if (vdotn >= 0)
        r += vdotn * _phi[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::ElementNeighbor:
      if (vdotn < 0)
        r += vdotn * _phi_neighbor[_j][_qp] * _test[_i][_qp];
      break;

    case Moose::NeighborElement:
      if (vdotn >= 0)
        r -= vdotn * _phi[_j][_qp] * _test_neighbor[_i][_qp];
      break;

    case Moose::NeighborNeighbor:
      if (vdotn < 0)
        r -= vdotn * _phi_neighbor[_j][_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}
