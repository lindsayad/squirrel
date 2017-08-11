#include "AdvectionSUPG.h"

template <>
InputParameters
validParams<AdvectionSUPG>()
{
  InputParameters params = validParams<Kernel>();

  params.addRequiredParam<RealVectorValue>("velocity", "velocity vector");
  params.addParam<MaterialPropertyName>("D_name", "D", "The name of the diffusivity");

  return params;
}

AdvectionSUPG::AdvectionSUPG(const InputParameters & parameters)
  : Kernel(parameters),
    _velocity(getParam<RealVectorValue>("velocity")),
    _D(getMaterialProperty<Real>("D_name"))
{
}

Real
AdvectionSUPG::computeQpResidual()
{
  Real alpha;
  if (_D[_qp] < std::numeric_limits<double>::epsilon())
    alpha = 1;
  else
  {
    Real Pe = _current_elem->hmax() * _velocity.norm() / (2. * _D[_qp]);
    alpha = 1. / std::tanh(Pe) - 1. / Pe;
  }
  Real PG_test =
      alpha * _current_elem->hmax() / 2. * _velocity / _velocity.norm() * _grad_test[_i][_qp];
  Real convective_part = PG_test * _velocity * _grad_u[_qp];

  return convective_part;
}

Real
AdvectionSUPG::computeQpJacobian()
{
  Real alpha;
  if (_D[_qp] < std::numeric_limits<double>::epsilon())
    alpha = 1;
  else
  {
    Real Pe = _current_elem->hmax() * _velocity.norm() / (2. * _D[_qp]);
    alpha = 1. / std::tanh(Pe) - 1. / Pe;
  }
  Real PG_test =
      alpha * _current_elem->hmax() / 2. * _velocity / _velocity.norm() * _grad_test[_i][_qp];
  Real convective_part = PG_test * _velocity * _grad_phi[_j][_qp];

  return convective_part;
}
