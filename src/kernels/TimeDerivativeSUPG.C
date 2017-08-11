#include "TimeDerivativeSUPG.h"

template <>
InputParameters
validParams<TimeDerivativeSUPG>()
{
  InputParameters params = validParams<TimeDerivative>();
  params.addRequiredParam<RealVectorValue>("velocity", "velocity vector");
  params.addParam<MaterialPropertyName>("D_name", "D", "The name of the diffusivity");
  return params;
}

TimeDerivativeSUPG::TimeDerivativeSUPG(const InputParameters & parameters)
  : TimeDerivative(parameters),
    _velocity(getParam<RealVectorValue>("velocity")),
    _D(getMaterialProperty<Real>("D_name"))
{
}

Real
TimeDerivativeSUPG::computeQpResidual()
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

  return TimeDerivative::computeQpResidual() + PG_test * _u_dot[_qp];
}

Real
TimeDerivativeSUPG::computeQpJacobian()
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

  return TimeDerivative::computeQpJacobian() + PG_test * _du_dot_du[_qp] * _phi[_j][_qp];
}
