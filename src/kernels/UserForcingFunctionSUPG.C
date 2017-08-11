#include "UserForcingFunctionSUPG.h"
#include "Function.h"

template <>
InputParameters
validParams<UserForcingFunctionSUPG>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<FunctionName>("function", "The forcing function");
  params.addRequiredParam<RealVectorValue>("velocity", "velocity vector");
  params.addParam<MaterialPropertyName>("D_name", "D", "The name of the diffusivity");
  return params;
}

UserForcingFunctionSUPG::UserForcingFunctionSUPG(const InputParameters & parameters)
  : Kernel(parameters),
    _func(getFunction("function")),
    _velocity(getParam<RealVectorValue>("velocity")),
    _D(getMaterialProperty<Real>("D_name"))
{
}

Real
UserForcingFunctionSUPG::f()
{
  return _func.value(_t, _q_point[_qp]);
}

Real
UserForcingFunctionSUPG::computeQpResidual()
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

  return -PG_test * f();
}
