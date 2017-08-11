#ifndef TIMEDERIVATIVESUPG_H
#define TIMEDERIVATIVESUPG_H

#include "TimeDerivative.h"

// Forward Declaration
class TimeDerivativeSUPG;

template <>
InputParameters validParams<TimeDerivativeSUPG>();

class TimeDerivativeSUPG : public TimeDerivative
{
public:
  TimeDerivativeSUPG(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const RealVectorValue & _velocity;
  const MaterialProperty<Real> & _D;
};

#endif // TIMEDERIVATIVESUPG_H
