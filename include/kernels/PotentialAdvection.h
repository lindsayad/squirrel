#ifndef POTENTIALADVECTION_H_
#define POTENTIALADVECTION_H_

#include "Kernel.h"

class PotentialAdvection;

template <>
InputParameters validParams<PotentialAdvection>();

class PotentialAdvection : public Kernel
{
public:
  PotentialAdvection(const InputParameters & parameters);
  virtual ~PotentialAdvection();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  const unsigned int _potential_id;
  const Real _sgn;
  VariableGradient _default;
  const VariableGradient & _grad_potential;
};

#endif // POTENTIALADVECTION_H
