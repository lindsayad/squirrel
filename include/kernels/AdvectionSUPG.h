/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ADVECTIONSUPG_H
#define ADVECTIONSUPG_H

#include "Kernel.h"

// Forward Declarations
class AdvectionSUPG;

template <>
InputParameters validParams<AdvectionSUPG>();

class AdvectionSUPG : public Kernel
{
public:
  AdvectionSUPG(const InputParameters & parameters);

  virtual ~AdvectionSUPG() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  const RealVectorValue & _velocity;
  const MaterialProperty<Real> & _D;
};

#endif
