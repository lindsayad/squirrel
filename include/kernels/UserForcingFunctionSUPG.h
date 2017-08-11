/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef USERFORCINGFUNCTIONSUPG_H
#define USERFORCINGFUNCTIONSUPG_H

#include "Kernel.h"

// Forward Declarations
class UserForcingFunctionSUPG;
class Function;

template <>
InputParameters validParams<UserForcingFunctionSUPG>();

/**
 * Define the Kernel for a user defined forcing function that looks like:
 *
 * test function * forcing function
 */
class UserForcingFunctionSUPG : public Kernel
{
public:
  UserForcingFunctionSUPG(const InputParameters & parameters);

protected:
  /**
   * Evaluate f at the current quadrature point.
   */
  Real f();

  /**
   * Computes test function * forcing function.
   */
  virtual Real computeQpResidual() override;

  Function & _func;
  const RealVectorValue & _velocity;
  const MaterialProperty<Real> & _D;
};

#endif // USERFORCINGFUNCTIONSUPG_H
