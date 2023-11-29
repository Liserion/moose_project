//created by Armin 29.10.2020

#pragma once

#include "IntegratedBC.h"

class ConstFluxForPhiSBC;

template <>
InputParameters validParams<ConstFluxForPhiSBC>();

class ConstFluxForPhiSBC:public IntegratedBC
{
public:
    ConstFluxForPhiSBC(const InputParameters &parameters);

protected:
    virtual Real computeQpResidual() override ;
    virtual Real computeQpJacobian() override { return 0.0;};

    const Real &_I;
    const Real &_ChargeTime;
};