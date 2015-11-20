/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_ves_biases_VesBias_h
#define __PLUMED_ves_biases_VesBias_h

#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "bias/Bias.h"

#define PLUMED_VARIATIONALBIAS_INIT(ao) Action(ao),VesBias(ao)

namespace PLMD{
  class CoeffsVector;
  class CoeffsMatrix;
namespace bias{

/**
\ingroup INHERIT
Abstract base class for implementing biases the extents the normal Bias.h class
to include functions related to the variational approach.
*/

class VesBias:
public Bias
{
private:
  CoeffsVector* coeffs_;
  CoeffsVector* gradient_;
  CoeffsMatrix* hessian_;
public:
  static void registerKeywords(Keywords&);
  VesBias(const ActionOptions&ao);
  CoeffsVector* getCoeffsPtr(){return coeffs_;}
  CoeffsVector* getGradientPtr(){return gradient_;}
  CoeffsMatrix* getHessianPtr(){return hessian_;}
  void updateGradientAndHessian();
  void clearGradientAndHessian();
};

}
}

#endif
