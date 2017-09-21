/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The ves-code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of ves-code module.

   The ves-code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The ves-code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the ves-code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "Optimizer.h"
#include "CoeffsVector.h"

#include "core/ActionRegister.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_OPTIMIZER OPT_STEEPEST_DECENT
/*
Steepest decent optimizer.

\par Examples

*/
//+ENDPLUMEDOC

class Opt_SteepestDecent : public Optimizer {

public:
  static void registerKeywords(Keywords&);
  explicit Opt_SteepestDecent(const ActionOptions&);
  void coeffsUpdate(const unsigned int c_id = 0);
};


PLUMED_REGISTER_ACTION(Opt_SteepestDecent,"OPT_STEEPEST_DECENT")


void Opt_SteepestDecent::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);
  Optimizer::useFixedStepSizeKeywords(keys);
  Optimizer::useMultipleWalkersKeywords(keys);
  Optimizer::useMaskKeywords(keys);
}


Opt_SteepestDecent::Opt_SteepestDecent(const ActionOptions&ao):
  PLUMED_VES_OPTIMIZER_INIT(ao)
{
  checkRead();
}


void Opt_SteepestDecent::coeffsUpdate(const unsigned int c_id) {
  Coeffs(c_id) += - StepSize(c_id) * CoeffsMask(c_id) * Gradient(c_id);
  //
  double aver_decay = 1.0 / ( getIterationCounterDbl() + 1.0 );
  AuxCoeffs(c_id) += aver_decay * ( Coeffs(c_id)-AuxCoeffs(c_id) );

}


}
}
