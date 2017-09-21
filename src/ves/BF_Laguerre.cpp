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

#include "BasisFunctions.h"

#include "core/ActionRegister.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_LAGUERRE
/*
Laguerre basis functions.

\par Examples

*/
//+ENDPLUMEDOC

class BF_Laguerre : public BasisFunctions {
  double scalingf_;
  virtual void setupLabels();
public:
  static void registerKeywords(Keywords&);
  explicit BF_Laguerre(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_Laguerre,"BF_LAGUERRE")


void BF_Laguerre::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","SCALING_FACTOR","scaling factor that is used to define the length scale of the basis functions. Depends also on the order employed. By default it is 1.0");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_Laguerre::BF_Laguerre(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  scalingf_(1.0)
{
  setNumberOfBasisFunctions(getOrder()+2);
  setIntrinsicInterval(intervalMin(),intervalMax());
  scalingf_ = 1.0;
  parse("SCALING_FACTOR",scalingf_);
  if(scalingf_!=1.0) {addKeywordToList("SCALING_FACTOR",scalingf_);}
  setNonPeriodic();
  setIntervalBounded();
  setType("Laguerre");
  setDescription("Laguerre functions");
  setupBF();
  checkRead();
}


void BF_Laguerre::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  argT = scalingf_*(argT-intervalMin());
  //
  std::vector<double> valuesL(getOrder()+1);
  std::vector<double> derivsL(getOrder()+1);
  //
  // calculate the Laguerre polynomials
  valuesL[0]=1.0;
  derivsL[0]=0.0;
  valuesL[1]=1.0-argT;
  derivsL[1]=-1.0;
  for(unsigned int i=1; i < getOrder(); i++) {
    double io = static_cast<double>(i);
    valuesL[i+1]  = ((2.0*io+1.0-argT)/(io+1.0))*valuesL[i] - (io/(io+1.0))*valuesL[i-1];
    derivsL[i+1]  = ((2.0*io+1.0-argT)/(io+1.0))*derivsL[i] - (1.0/(io+1.0))*valuesL[i] - (io/(io+1.0))*derivsL[i-1];
  }
  // calculate the Laguerre functions, the constant has index 0, the index is then shifted
  // index 1: exp(-x/2)*L0(x) = exp(-x/2), index 2: exp(-x/2)*L1(x), etc.
  values[0]=1.0;
  derivs[0]=0.0;
  double vexp = exp(-0.5*argT);
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    values[i] = vexp*valuesL[i-1];
    derivs[i] = scalingf_*vexp*(-0.5*valuesL[i-1]+derivsL[i-1]);
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


void BF_Laguerre::setupLabels() {
  setLabel(0,"1");
  for(unsigned int i=1; i < getNumberOfBasisFunctions() ; i++) {
    std::string is; Tools::convert(i-1,is);
    setLabel(i,"l"+is+"(s)");
  }
}


}
}
