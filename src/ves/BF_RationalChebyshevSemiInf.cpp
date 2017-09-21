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

//+PLUMEDOC VES_BASISF BF_CHEBYSHEV_RATIONAL_SEMI_INFINITE
/*
Rational Chebyshev basis functions on a semi-infinite interval

\par Examples

*/
//+ENDPLUMEDOC


class BF_RationalChebyshevSemiInf : public BasisFunctions {
  double mapf_;
public:
  static void registerKeywords(Keywords&);
  explicit BF_RationalChebyshevSemiInf(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_RationalChebyshevSemiInf,"BF_CHEBYSHEV_RATIONAL_SEMI_INFINITE")


void BF_RationalChebyshevSemiInf::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("compulsory","MAP_PARAMETER","Constant map parameter that defines the interval on which is the bias is active (should match the width of the function to be expanded).");
}

BF_RationalChebyshevSemiInf::BF_RationalChebyshevSemiInf(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  mapf_(0.0)
{
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  //
  parse("MAP_PARAMETER",mapf_); addKeywordToList("MAP_PARAMETER",mapf_);
  if(mapf_ <= 0.0) {plumed_merror("MAP_PARAMETER should be larger than 0");}
  //
  setNonPeriodic();
  setIntervalBounded();
  setType("rational-chebyshev-semi-inf");
  setDescription("Rational Chebyshev basis functions on a semi-infinite interval");
  setLabelPrefix("TL");
  setupBF();
  checkRead();
}


void BF_RationalChebyshevSemiInf::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  argT -= intervalMin();
  double derivf = (2.0*mapf_)/((argT+mapf_)*(argT+mapf_));
  argT = (argT-mapf_)/(argT+mapf_);
  //
  std::vector<double> derivsT(derivs.size());
  //
  values[0]=1.0;
  derivsT[0]=0.0;
  derivs[0]=0.0;
  values[1]=argT;
  derivsT[1]=1.0;
  derivs[1]=derivf;
  for(unsigned int i=1; i < getOrder(); i++) {
    values[i+1]  = 2.0*argT*values[i]-values[i-1];
    derivsT[i+1] = 2.0*values[i]+2.0*argT*derivsT[i]-derivsT[i-1];
    derivs[i+1]  = derivf*derivsT[i+1];
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}




}
}
