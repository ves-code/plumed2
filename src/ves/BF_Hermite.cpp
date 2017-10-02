/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "BasisFunctions.h"

#include "core/ActionRegister.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_HERMITE
/*
Hermite basis functions.

\par Examples

*/
//+ENDPLUMEDOC

class BF_Hermite : public BasisFunctions {
  double scalingf_;
  double center_;
  std::vector<double> normf_;
  virtual void setupLabels();
public:
  static void registerKeywords(Keywords&);
  explicit BF_Hermite(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_Hermite,"BF_HERMITE")


void BF_Hermite::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","SCALING_FACTOR","scaling factor that is used to define the length scale of the basis functions. Depends also on the order employed. By default it is 1.0");
  keys.add("optional","CENTER","the location of the center of the Hermite functions. By default it is 0.0");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_Hermite::BF_Hermite(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  scalingf_(1.0),
  center_(0.0),
  normf_(0)
{
  setNumberOfBasisFunctions(getOrder()+2);
  setIntrinsicInterval(intervalMin(),intervalMax());
  scalingf_=1.0;
  parse("SCALING_FACTOR",scalingf_);
  if(scalingf_!=1.0) {addKeywordToList("SCALING_FACTOR",scalingf_);}
  center_=0.0;
  parse("CENTER",center_);
  if(center_!=0.0) {addKeywordToList("CENTER",center_);}
  //
  // To normalize with 1.0/sqrt(sqrt(pi)*2^n*n!)
  normf_.resize(getOrder()+1);
  normf_[0] = 1.0/(sqrt( sqrt(pi) ));
  for(unsigned int i=1; i<getOrder()+1; i++) {
    double io = static_cast<double>(i);
    normf_[i] = normf_[i-1]*(1.0/sqrt(io*2.0));
  }
  //
  setNonPeriodic();
  setIntervalBounded();
  setType("Hermite");
  setDescription("Hermite functions");
  setupBF();
  checkRead();
}


void BF_Hermite::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  argT = scalingf_*(argT-center_);
  //
  std::vector<double> valuesH(getOrder()+1);
  std::vector<double> derivsH(getOrder()+1);
  //
  // calculate the Hermite polynomials
  valuesH[0]=1.0;
  derivsH[0]=0.0;
  valuesH[1]=2.0*argT;
  derivsH[1]=2.0;
  for(unsigned int i=1; i < getOrder(); i++) {
    double io = static_cast<double>(i);
    valuesH[i+1]  = 2.0*argT*valuesH[i] - 2.0*io*valuesH[i-1];
    derivsH[i+1]  = 2.0*argT*derivsH[i] + 2.0*valuesH[i] - 2.0*io*derivsH[i-1];
  }
  // calculate the Hermite functions, the constant has index 0, the index is then shifted
  // index 1: exp(-x^2/2)*H0(x) = exp(-x^2/2), index 2: exp(-x^2/2)*H1(x), etc.
  values[0]=1.0;
  derivs[0]=0.0;
  double vexp = exp(-0.5*argT*argT);
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    values[i] = normf_[i-1] * vexp*valuesH[i-1];
    derivs[i] = normf_[i-1] * scalingf_*vexp*(-argT*valuesH[i-1]+derivsH[i-1]);
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


void BF_Hermite::setupLabels() {
  setLabel(0,"1");
  for(unsigned int i=1; i < getNumberOfBasisFunctions() ; i++) {
    std::string is; Tools::convert(i-1,is);
    setLabel(i,"h"+is+"(s)");
  }
}


}
}
