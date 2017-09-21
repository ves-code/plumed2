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

#include "TargetDistribution.h"
#include "GridIntegrationWeights.h"

#include "core/ActionRegister.h"
#include "tools/Grid.h"

#include "GridProjWeights.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_TARGETDIST TD_MARGINAL_WELLTEMPERED
/*
One-dimensional marginal well-tempered target distribution (dynamic).

\par Examples

*/
//+ENDPLUMEDOC

class TD_MarginalWellTempered: public TargetDistribution {
private:
  double bias_factor_;
  std::vector<std::string> proj_args;
public:
  static void registerKeywords(Keywords&);
  explicit TD_MarginalWellTempered(const ActionOptions& ao);
  void updateGrid();
  double getValue(const std::vector<double>&) const;
  ~TD_MarginalWellTempered() {}
};


PLUMED_REGISTER_ACTION(TD_MarginalWellTempered,"TD_MARGINAL_WELLTEMPERED")


void TD_MarginalWellTempered::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","BIAS_FACTOR","The bias factor to be used for the well tempered distribution");
  keys.add("compulsory","PROJ_ARGS","The argument to be used for the marginal.");
}


TD_MarginalWellTempered::TD_MarginalWellTempered(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  bias_factor_(0.0),
  proj_args(0)
{
  parse("BIAS_FACTOR",bias_factor_);
  if(bias_factor_<=1.0) {
    plumed_merror(getName()+" target distribution: the value of the bias factor doesn't make sense, it should be larger than 1.0");
  }
  parseVector("MARGINAL_ARG",proj_args);
  if(proj_args.size()!=1) {
    plumed_merror(getName()+" target distribution: currently only supports one marginal argument in MARGINAL_ARG");
  }
  setDynamic();
  setFesGridNeeded();
  doNotAllowBiasCutoff();
  checkRead();
}


double TD_MarginalWellTempered::getValue(const std::vector<double>& argument) const {
  plumed_merror("getValue not implemented for TD_MarginalWellTempered");
  return 0.0;
}


void TD_MarginalWellTempered::updateGrid() {
  double beta_prime = getBeta()/bias_factor_;
  plumed_massert(getFesGridPntr()!=NULL,"the FES grid has to be linked to use TD_MarginalWellTempered!");
  //
  FesWeight* Fw = new FesWeight(getBeta());
  Grid fes_proj = getFesGridPntr()->project(proj_args,Fw);
  delete Fw;
  plumed_massert(fes_proj.getSize()==targetDistGrid().getSize(),"problem with FES projection - inconsistent grids");
  plumed_massert(fes_proj.getDimension()==1,"problem with FES projection - projected grid is not one-dimensional");
  //
  std::vector<double> integration_weights = GridIntegrationWeights::getIntegrationWeights(getTargetDistGridPntr());
  double norm = 0.0;
  for(Grid::index_t l=0; l<targetDistGrid().getSize(); l++) {
    double value = beta_prime * fes_proj.getValue(l);
    logTargetDistGrid().setValue(l,value);
    value = exp(-value);
    norm += integration_weights[l]*value;
    targetDistGrid().setValue(l,value);
  }
  targetDistGrid().scaleAllValuesAndDerivatives(1.0/norm);
  logTargetDistGrid().setMinToZero();
}

}
}
