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

#include "tools/File.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_OPTIMIZER FES_DUMPER
/*
Tool to output biases and FESs from coefficients.

\par Examples

*/
//+ENDPLUMEDOC

class Opt_FesDumper : public Optimizer {

public:
  static void registerKeywords(Keywords&);
  explicit Opt_FesDumper(const ActionOptions&);
  void update() {};
  void coeffsUpdate(const unsigned int c_id = 0) {};
};


PLUMED_REGISTER_ACTION(Opt_FesDumper,"FES_DUMPER")


void Opt_FesDumper::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);
  keys.remove("COEFFS_FILE");
  keys.remove("COEFFS_OUTPUT");
  keys.remove("COEFFS_FMT");
  keys.remove("GRADIENT_FILE");
  keys.remove("GRADIENT_OUTPUT");
  keys.remove("GRADIENT_FMT");
  keys.remove("COEFFS_SET_ID_PREFIX");
  keys.remove("INITIAL_COEFFS");
  keys.remove("STRIDE");
  keys.remove("TARGETDIST_AVERAGES_FILE");
  keys.remove("TARGETDIST_AVERAGES_OUTPUT");
  keys.remove("TARGETDIST_AVERAGES_FMT");
  keys.remove("RESTART");
  keys.add("compulsory","COEFFS_INPUT","coeffs.data","the name of input coefficient file");
  //
}


Opt_FesDumper::Opt_FesDumper(const ActionOptions&ao):
  PLUMED_VES_OPTIMIZER_INIT(ao)
{
  turnOffHessian();
  turnOffCoeffsOutputFiles();

  plumed_massert(numberOfBiases()==1,"FES_DUMPER only works with one bias for now");
  plumed_massert(numberOfCoeffsSets()==1,"FES_DUMPER only works with one coefficient for now");
  std::string fname = "coeffs.data";
  parse("COEFFS_INPUT",fname);
  checkRead();
  IFile ifile;
  ifile.open(fname);
  while(ifile) {
    getBiasPntrs()[0]->resetBiasFileOutput();
    getBiasPntrs()[0]->resetFesFileOutput();

    getCoeffsPntrs()[0]->readOneSetFromFile(ifile);

    setIterationCounter(getCoeffsPntrs()[0]->getIterationCounter());
    if(isBiasOutputActive() && getIterationCounter()%getBiasOutputStride()==0) {
      writeBiasOutputFiles();
    }
    if(isFesOutputActive() && getIterationCounter()%getFesOutputStride()==0) {
      writeFesOutputFiles();
    }
    if(isFesProjOutputActive() && getIterationCounter()%getFesProjOutputStride()==0) {
      writeFesProjOutputFiles();
    }
  }
  log.printf("Stopping");
  plumed.stop();
}


}
}
