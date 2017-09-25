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

#include "CoeffsVector.h"
#include "VesTools.h"
#include "VesBias.h"


#include "tools/File.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_UTILS VES_OUTPUT_FES
/*
Tool to output biases and FESs from previously obtained coefficients.

\par Examples

*/
//+ENDPLUMEDOC

class Opt_FesDumper : public Action {

public:
  static void registerKeywords(Keywords&);
  explicit Opt_FesDumper(const ActionOptions&);
  void update() {}
  void calculate() {}
  void apply() {}
};


PLUMED_REGISTER_ACTION(Opt_FesDumper,"VES_OUTPUT_FES")


void Opt_FesDumper::registerKeywords(Keywords& keys) {
  keys.add("compulsory","BIAS","the label of the VES bias for to output the FESs and the bias files");
  keys.add("compulsory","COEFFS_INPUT","the name of input coefficient file");
  keys.add("optional","BIAS_OUTPUT","how often the bias(es) should be written out to file. Note that the value is given in terms of coefficent iterations.");
  keys.add("optional","FES_OUTPUT","how often the FES(s) should be written out to file. Note that the value is given in terms of coefficent iterations.");
  keys.add("optional","FES_PROJ_OUTPUT","how often the projections of the FES(s) should be written out to file. Note that the value is given in terms of coefficent iterations.");
  //
}


Opt_FesDumper::Opt_FesDumper(const ActionOptions&ao):
  Action(ao)
{

  std::vector<std::string> bias_labels;
  parseVector("BIAS",bias_labels);
  if(bias_labels.size()>1) {
    plumed_merror(getName()+" only support one VES bias");
  }

  std::string error_msg = "";
  std::vector<VesBias*> bias_pntrs = VesTools::getPointersFromLabels<VesBias*>(bias_labels,plumed.getActionSet(),error_msg);
  if(error_msg.size()>0) {plumed_merror("Error in keyword BIAS of "+getName()+": "+error_msg);}

  for(unsigned int i=0; i<bias_pntrs.size(); i++) {
    if(bias_pntrs[i]->numberOfCoeffsSets()>1) {
      plumed_merror(getName()+" at the moment supports only VES biases with a single coefficient set");
    }
  }

  std::vector<std::string> coeffs_fnames;
  parseVector("COEFFS_INPUT",coeffs_fnames);
  if(coeffs_fnames.size()!=bias_pntrs.size()) {
    plumed_merror(getName()+": there have to be as many coefficient file given in COEFFS_INPUT as VES biases given in BIAS");
  }

  unsigned int bias_output_stride = 0;
  parse("BIAS_OUTPUT",bias_output_stride);

  unsigned int fes_output_stride = 0;
  parse("FES_OUTPUT",fes_output_stride);

  unsigned int fesproj_output_stride = 0;
  parse("FES_PROJ_OUTPUT",fesproj_output_stride);

  if(bias_output_stride == 0 && fes_output_stride == 0 && fesproj_output_stride == 0) {
    plumed_merror(getName()+": you are not telling the action to do anything, you need to use one of the keywords BIAS_OUTPUT, FES_OUTPUT, or FES_PROJ_OUTPUT");
  }

  checkRead();

  for(unsigned int i=0; i<bias_pntrs.size(); i++) {

    if(bias_pntrs[i]->dynamicTargetDistribution()) {
      plumed_merror(getName()+" does not support dynamic target distributions at the moment");
    }

    if(bias_output_stride>0) {
      bias_pntrs[i]->enableBiasFileOutput();
      bias_pntrs[i]->setupBiasFileOutput();
    }

    if(fes_output_stride>0) {
      bias_pntrs[i]->enableFesFileOutput();
      bias_pntrs[i]->setupFesFileOutput();
    }

    if(fesproj_output_stride>0) {
      bias_pntrs[i]->enableFesProjFileOutput();
      bias_pntrs[i]->setupFesProjFileOutput();
    }

    bias_pntrs[i]->enableIterationNumberInFilenames();

    IFile ifile;
    ifile.open(coeffs_fnames[i]);

    while(ifile) {

      bias_pntrs[i]->resetBiasFileOutput();
      bias_pntrs[i]->resetFesFileOutput();

      bias_pntrs[i]->getCoeffsPntrs()[0]->readOneSetFromFile(ifile);
      unsigned int iteration = bias_pntrs[i]->getCoeffsPntrs()[0]->getIterationCounter();

      if(bias_output_stride>0 && iteration%bias_output_stride==0) {
        bias_pntrs[i]->writeBiasToFile();
      }

      if(fes_output_stride>0 && iteration%fes_output_stride==0) {
        bias_pntrs[i]->writeFesToFile();
      }

      if(fesproj_output_stride>0 && iteration%fesproj_output_stride==0) {
        bias_pntrs[i]->writeFesProjToFile();
      }

    }

  }

  log.printf("Stopping");
  plumed.stop();
}


}
}
