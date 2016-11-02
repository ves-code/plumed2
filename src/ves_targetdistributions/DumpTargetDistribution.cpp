/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2016 The ves-code team
   (see the PEOPLE-VES file at the root of the distribution for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of ves-code, version 1.

   ves-code is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ves-code is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with ves-code.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "TargetDistribution.h"
#include "TargetDistributionRegister.h"
#include "ves_tools/GridIntegrationWeights.h"

#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "core/PlumedMain.h"
#include "tools/File.h"
#include "tools/Grid.h"




// using namespace std;

namespace PLMD{
namespace generic{

//+PLUMEDOC FUNCTION DUMP_BASISFUNCTIONS
/*

*/
//+ENDPLUMEDOC


class DumpTargetDistribution :
  public Action
{
public:
  explicit DumpTargetDistribution(const ActionOptions&);
  void calculate(){}
  void apply(){}
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(DumpTargetDistribution,"DUMP_TARGET_DISTRIBUTION")

void DumpTargetDistribution::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  keys.add("compulsory","GRID_MIN","the lower bounds for the grid");
  keys.add("compulsory","GRID_MAX","the upper bounds for the grid");
  keys.add("compulsory","GRID_BINS","the number of bins used for the grid.");
  keys.add("optional","GRID_PERIODICITY","specfiy if the individual arguments should be made periodic (YES) or not (NO). By default all arguments are taken as not periodic.");
  keys.add("compulsory","TARGETDIST_FILE","filename of the file for writing the target distribution");
  keys.add("optional","LOG_TARGETDIST_FILE","filename of the file for writing the log of the target distribution");
  keys.add("compulsory","TARGET_DISTRIBUTION","the target distribution to be used.");
}

DumpTargetDistribution::DumpTargetDistribution(const ActionOptions&ao):
Action(ao)
{

  std::string targetdist_fname;
  parse("TARGETDIST_FILE",targetdist_fname);
  std::string log_targetdist_fname;
  parse("LOG_TARGETDIST_FILE",log_targetdist_fname);
  if(targetdist_fname==log_targetdist_fname){
    plumed_merror("error in DUMP_TARGET_DISTRIBUTION: TARGETDIST_FILE and LOG_TARGETDIST_FILE cannot be the same");
  }

  std::vector<unsigned int> grid_bins;
  parseVector("GRID_BINS",grid_bins);
  unsigned int nargs = grid_bins.size();

  std::vector<std::string> grid_min(nargs);
  parseVector("GRID_MIN",grid_min);
  std::vector<std::string> grid_max(nargs);
  parseVector("GRID_MAX",grid_max);

  std::vector<std::string> grid_periodicity(nargs);
  parseVector("GRID_PERIODICITY",grid_periodicity);
  if(grid_periodicity.size()==0){grid_periodicity.assign(nargs,"NO");}

  plumed_massert(grid_min.size()==nargs,"mismatch between number of values given for grid parameters");
  plumed_massert(grid_max.size()==nargs,"mismatch between number of values given for grid parameters");
  plumed_massert(grid_periodicity.size()==nargs,"mismatch between number of values given for grid parameters");

  std::string targetdist_keyword;
  parse("TARGET_DISTRIBUTION",targetdist_keyword);
  checkRead();
  //
  std::vector<Value*> arguments(nargs);
  for(unsigned int i=0; i < nargs; i++) {
    std::string is; Tools::convert(i+1,is);
    if(nargs==1){is="";}
    arguments[i]= new Value(NULL,"arg"+is,false);
    if(grid_periodicity[i]=="YES"){
      arguments[i]->setDomain(grid_min[i],grid_max[i]);
    }
    else if(grid_periodicity[i]=="NO"){
      arguments[i]->setNotPeriodic();
    }
    else{
      plumed_merror("wrong value given in GRID_PERIODICITY, either specfiy YES or NO");
    }
  }

  std::vector<std::string> words = Tools::getWords(targetdist_keyword);
  TargetDistribution* targetdist_pntr = targetDistributionRegister().create(words);
  if(targetdist_pntr->isDynamic()){
    plumed_merror("DUMP_TARGET_DISTRIBUTION only works for static target distributions");
  }
  targetdist_pntr->setupGrids(arguments,grid_min,grid_max,grid_bins);
  targetdist_pntr->update();
  Grid* targetdist_grid_pntr = targetdist_pntr->getTargetDistGridPntr();
  Grid* log_targetdist_grid_pntr = targetdist_pntr->getLogTargetDistGridPntr();

  double sum_grid = TargetDistribution::integrateGrid(targetdist_grid_pntr);
  log.printf("  target distribution integrated over the grid: %16.12f\n",sum_grid);
  log.printf("                                                (%30.16e)\n",sum_grid);
  //
  OFile ofile;
  ofile.link(*this);
  ofile.enforceBackup();
  ofile.open(targetdist_fname);
  targetdist_grid_pntr->writeToFile(ofile);
  ofile.close();
  if(log_targetdist_fname.size()>0){
    OFile ofile2;
    ofile2.link(*this);
    ofile2.enforceBackup();
    ofile2.open(log_targetdist_fname);
    log_targetdist_grid_pntr->writeToFile(ofile2);
    ofile2.close();
  }

  //
  delete targetdist_pntr;
  for(unsigned int i=0; i < nargs; i++) {
    delete arguments[i];
  }
  arguments.clear();


}





}
}
