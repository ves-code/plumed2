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
#include "LinearBasisSetExpansion.h"
#include "VesBias.h"
#include "ves_tools/CoeffsVector.h"
#include "ves_tools/VesTools.h"
#include "ves_basisfunctions/BasisFunctions.h"
#include "ves_targetdistributions/TargetDistribution.h"
#include "ves_targetdistributions/TargetDistributionRegister.h"


#include "tools/Keywords.h"
#include "tools/Grid.h"
#include "ves_tools/GridProjWeights.h"
#include "tools/Communicator.h"



namespace PLMD{
namespace bias{

//+PLUMEDOC VES LinearBasisSetExpansion
/*
*/
//+ENDPLUMEDOC

void LinearBasisSetExpansion::registerKeywords(Keywords& keys) {
}


LinearBasisSetExpansion::LinearBasisSetExpansion(
  const std::string& label,
  const double beta_in,
  Communicator& cc,
  std::vector<Value*> args_pntrs_in,
  std::vector<BasisFunctions*> basisf_pntrs_in,
  CoeffsVector* bias_coeffs_pntr_in):
label_(label),
action_pntr_(NULL),
vesbias_pntr_(NULL),
mycomm_(cc),
serial_(false),
beta_(beta_in),
kbt_(1.0/beta_),
args_pntrs_(args_pntrs_in),
nargs_(args_pntrs_.size()),
basisf_pntrs_(basisf_pntrs_in),
nbasisf_(basisf_pntrs_.size()),
bias_coeffs_pntr_(bias_coeffs_pntr_in),
ncoeffs_(0),
targetdist_averages_pntr_(NULL),
fes_wt_coeffs_pntr_(NULL),
bias_fes_scalingf_(-1.0),
log_ps_fes_contribution_(false),
welltemp_biasf_(-1.0),
inv_welltemp_biasf_(-1.0),
welltemp_beta_prime_(0.0),
bias_cutoff_active_(false),
grid_bins_(nargs_,100),
targetdist_grid_label_("targetdist"),
step_of_last_biasgrid_update(-1000),
step_of_last_fesgrid_update(-1000),
bias_grid_pntr_(NULL),
bias_withoutcutoff_grid_pntr_(NULL),
fes_grid_pntr_(NULL),
log_ps_grid_pntr_(NULL),
dynamic_ps_grid_pntr_(NULL)
{
  plumed_massert(args_pntrs_.size()==basisf_pntrs_.size(),"number of arguments and basis functions do not match");
  for(unsigned int k=0;k<nargs_;k++){nbasisf_[k]=basisf_pntrs_[k]->getNumberOfBasisFunctions();}
  //
  if(bias_coeffs_pntr_==NULL){
    bias_coeffs_pntr_ = new CoeffsVector(label_+".coeffs",args_pntrs_,basisf_pntrs_,mycomm_,true);
  }
  plumed_massert(bias_coeffs_pntr_->numberOfDimensions()==basisf_pntrs_.size(),"dimension of coeffs does not match with number of basis functions ");
  //
  ncoeffs_ = bias_coeffs_pntr_->numberOfCoeffs();
  targetdist_averages_pntr_ = new CoeffsVector(*bias_coeffs_pntr_);

  std::string targetdist_averages_label = bias_coeffs_pntr_->getLabel();
  if(targetdist_averages_label.find("coeffs")!=std::string::npos){
    targetdist_averages_label.replace(targetdist_averages_label.find("coeffs"), std::string("coeffs").length(), "targetdist_averages");
  }
  else {
    targetdist_averages_label += "_targetdist_averages";
  }
  targetdist_averages_pntr_->setLabels(targetdist_averages_label);
  //
}

LinearBasisSetExpansion::~LinearBasisSetExpansion() {
  if(targetdist_averages_pntr_!=NULL){
    delete targetdist_averages_pntr_;
  }
  if(fes_wt_coeffs_pntr_!=NULL){
    delete fes_wt_coeffs_pntr_;
  }
  if(bias_grid_pntr_!=NULL){
    delete bias_grid_pntr_;
  }
  if(bias_withoutcutoff_grid_pntr_!=NULL){
    delete bias_withoutcutoff_grid_pntr_;
  }
  if(fes_grid_pntr_!=NULL){
    delete fes_grid_pntr_;
  }
  if(log_ps_grid_pntr_!=NULL){
    delete log_ps_grid_pntr_;
  }
  if(dynamic_ps_grid_pntr_!=NULL){
    delete dynamic_ps_grid_pntr_;
  }
}


bool LinearBasisSetExpansion::isStaticTargetDistFileOutputActive() const {
  bool output_static_targetdist_files=true;
  if(vesbias_pntr_!=NULL){
    output_static_targetdist_files = vesbias_pntr_->isStaticTargetDistFileOutputActive();
  }
  return output_static_targetdist_files;
}


void LinearBasisSetExpansion::linkVesBias(bias::VesBias* vesbias_pntr_in) {
  vesbias_pntr_ = vesbias_pntr_in;
  action_pntr_ = static_cast<Action*>(vesbias_pntr_in);

}


void LinearBasisSetExpansion::linkAction(Action* action_pntr_in) {
  action_pntr_ = action_pntr_in;
}


void LinearBasisSetExpansion::setGridBins(const std::vector<unsigned int>& grid_bins_in) {
  plumed_massert(grid_bins_in.size()==nargs_,"the number of grid bins given doesn't match the number of arguments");
  grid_bins_=grid_bins_in;
}


void LinearBasisSetExpansion::setGridBins(const unsigned int nbins) {
  std::vector<unsigned int> grid_bins_in(nargs_,nbins);
  grid_bins_=grid_bins_in;
}


Grid* LinearBasisSetExpansion::setupGeneralGrid(const std::string& label_suffix, const std::vector<unsigned int>& nbins, const bool usederiv) {
  plumed_assert(nbins.size()==nargs_);
  std::vector<std::string> min(nargs_);
  std::vector<std::string> max(nargs_);
  for(unsigned int k=0;k<nargs_;k++){
    Tools::convert(basisf_pntrs_[k]->intervalMin(),min[k]);
    Tools::convert(basisf_pntrs_[k]->intervalMax(),max[k]);
  }
  Grid* grid_pntr = new Grid(label_+"."+label_suffix,args_pntrs_,min,max,nbins,false,usederiv);
  return grid_pntr;
}


Grid* LinearBasisSetExpansion::setupOneDimensionalMarginalGrid(const std::string& label_suffix, const unsigned int nbins, const unsigned int dim) {
  plumed_assert(dim<nargs_);
  std::vector<Value*> arg(1);
  arg[0]=args_pntrs_[dim];
  std::vector<std::string> min(1);
  std::vector<std::string> max(1);
  Tools::convert(basisf_pntrs_[dim]->intervalMin(),min[0]);
  Tools::convert(basisf_pntrs_[dim]->intervalMax(),max[0]);
  std::vector<unsigned int> nbinsv(1);
  nbinsv[0]=nbins;
  //
  Grid* grid_pntr = new Grid(label_+"."+label_suffix,arg,min,max,nbinsv,false,false);
  return grid_pntr;
}


void LinearBasisSetExpansion::setupBiasGrid(const bool usederiv) {
  plumed_massert(bias_grid_pntr_==NULL,"setupBiasGrid should only be called once: the bias grid has already been defined.");
  bias_grid_pntr_ = setupGeneralGrid("bias",grid_bins_,usederiv);
}


void LinearBasisSetExpansion::setupFesGrid() {
  plumed_massert(fes_grid_pntr_==NULL,"setupFesGrid should only be called once: the fes grid has already been defined.");
  if(bias_grid_pntr_==NULL){
    setupBiasGrid(false);
  }
  fes_grid_pntr_ = setupGeneralGrid("fes",grid_bins_,false);
}


void LinearBasisSetExpansion::setupFesProjGrid() {
  if(fes_grid_pntr_==NULL){
    setupFesGrid();
  }
}


void LinearBasisSetExpansion::updateBiasGrid() {
  plumed_massert(bias_grid_pntr_!=NULL,"the bias grid is not defined");
  if(getStepOfLastBiasGridUpdate() == action_pntr_->getStep()){
    return;
  }
  for(unsigned int l=0; l<bias_grid_pntr_->getSize(); l++){
    std::vector<double> forces(nargs_);
    std::vector<double> coeffsderivs_values(ncoeffs_);
    std::vector<double> args = bias_grid_pntr_->getPoint(l);
    double bias=getBiasAndForces(args,forces,coeffsderivs_values);
    //
    if(biasCutoffActive()){
      vesbias_pntr_->applyBiasCutoff(bias,forces);
    }
    //
    if(bias_grid_pntr_->hasDerivatives()){
      bias_grid_pntr_->setValueAndDerivatives(l,bias,forces);
    }
    else{
      bias_grid_pntr_->setValue(l,bias);
    }
    //
  }
  if(vesbias_pntr_!=NULL){
    vesbias_pntr_->setCurrentBiasMaxValue(bias_grid_pntr_->getMaxValue());
  }
  setStepOfLastBiasGridUpdate(action_pntr_->getStep());
}


void LinearBasisSetExpansion::updateBiasWithoutCutoffGrid() {
  plumed_massert(bias_withoutcutoff_grid_pntr_!=NULL,"the bias without cutoff grid is not defined");
  plumed_massert(biasCutoffActive(),"the bias cutoff has to be active");
  plumed_massert(vesbias_pntr_!=NULL,"has to be linked to a VesBias to work");
  //
  for(unsigned int l=0; l<bias_withoutcutoff_grid_pntr_->getSize(); l++){
    std::vector<double> forces(nargs_);
    std::vector<double> coeffsderivs_values(ncoeffs_);
    std::vector<double> args = bias_withoutcutoff_grid_pntr_->getPoint(l);
    double bias=getBiasAndForces(args,forces,coeffsderivs_values);
    if(bias_withoutcutoff_grid_pntr_->hasDerivatives()){
      bias_withoutcutoff_grid_pntr_->setValueAndDerivatives(l,bias,forces);
    }
    else{
      bias_withoutcutoff_grid_pntr_->setValue(l,bias);
    }
  }
  //
  double bias_max = bias_withoutcutoff_grid_pntr_->getMaxValue();
  double bias_min = bias_withoutcutoff_grid_pntr_->getMinValue();
  double shift = 0.0;
  bool bias_shifted=false;
  if(bias_min < 0.0){
    shift += -bias_min;
    bias_shifted=true;
    BiasCoeffs()[0] -= bias_min;
    bias_max -= bias_min;
  }
  if(bias_max > vesbias_pntr_->getBiasCutoffValue()){
    shift += -(bias_max-vesbias_pntr_->getBiasCutoffValue());
    bias_shifted=true;
    BiasCoeffs()[0] -= (bias_max-vesbias_pntr_->getBiasCutoffValue());
    bias_max -= (bias_max-vesbias_pntr_->getBiasCutoffValue());
  }
  if(bias_shifted){
    // this should be done inside a grid function really,
    // need to define my grid class for that
    for(unsigned int l=0; l<bias_withoutcutoff_grid_pntr_->getSize(); l++){
      double value = bias_withoutcutoff_grid_pntr_->getValue(l) + shift;
      bias_withoutcutoff_grid_pntr_->setValue(l,value);
    }
  }
  vesbias_pntr_->setCurrentBiasMaxValue(bias_max);
}


void LinearBasisSetExpansion::updateFesGrid() {
  plumed_massert(fes_grid_pntr_!=NULL,"the FES grid is not defined");
  updateBiasGrid();
  if(getStepOfLastFesGridUpdate() == action_pntr_->getStep()){
    return;
  }
  //
  double fes_scalingf = getBiasToFesScalingFactor();
  for(unsigned int l=0; l<fes_grid_pntr_->getSize(); l++){
    double fes_value = fes_scalingf*bias_grid_pntr_->getValue(l);
    if(log_ps_fes_contribution_){
      fes_value += log_ps_grid_pntr_->getValue(l);
    }
    fes_grid_pntr_->setValue(l,fes_value);
  }
  fes_grid_pntr_->setMinToZero();
  setStepOfLastFesGridUpdate(action_pntr_->getStep());
}


void LinearBasisSetExpansion::writeBiasGridToFile(const std::string& filepath, const bool append_file) const {
  plumed_massert(bias_grid_pntr_!=NULL,"the bias grid is not defined");
  OFile file;
  file.link(*action_pntr_);
  if(append_file){file.enforceRestart();}
  file.open(filepath);
  bias_grid_pntr_->writeToFile(file);
  file.close();
  //
  if(biasCutoffActive()){
    OFile file2;
    file2.link(*action_pntr_);
    if(append_file){file2.enforceRestart();}
    std::string filename = FileBase::appendSuffix(filepath,".without-cutoff");
    file2.open(filename);
    bias_withoutcutoff_grid_pntr_->writeToFile(file2);
    file2.close();
  }
}


void LinearBasisSetExpansion::writeFesGridToFile(const std::string& filepath, const bool append_file) const {
  plumed_massert(fes_grid_pntr_!=NULL,"the FES grid is not defined");
  OFile file;
  file.link(*action_pntr_);
  if(append_file){file.enforceRestart();}
  file.open(filepath);
  fes_grid_pntr_->writeToFile(file);
  file.close();
}


void LinearBasisSetExpansion::writeFesProjGridToFile(const std::vector<std::string>& proj_arg, const std::string& filepath, const bool append_file) const {
  plumed_massert(fes_grid_pntr_!=NULL,"the FES grid is not defined");
  FesWeight* Fw = new FesWeight(beta_);
  Grid proj_grid = fes_grid_pntr_->project(proj_arg,Fw);
  proj_grid.setMinToZero();
  OFile file;
  file.link(*action_pntr_);
  if(append_file){file.enforceRestart();}
  file.open(filepath);
  proj_grid.writeToFile(file);
  file.close();
  delete Fw;
}


double LinearBasisSetExpansion::getBiasAndForces(const std::vector<double>& args_values, std::vector<double>& forces, std::vector<double>& coeffsderivs_values, std::vector<BasisFunctions*>& basisf_pntrs_in, CoeffsVector* coeffs_pntr_in, Communicator* comm_in) {
  unsigned int nargs = args_values.size();
  plumed_assert(coeffs_pntr_in->numberOfDimensions()==nargs);
  plumed_assert(basisf_pntrs_in.size()==nargs);
  plumed_assert(forces.size()==nargs);
  plumed_assert(coeffsderivs_values.size()==coeffs_pntr_in->numberOfCoeffs());

  std::vector<double> args_values_trsfrm(nargs);
  std::vector<bool>   inside_interval(nargs,true);
  //
  std::vector< std::vector <double> > bf_values;
  std::vector< std::vector <double> > bf_derivs;
  //
  for(unsigned int k=0;k<nargs;k++){
    std::vector<double> tmp_val(basisf_pntrs_in[k]->getNumberOfBasisFunctions());
    std::vector<double> tmp_der(tmp_val.size());
    bool inside=true;
    basisf_pntrs_in[k]->getAllValues(args_values[k],args_values_trsfrm[k],inside,tmp_val,tmp_der);
    inside_interval[k]=inside;
    bf_values.push_back(tmp_val);
    bf_derivs.push_back(tmp_der);
    forces[k]=0.0;
  }
  //
  size_t stride=1;
  size_t rank=0;
  if(comm_in!=NULL)
  {
    stride=comm_in->Get_size();
    rank=comm_in->Get_rank();
  }
  // loop over coeffs
  double bias=0.0;
  for(size_t i=rank;i<coeffs_pntr_in->numberOfCoeffs();i+=stride){
    std::vector<unsigned int> indices=coeffs_pntr_in->getIndices(i);
    double coeff = coeffs_pntr_in->getValue(i);
    double bf_curr=1.0;
    for(unsigned int k=0;k<nargs;k++){
      bf_curr*=bf_values[k][indices[k]];
    }
    bias+=coeff*bf_curr;
    coeffsderivs_values[i] = bf_curr;
    for(unsigned int k=0;k<nargs;k++){
      forces[k]-=coeff*bf_curr*(bf_derivs[k][indices[k]]/bf_values[k][indices[k]]);
    }
  }
  //
  if(comm_in!=NULL){
    comm_in->Sum(bias);
    comm_in->Sum(forces);
  }
  return bias;
}


void LinearBasisSetExpansion::getBasisSetValues(const std::vector<double>& args_values, std::vector<double>& basisset_values, std::vector<BasisFunctions*>& basisf_pntrs_in, CoeffsVector* coeffs_pntr_in, Communicator* comm_in) {
  unsigned int nargs = args_values.size();
  plumed_assert(coeffs_pntr_in->numberOfDimensions()==nargs);
  plumed_assert(basisf_pntrs_in.size()==nargs);

  std::vector<double> args_values_trsfrm(nargs);
  std::vector< std::vector <double> > bf_values;
  //
  for(unsigned int k=0;k<nargs;k++){
    std::vector<double> tmp_val(basisf_pntrs_in[k]->getNumberOfBasisFunctions());
    std::vector<double> tmp_der(tmp_val.size());
    bool inside=true;
    basisf_pntrs_in[k]->getAllValues(args_values[k],args_values_trsfrm[k],inside,tmp_val,tmp_der);
    bf_values.push_back(tmp_val);
  }
  //
  size_t stride=1;
  size_t rank=0;
  if(comm_in!=NULL)
  {
    stride=comm_in->Get_size();
    rank=comm_in->Get_rank();
  }
  // loop over basis set
  for(size_t i=rank;i<coeffs_pntr_in->numberOfCoeffs();i+=stride){
    std::vector<unsigned int> indices=coeffs_pntr_in->getIndices(i);
    double bf_curr=1.0;
    for(unsigned int k=0;k<nargs;k++){
      bf_curr*=bf_values[k][indices[k]];
    }
    basisset_values[i] = bf_curr;
  }
  //
  if(comm_in!=NULL){
    comm_in->Sum(basisset_values);
  }
}


void LinearBasisSetExpansion::setupUniformTargetDistribution() {
  std::vector<TargetDistribution*> targetdist_dummy(nargs_,NULL);
  setupTargetDistribution(targetdist_dummy);
}


void LinearBasisSetExpansion::setupTargetDistribution(const std::vector<std::string>& targetdist_keywords) {
  if(targetdist_keywords.size()!=1 && targetdist_keywords.size()!=nargs_){
    plumed_merror("the number of target distribution keywords needs to be either 1 or equal to the number of arguments");
  }
  std::vector<TargetDistribution*> targetdist_pntrs(targetdist_keywords.size(),NULL);
  for(unsigned int k=0;k<targetdist_keywords.size();k++){
    std::vector<std::string> words = Tools::getWords(targetdist_keywords[k]);
    if(words[0]=="UNIFORM"){
      targetdist_pntrs[k] = NULL;
    }
    else{
      targetdist_pntrs[k] = targetDistributionRegister().create(words);
    }
  }
  setupTargetDistribution(targetdist_pntrs);
  for(unsigned int k=0;k<targetdist_pntrs.size();k++){
    if(targetdist_pntrs[k]!=NULL){
      delete targetdist_pntrs[k];
    }
  }
}


void LinearBasisSetExpansion::setupTargetDistribution(const std::vector<TargetDistribution*>& targetdist_pntrs) {
  if(targetdist_pntrs.size()!=1 && targetdist_pntrs.size()!=nargs_){
    plumed_merror("the number of target distribution pointers needs to be either 1 or equal to the number of arguments");
  }
  if(nargs_==1){
    setupOneDimensionalTargetDistribution(targetdist_pntrs);
  }
  else if(nargs_>1 && targetdist_pntrs.size()==1){
    setupNonSeperableTargetDistribution(targetdist_pntrs[0]);
  }
  else{
    setupSeperableTargetDistribution(targetdist_pntrs);
  }
}


void LinearBasisSetExpansion::setupSeperableTargetDistribution(const std::vector<TargetDistribution*>& targetdist_pntrs) {
  plumed_massert(targetdist_pntrs.size()==nargs_,"number of target distribution does not match the number of basis functions");
  plumed_massert(nargs_>1,"setupSeperableTargetDistribution should not be used for one-dimensional cases");
  //
  std::vector< std::vector <double> > bf_integrals(0);
  for(unsigned int k=0;k<nargs_;k++){
    if(targetdist_pntrs[k]!=NULL){
      bf_integrals.push_back(basisf_pntrs_[k]->getTargetDistributionIntegrals(targetdist_pntrs[k]));
    }
    else{
      bf_integrals.push_back(basisf_pntrs_[k]->getUniformIntegrals());
    }
  }
  //
  for(size_t i=0;i<ncoeffs_;i++){
    std::vector<unsigned int> indices=bias_coeffs_pntr_->getIndices(i);
    double value = 1.0;
    for(unsigned int k=0;k<nargs_;k++){
      value*=bf_integrals[k][indices[k]];
    }
    targetdist_averages_pntr_->setValue(i,value);
  }
  //
  bool all_targetdist_uniform = true;
  for(unsigned int k=0;k<nargs_;k++){
    if(targetdist_pntrs[k]!=NULL){
      all_targetdist_uniform = false;
      if(isStaticTargetDistFileOutputActive()){
        Grid* grid1d_pntr_tmp = setupOneDimensionalMarginalGrid(targetdist_grid_label_+"_proj",grid_bins_[k],k);
        targetdist_pntrs[k]->calculateDistributionOnGrid(grid1d_pntr_tmp);
        writeTargetDistGridToFile(grid1d_pntr_tmp,"proj-" + args_pntrs_[k]->getName());
        delete grid1d_pntr_tmp;
      }
    }
  }
  //
  if(!all_targetdist_uniform){
    std::vector<TargetDistribution*> targetdist_pntrs_modifed(nargs_);
    for(unsigned int k=0;k<nargs_;k++){
      targetdist_pntrs_modifed[k]=targetdist_pntrs[k];
      if(targetdist_pntrs[k]==NULL){
        std::vector<std::string> words;
        std::string tmp1;
        words.push_back("UNIFORM");
        VesTools::convertDbl2Str(basisf_pntrs_[k]->intervalMin(),tmp1);
        words.push_back("MINIMA=" + tmp1);
        VesTools::convertDbl2Str(basisf_pntrs_[k]->intervalMax(),tmp1);
        words.push_back("MAXIMA=" + tmp1);
        targetdist_pntrs_modifed[k] = targetDistributionRegister().create(words);
      }
    }
    Grid* ps_grid_pntr = setupGeneralGrid(targetdist_grid_label_,grid_bins_,false);
    TargetDistribution::calculateSeperableDistributionOnGrid(ps_grid_pntr,targetdist_pntrs_modifed);
    if(isStaticTargetDistFileOutputActive()){
      writeTargetDistGridToFile(ps_grid_pntr);
    }
    //
    log_ps_grid_pntr_ = setupGeneralGrid(targetdist_grid_label_+"log-beta",grid_bins_,false);
    VesTools::copyGridValues(ps_grid_pntr,log_ps_grid_pntr_);
    delete ps_grid_pntr;
    log_ps_grid_pntr_->logAllValuesAndDerivatives( (-1.0/beta_) );
    log_ps_fes_contribution_ = true;
    if(isStaticTargetDistFileOutputActive()){
      writeTargetDistGridToFile(log_ps_grid_pntr_,"log-beta");
    }
    //
    for(unsigned int k=0;k<nargs_;k++){
      if(targetdist_pntrs[k]==NULL && targetdist_pntrs_modifed[k]!=NULL){
        delete targetdist_pntrs_modifed[k];
      }
    }
  }
}


void LinearBasisSetExpansion::setupOneDimensionalTargetDistribution(const std::vector<TargetDistribution*>& targetdist_pntrs) {
  plumed_massert(nargs_==1,"setupOneDimensionalTargetDistribution should only be used for one-dimensional cases");
  plumed_massert(targetdist_pntrs.size()==nargs_,"number of target distribution does not match the number of basis functions");
  //
  std::vector<double> bf_integrals;
  if(targetdist_pntrs[0]!=NULL){
    bf_integrals = basisf_pntrs_[0]->getTargetDistributionIntegrals(targetdist_pntrs[0]);
  }
  else{
    bf_integrals = basisf_pntrs_[0]->getUniformIntegrals();
  }
  plumed_massert(bf_integrals.size()==ncoeffs_,"something wrong in setupOneDimensionalTargetDistribution");
  targetdist_averages_pntr_->setValues(bf_integrals);
  //
  if(targetdist_pntrs[0]!=NULL){
    Grid* ps_grid_pntr = setupGeneralGrid(targetdist_grid_label_,grid_bins_,false);
    targetdist_pntrs[0]->calculateDistributionOnGrid(ps_grid_pntr);
    if(isStaticTargetDistFileOutputActive()){
      writeTargetDistGridToFile(ps_grid_pntr);
    }
    //
    log_ps_grid_pntr_ = setupGeneralGrid(targetdist_grid_label_+"log-beta",grid_bins_,false);
    VesTools::copyGridValues(ps_grid_pntr,log_ps_grid_pntr_);
    delete ps_grid_pntr;
    log_ps_grid_pntr_->logAllValuesAndDerivatives( (-1.0/beta_) );
    log_ps_fes_contribution_ = true;
    if(isStaticTargetDistFileOutputActive()){
      writeTargetDistGridToFile(log_ps_grid_pntr_,"log-beta");
    }
  }
}


void LinearBasisSetExpansion::setupNonSeperableTargetDistribution(const TargetDistribution* targetdist_pntr) {
  if(targetdist_pntr==NULL){
    setupUniformTargetDistribution();
    return;
  }
  plumed_massert(nargs_>1,"setupNonSeperableTargetDistribution should not be used for one-dimensional cases");
  Grid* ps_grid_pntr = setupGeneralGrid(targetdist_grid_label_,grid_bins_,false);
  targetdist_pntr->calculateDistributionOnGrid(ps_grid_pntr);
  calculateTargetDistAveragesFromGrid(ps_grid_pntr);
  if(isStaticTargetDistFileOutputActive()){
    writeTargetDistGridToFile(ps_grid_pntr);
  }
  log_ps_grid_pntr_ = setupGeneralGrid(targetdist_grid_label_+"_log-beta",grid_bins_,false);
  VesTools::copyGridValues(ps_grid_pntr,log_ps_grid_pntr_);
  log_ps_grid_pntr_->logAllValuesAndDerivatives( (-1.0/beta_) );
  log_ps_fes_contribution_ = true;
  if(isStaticTargetDistFileOutputActive()){
    writeTargetDistGridToFile(log_ps_grid_pntr_,"log-beta");
    for(unsigned int k=0; k<args_pntrs_.size(); k++){
      Grid proj_grid = TargetDistribution::getMarginalGrid(ps_grid_pntr,args_pntrs_[k]->getName());
      writeTargetDistGridToFile(&proj_grid,"proj-"+args_pntrs_[k]->getName());
    }
  }
  delete ps_grid_pntr;
}


void LinearBasisSetExpansion::setupWellTemperedTargetDistribution(const double biasf) {
  plumed_massert(welltemp_biasf_<0.0,"setupWellTemperedTargetDistribution should only be run once.");
  plumed_massert(biasf>1.0,"setupWellTemperedTargetDistribution: the value of the bias factor doesn't make sense, it should be larger than 1.0");
  plumed_massert(fes_wt_coeffs_pntr_==NULL,"setupWellTemperedTargetDistribution should only be called once: the CoeffsVector for the well-tempered FES has already been defined");
  plumed_massert(dynamic_ps_grid_pntr_==NULL,"setupWellTemperedTargetDistribution should only be called once: the grid for the well-tempered p(s) has already been defined");
  //
  welltemp_biasf_=biasf;
  inv_welltemp_biasf_ = 1.0/welltemp_biasf_;
  welltemp_beta_prime_ = beta_/welltemp_biasf_;
  bias_fes_scalingf_ = -(welltemp_biasf_)/(welltemp_biasf_-1.0);
  //
  fes_wt_coeffs_pntr_ = new CoeffsVector(*bias_coeffs_pntr_);
  std::string fes_wt_label = bias_coeffs_pntr_->getLabel();
  if(fes_wt_label.find("coeffs")!=std::string::npos){
    fes_wt_label.replace(fes_wt_label.find("coeffs"), std::string("coeffs").length(), "fes_wt_coeffs");
  }
  else {
    fes_wt_label += "_fes_wt";
  }
  fes_wt_coeffs_pntr_->setLabels(fes_wt_label);
  //
  dynamic_ps_grid_pntr_ = setupGeneralGrid("ps_wt",grid_bins_,false);
}

void LinearBasisSetExpansion::setWellTemperedBiasFactor(const double biasf) {
  plumed_massert(welltemp_biasf_>1.0,"setupWellTemperedTargetDistribution has to be run before changing the bias factor using setWellTemperedBiasFactor.");
  plumed_massert(biasf>1.0,"setWellTemperedBiasFactor: the value of the bias factor doesn't make sense, it should be larger than 1.0");
  welltemp_biasf_=biasf;
  inv_welltemp_biasf_ = 1.0/welltemp_biasf_;
  welltemp_beta_prime_ = beta_/welltemp_biasf_;
  bias_fes_scalingf_ = -(welltemp_biasf_)/(welltemp_biasf_-1.0);
}


void LinearBasisSetExpansion::updateWellTemperedPsGrid() {
  double norm = 0.0;
  size_t stride=mycomm_.Get_size();
  size_t rank=mycomm_.Get_rank();
  for(unsigned int l=rank; l<dynamic_ps_grid_pntr_->getSize(); l+=stride){
    std::vector<double> args = dynamic_ps_grid_pntr_->getPoint(l);
    double value = -welltemp_beta_prime_ * getFES_WellTempered(args,false);
    value = exp(value);
    norm += value;
    dynamic_ps_grid_pntr_->setValue(l,value);
  }
  mycomm_.Sum(norm);
  norm = 1.0/(dynamic_ps_grid_pntr_->getBinVolume()*norm);
  dynamic_ps_grid_pntr_->scaleAllValuesAndDerivatives(norm);
  dynamic_ps_grid_pntr_->mpiSumValuesAndDerivatives(mycomm_);
}


void LinearBasisSetExpansion::updateWellTemperedTargetDistribution() {
  plumed_assert(welltemp_biasf_>1.0);
  // Update FES WT coeffs
  FesWTCoeffs() = -BiasCoeffs() + inv_welltemp_biasf_*FesWTCoeffs();
  // Update p(s) grid
  updateWellTemperedPsGrid();
  // calcuate coeffs derivs from grid
  calculateTargetDistAveragesFromGrid(dynamic_ps_grid_pntr_);
}


void LinearBasisSetExpansion::calculateTargetDistAveragesFromGrid(const Grid* ps_grid_pntr, const bool normalize_dist) {
  plumed_assert(ps_grid_pntr!=NULL);
  std::vector<double> targetdist_averages(ncoeffs_,0.0);
  double sum_grid = 0.0;
  double binVol = ps_grid_pntr->getBinVolume();
  size_t stride=mycomm_.Get_size();
  size_t rank=mycomm_.Get_rank();
  for(unsigned int l=rank; l<ps_grid_pntr->getSize(); l+=stride){
    std::vector<double> args_values = ps_grid_pntr->getPoint(l);
    std::vector<double> basisset_values(ncoeffs_);
    getBasisSetValues(args_values,basisset_values,false);
    double weight = ps_grid_pntr->getValue(l)*binVol;
    sum_grid += weight;
    for(unsigned int i=0; i<ncoeffs_; i++){
      targetdist_averages[i] += weight*basisset_values[i];
    }
  }
  if(normalize_dist){
    mycomm_.Sum(sum_grid);
    for(unsigned int i=0; i<ncoeffs_; i++){
      targetdist_averages[i] /= sum_grid;
    }
  }
  mycomm_.Sum(targetdist_averages);
  // the overall constant;
  targetdist_averages[0] = 1.0;
  TargetDistAverages() = targetdist_averages;
}


void LinearBasisSetExpansion::updateBiasCutoffPsGrid() {
  plumed_massert(vesbias_pntr_!=NULL,"has to be linked to a VesBias to work");
  plumed_massert(bias_withoutcutoff_grid_pntr_!=NULL,"the bias grid has to be defined");
  plumed_massert(dynamic_ps_grid_pntr_!=NULL,"the p(s) grid has to be defined");
  double norm = 0.0;
  size_t stride=mycomm_.Get_size();
  size_t rank=mycomm_.Get_rank();
  for(unsigned int l=rank; l<dynamic_ps_grid_pntr_->getSize(); l+=stride){
    std::vector<double> args = dynamic_ps_grid_pntr_->getPoint(l);
    double bias = bias_withoutcutoff_grid_pntr_->getValue(l);
    double deriv_factor_sf = 0.0;
    // this comes from the p(s)
    double value = vesbias_pntr_->getBiasCutoffSwitchingFunction(bias,deriv_factor_sf);
    norm += value;
    // this comes from the derivative of V(s)
    value *= deriv_factor_sf;
    dynamic_ps_grid_pntr_->setValue(l,value);
  }
  mycomm_.Sum(norm);
  norm = 1.0/(dynamic_ps_grid_pntr_->getBinVolume()*norm);
  dynamic_ps_grid_pntr_->scaleAllValuesAndDerivatives(norm);
  dynamic_ps_grid_pntr_->mpiSumValuesAndDerivatives(mycomm_);
}


void LinearBasisSetExpansion::updateBiasCutoffTargetDistribution() {
  plumed_massert(biasCutoffActive(),"the bias cutoff is not active");
  updateBiasWithoutCutoffGrid();
  updateBiasCutoffPsGrid();
  calculateTargetDistAveragesFromGrid(dynamic_ps_grid_pntr_);
}


void LinearBasisSetExpansion::setupBiasCutoffTargetDistribution() {
  plumed_massert(dynamic_ps_grid_pntr_==NULL,"setupBiasCutoffTargetDistribution should only be called once: the grid for the p(s) has already been defined");
  plumed_massert(bias_withoutcutoff_grid_pntr_==NULL,"setupBiasGrid should only be called once: the bias with cutoff grid has already been defined.");
  plumed_massert(vesbias_pntr_!=NULL,"has to be linked to a VesBias to work");
  plumed_massert(!bias_cutoff_active_,"setupBiasCutoffTargetDistribution should only be called once");
  //
  bias_cutoff_active_=true;
  dynamic_ps_grid_pntr_ = setupGeneralGrid("ps_cutoff",grid_bins_,false);
  bias_withoutcutoff_grid_pntr_ = setupGeneralGrid("bias_withoutcutoff",grid_bins_,true);
}


void LinearBasisSetExpansion::writeTargetDistGridToFile(Grid* grid_pntr, const std::string& suffix) const {
  plumed_massert(grid_pntr!=NULL,"the grid is not defined");
  //
  std::string filepath = "targetdist.data";
  if(suffix.size()>0){
    filepath = FileBase::appendSuffix(filepath,"."+suffix);
  }
  if(vesbias_pntr_!=NULL){
    filepath = vesbias_pntr_->getCurrentTargetDistOutputFilename(suffix);
  }
  OFile file;
  file.link(*action_pntr_);
  // file.setBackupString("bck");
  file.open(filepath);
  grid_pntr->writeToFile(file);
  file.close();
}


void LinearBasisSetExpansion::writeDynamicTargetDistGridToFile(const bool do_projections, const std::string& suffix) const {
  writeTargetDistGridToFile(dynamic_ps_grid_pntr_);
  if(do_projections){
    for(unsigned int k=0; k<args_pntrs_.size(); k++){
      Grid proj_grid = TargetDistribution::getMarginalGrid(dynamic_ps_grid_pntr_,args_pntrs_[k]->getName());
      writeTargetDistGridToFile(&proj_grid,"proj-"+args_pntrs_[k]->getName());
    }
  }
}











}

}
