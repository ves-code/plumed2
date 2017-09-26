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
#include "CoeffsMatrix.h"

#include "core/ActionRegister.h"
#include "core/PlumedMain.h"



namespace PLMD {
namespace ves {

//+PLUMEDOC VES_OPTIMIZER OPT_AVERAGED_SGD
/*
Averaged stochastic gradient decent with fixed step size.

\par Algorithim

This optimizer updates the coefficients according to the averaged stochastic gradient decent algorithm described in ref \cite Bach-NIPS-2013. This algorithm considers two sets of coefficients, the so-called instantaneous coefficients that are updated according to the recursion formula given by
\f[
\boldsymbol{\alpha}^{(n+1)} = \boldsymbol{\alpha}^{(n)} -
\mu \left[
\nabla \Omega(\bar{\boldsymbol{\alpha}}^{(n)}) +
\mathbf{H}(\bar{\boldsymbol{\alpha}}^{(n)})
[\boldsymbol{\alpha}^{(n)}-\bar{\boldsymbol{\alpha}}^{(n)}]
\right],
\f]
where \f$\mu\f$ is a fixed step size and the gradient \f$ \nabla\Omega(\bar{\boldsymbol{\alpha}}^{(n)})\f$ and the Hessian \f$\mathbf{H}(\bar{\boldsymbol{\alpha}}^{(n)})\f$ depend on the averaged coefficients defined as
\f[
\bar{\boldsymbol{\alpha}}^{(n)} = \frac{1}{n+1} \sum_{k=0}^{n} \boldsymbol{\alpha}^{(k)}.
\f]
This means that the bias acting on the system depends on the averaged coefficients \f$\bar{\boldsymbol{\alpha}}^{(n)}\f$ which leads to a smooth convergence of the bias and the estimated free energy surface. Furthermore, this allows for a rather short sampling time for each iteration, for classical MD simulations typical sampling times are on the order of few ps (around 1000-4000 MD steps).

Currently it is only supported to employ the diagonal part of the Hessian which is generally sufficient. Support for employing the full Hessian will be added later on.

The VES bias that is to be optimized should be specified using the
BIAS keyword.
The fixed step size \f$\mu\f$ is given using the STEPSIZE keyword.
The frequency of updating the coefficients is given using the
STRIDE keyword where the value is given in the number of MD steps.
For example, if the MD time step is 0.02 ps and STRIDE=2000 will the
coefficients be updated every 4 ps.
The coefficients will be outputted to the file given by the
COEFFS_FILE keyword. How often the coefficients are written
to this file is controlled by the COEFFS_OUTPUT keyword.


\par Multiple walkers

This optimizer supports the usage of multiple walkers where different copies of the system share the same bias potential (i.e. coefficients) and cooperatively sample the averages needed for the gradient and Hessian. This can significantly help with convergence in difficult cases. It is of course best to start the different copies from different positions in CV space. To activate this option you just need to add the MULTIPLE_WALKERS flag. Note that this is only supported if the MD code support running multiple replicas connected via MPI.



\par Examples

*/
//+ENDPLUMEDOC

class Opt_BachAveragedSGD : public Optimizer {
private:
  std::vector<CoeffsVector*> combinedgradient_pntrs_;
  unsigned int combinedgradient_wstride_;
  std::vector<OFile*> combinedgradientOFiles_;
  double decaying_aver_tau_;
private:
  CoeffsVector& CombinedGradient(const unsigned int c_id) const {return *combinedgradient_pntrs_[c_id];}
  double getAverDecay() const;
public:
  static void registerKeywords(Keywords&);
  explicit Opt_BachAveragedSGD(const ActionOptions&);
  ~Opt_BachAveragedSGD();
  void coeffsUpdate(const unsigned int c_id = 0);
};


PLUMED_REGISTER_ACTION(Opt_BachAveragedSGD,"OPT_AVERAGED_SGD")


void Opt_BachAveragedSGD::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);
  Optimizer::useFixedStepSizeKeywords(keys);
  Optimizer::useMultipleWalkersKeywords(keys);
  Optimizer::useHessianKeywords(keys);
  Optimizer::useMaskKeywords(keys);
  Optimizer::useRestartKeywords(keys);
  // Optimizer::useMonitorAveragesKeywords(keys);
  Optimizer::useDynamicTargetDistributionKeywords(keys);
  keys.add("hidden","COMBINED_GRADIENT_FILE","the name of output file for the combined gradient (gradient + Hessian term)");
  keys.add("hidden","COMBINED_GRADIENT_OUTPUT","how often the combined gradient should be written to file. This parameter is given as the number of bias iterations. It is by default 100 if COMBINED_GRADIENT_FILE is specficed");
  keys.add("hidden","COMBINED_GRADIENT_FMT","specify format for combined gradient file(s) (useful for decrease the number of digits in regtests)");
  keys.add("optional","EXP_DECAYING_AVER","calculate the averaged coefficients using exponentially decaying averaging using the decaying constant given here in the number of iterations");
}


Opt_BachAveragedSGD::~Opt_BachAveragedSGD() {
  for(unsigned int i=0; i<combinedgradient_pntrs_.size(); i++) {
    delete combinedgradient_pntrs_[i];
  }
  for(unsigned int i=0; i<combinedgradientOFiles_.size(); i++) {
    combinedgradientOFiles_[i]->close();
    delete combinedgradientOFiles_[i];
  }
}


Opt_BachAveragedSGD::Opt_BachAveragedSGD(const ActionOptions&ao):
  PLUMED_VES_OPTIMIZER_INIT(ao),
  combinedgradient_pntrs_(0),
  combinedgradient_wstride_(100),
  combinedgradientOFiles_(0),
  decaying_aver_tau_(0.0)
{
  log.printf("  Averaged stochastic gradient decent, see and cite ");
  log << plumed.cite("Bach and Moulines, NIPS 26, 773-781 (2013)");
  log.printf("\n");
  unsigned int decaying_aver_tau_int=0;
  parse("EXP_DECAYING_AVER",decaying_aver_tau_int);
  if(decaying_aver_tau_int>0) {
    decaying_aver_tau_ = static_cast<double>(decaying_aver_tau_int);
    log.printf("  Coefficients calculated using an exponentially decaying average with a decaying constant of %u iterations, see and cite ",decaying_aver_tau_int);
    log << plumed.cite("Invernizzi, Valsson, and Parrinello, Proc. Natl. Acad. Sci. USA 114, 3370-3374 (2017)");
    log.printf("\n");
  }
  //
  std::vector<std::string> combinedgradient_fnames;
  parseFilenames("COMBINED_GRADIENT_FILE",combinedgradient_fnames);
  parse("COMBINED_GRADIENT_OUTPUT",combinedgradient_wstride_);
  setupOFiles(combinedgradient_fnames,combinedgradientOFiles_,useMultipleWalkers());
  std::string combinedgradient_fmt="";
  parse("COMBINED_GRADIENT_FMT",combinedgradient_fmt);
  if(combinedgradient_fnames.size()>0) {
    for(unsigned int i=0; i<numberOfCoeffsSets(); i++) {
      CoeffsVector* combinedgradient_tmp = new CoeffsVector(*getGradientPntrs()[i]);
      std::string label = getGradientPntrs()[i]->getLabel();
      if(label.find("gradient")!=std::string::npos) {
        label.replace(label.find("gradient"), std::string("gradient").length(), "combined_gradient");
      }
      else {
        label += "_combined";
      }
      combinedgradient_tmp->setLabels(label);
      if(combinedgradient_fmt.size()>0) {
        combinedgradient_tmp->setOutputFmt(combinedgradient_fmt);
      }
      combinedgradient_pntrs_.push_back(combinedgradient_tmp);
    }
    //
    if(numberOfCoeffsSets()==1) {
      log.printf("  Combined gradient (gradient + Hessian term) will be written out to file %s every %u iterations\n",combinedgradientOFiles_[0]->getPath().c_str(),combinedgradient_wstride_);
    }
    else {
      log.printf("  Combined gradient (gradient + Hessian term) will be written out to the following files every %u iterations:\n",combinedgradient_wstride_);
      for(unsigned int i=0; i<combinedgradientOFiles_.size(); i++) {
        log.printf("   coefficient set %u: %s\n",i,combinedgradientOFiles_[i]->getPath().c_str());
      }
    }
  }
  //

  turnOnHessian();
  checkRead();
}


void Opt_BachAveragedSGD::coeffsUpdate(const unsigned int c_id) {
  //
  if(combinedgradientOFiles_.size()>0 && (getIterationCounter()+1)%combinedgradient_wstride_==0) {
    CombinedGradient(c_id).setValues( ( Gradient(c_id) + Hessian(c_id)*(AuxCoeffs(c_id)-Coeffs(c_id)) ) );
    combinedgradient_pntrs_[c_id]->setIterationCounterAndTime(getIterationCounter()+1,getTime());
    combinedgradient_pntrs_[c_id]->writeToFile(*combinedgradientOFiles_[c_id]);
  }
  //
  double aver_decay = getAverDecay();
  AuxCoeffs(c_id) += - StepSize(c_id)*CoeffsMask(c_id) * ( Gradient(c_id) + Hessian(c_id)*(AuxCoeffs(c_id)-Coeffs(c_id)) );
  //AuxCoeffs() = AuxCoeffs() - StepSize() * ( Gradient() + Hessian()*(AuxCoeffs()-Coeffs()) );
  Coeffs(c_id) += aver_decay * ( AuxCoeffs(c_id)-Coeffs(c_id) );
}


inline
double Opt_BachAveragedSGD::getAverDecay() const {
  double aver_decay = 1.0 / ( getIterationCounterDbl() + 1.0 );
  if(decaying_aver_tau_ > 0.0 && (getIterationCounterDbl() + 1.0) > decaying_aver_tau_) {
    aver_decay = 1.0 / decaying_aver_tau_;
  }
  return aver_decay;
}


}
}
