/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

using namespace std;

namespace PLMD {
namespace vatom {

//+PLUMEDOC VATOM ROUSE_MODE_COORDINATE
/*
Calculate the Rouse mode coordinate for a group of atoms. 

This is done by using the formula 
\f[
\mathbf{X}_p = \sqrt{ \frac{c_{p}}{2} } 
\sum_{i=0}^{N-1} 
\mathbf{R}_i 
\cos 
\left[ 
\frac{p \pi}{N}
\left(
i+\frac{1}{2}
\right)
\right]
\f]
where \f$N\f$ is the number of atoms considered, \f$p=0,\ldots,N-1\f$ is the 
Rouse mode index, \f$\mathbf{R}_i\f$ is the position of the \f$i\f$-th atom, 
and \f$c_{p}\f$ is a constant equal 1 for \f$p=0\f$ and 2 otherwise. 

It is assumed that the atoms given in ATOMS are in the right order. 

The computed Rouse mode coordinate is stored as a virtual atom that can be accessed in
an atom list through the label for the ROUSE_MODE_COORDINATE action that creates it.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically 
by considering the ordered list of atoms and rebuilding PBCs with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES. Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

The following input instructs plumed to calculate the Rouse mode coordinate for mode 1:
\plumedfile
r1: ROUSE_MODE_COORDINATE ATOMS=1-128 MODE=1
\endplumedfile

The following input instructs plumed to calculate a collective variable that is the norm of the Rouse mode coordinate:
\plumedfile
r1: ROUSE_MODE_COORDINATE ATOMS=1-128 MODE=1
p_r1: POSITION ATOM=r1
n_r1: CUSTOM ARG=p_r1.x,p_r1.y,p_r1.z FUNC=sqrt(x^2+y^2+z^2) PERIODIC=NO
\endplumedfile





*/
//+ENDPLUMEDOC


class RouseModeCoordinate:
  public ActionWithVirtualAtom
{
  bool nopbc;
  double prefactor_;
  double cosine_factor_;
public:
  explicit RouseModeCoordinate(const ActionOptions&ao);
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(RouseModeCoordinate,"ROUSE_MODE_COORDINATE")

void RouseModeCoordinate::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions");
  keys.add("compulsory","MODE","1","The Rouse mode to be calculated");
}

RouseModeCoordinate::RouseModeCoordinate(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  nopbc(false),
  prefactor_(1.0),
  cosine_factor_(0.0)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("at least one atom should be specified");
  parseFlag("NOPBC",nopbc);
  unsigned int rouse_mode = 1;
  parse("MODE",rouse_mode);
  checkRead();
  
  unsigned int N = atoms.size();
  if(rouse_mode>N-1){plumed_merror("the Rouse mode index cannot be larger than N-1 where N is the number of atoms considered");}
  cosine_factor_ = (rouse_mode*pi) / N;  
  if(rouse_mode==0){
    prefactor_ = sqrt(0.5);
  }
  
  log.printf("  Rouse mode %u calculated for atoms",rouse_mode);
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(i%25==0) log<<"\n";
    log.printf(" %d",atoms[i].serial());
  }
  log.printf("\n");
  if(nopbc) {
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }
  requestAtoms(atoms);
}

void RouseModeCoordinate::calculate() {
  Vector coord;
  if(!nopbc) makeWhole();
  
  vector<Tensor> deriv(getNumberOfAtoms());
  
  for(unsigned i=0; i<getNumberOfAtoms(); i++) {
    double cos_tmp = cos(cosine_factor_*(i+0.5));
    coord += prefactor_*cos_tmp*getPosition(i);
    deriv[i] = prefactor_*cos_tmp*Tensor::identity();
  }
  setPosition(coord);
  setAtomsDerivatives(deriv);
  
  setMass(1.0);  
  setCharge(0.0);
}

}
}
