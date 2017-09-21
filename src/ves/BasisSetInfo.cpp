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

#include "BasisSetInfo.h"
#include "BasisFunctions.h"


namespace PLMD {
namespace ves {

BasisSetInfo::BasisSetInfo(
  const std::string& basisset_label,
  std::vector<BasisFunctions*>& basisf):
  basisset_label_(basisset_label),
  basisf_(basisf),
  basisset_dim_(basisf.size()),
  basisset_size_(1),
  bf_keywords_(basisset_dim_),
  bf_types_(basisset_dim_),
  bf_orders_(basisset_dim_),
  bf_sizes_(basisset_dim_),
  bf_intervals_min_(basisset_dim_),
  bf_intervals_max_(basisset_dim_),
  bf_intervals_range_(basisset_dim_),
  bf_intervals_mean_(basisset_dim_),
  bf_intervals_bounded_(basisset_dim_),
  bf_intervals_periodic_(basisset_dim_),
  basisset_volume_(1.0)
{
  // initialize();
}

//
// void BasisSetInfo::initialize() {
//   plumed_assert(basisf_.size()>0);
//   //
//   basisset_size_=1;
//   basisset_volume_=1.0;
//   for(unsigned int k=0;k<basisset_dim_;k++){
//     bf_types_[k]=basisf_[k]->getType();
//     bf_orders_[k]=basisf_[k]->getOrder();
//     bf_sizes_[k]=basisf_[k]->getNumberOfBasisFunctions();
//     basisset_size_*=bf_sizes_[k];
//     bf_intervals_min_[k]=basisf_[k]->intervalMin();
//     bf_intervals_max_[k]=basisf_[k]->intervalMax();
//     bf_intervals_range_[k]=basisf_[k]->intervalRange();
//     bf_intervals_mean_[k]=basisf_[k]->intervalMean();
//     bf_intervals_bounded_[k]=basisf_[k]->intervalBounded();
//     basisset_volume_*=bf_intervals_range_[k];
//     bf_intervals_periodic_[k]=basisf_[k]->arePeriodic();
//     bf_keywords_[k]=basisf_[k]->getKeywordString();
//   }
// }
//
//
// unsigned int BasisSetInfo::getDimension() const {
//   return basisset_dim_;
// }
//
//
// size_t BasisSetInfo::getSize() const {
//   return basisset_size_;
// }
//
//
// double BasisSetInfo::getVolume() const {
//   return basisset_volume_;
// }
//
//
// std::vector<std::string> BasisSetInfo::getTypes() const {
//   return bf_types_;
// }
//
//
// std::vector<unsigned int> BasisSetInfo::getOrders() const {
//   return bf_orders_;
// }
//
//
// std::vector<unsigned int> BasisSetInfo::getSizes() const {
//   return bf_sizes_;
// }
//
//
// std::vector<double> BasisSetInfo::getMinima() const {
//   return bf_intervals_min_;
// }
//
//
// std::vector<double> BasisSetInfo::getMaxima() const {
//   return bf_intervals_max_;
// }
//
//
// std::vector<std::string> BasisSetInfo::getKeywords() const {
//   return bf_keywords_;
// }
//
//
// std::string BasisSetInfo::getBasisSetDescription(std::vector<unsigned int>& indices) const {
//   plumed_massert(indices.size()==basisset_dim_,"incorrect number of indicies");
//   std::string desc;
//   desc=basisf_[0]->getBasisFunctionLabel(indices[0]);
//   for(unsigned int k=1;k<basisset_dim_;k++){ desc+="*"+basisf_[k]->getBasisFunctionLabel(indices[k]); }
//   return desc;
// }


}
}
