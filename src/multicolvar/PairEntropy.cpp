/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR PAIRENTROPY
/*
Calculate the coordination numbers of atoms so that you can then calculate functions of the distribution of
coordination numbers such as the minimum, the number less than a certain quantity and so on.   

To make the calculation of coordination numbers differentiable the following function is used:

\f[
s = \frac{ 1 - \left(\frac{r-d_0}{r_0}\right)^n } { 1 - \left(\frac{r-d_0}{r_0}\right)^m }
\f]

\par Examples

The following input tells plumed to calculate the coordination numbers of atoms 1-100 with themselves.
The minimum coordination number is then calculated.
\verbatim
PAIRENTROPY SPECIES=1-100 R_0=1.0 MIN={BETA=0.1}
\endverbatim

The following input tells plumed to calculate how many atoms from 1-100 are within 3.0 of each of the atoms
from 101-110.  In the first 101 is the central atom, in the second 102 is the central atom and so on.  The 
number of coordination numbers more than 6 is then computed.
\verbatim
PAIRENTROPY SPECIESA=101-110 SPECIESB=1-100 R_0=3.0 MORE_THAN={RATIONAL R_0=6.0 NN=6 MM=12 D_0=0}
\endverbatim

*/
//+ENDPLUMEDOC


class PairEntropy : public MultiColvarBase {
private:
//  double nl_cut;
  double rcut2;
  double invSqrt2piSigma, sigmaSqr2;
  double maxr, nhist,sigma;
  double deltar;
  unsigned deltaBin;
public:
  static void registerKeywords( Keywords& keys );
  explicit PairEntropy(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ; 
  virtual double kernel(double distance, double&dfunc)const;
  virtual double integrate(vector<double> integrand, double delta)const;
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(PairEntropy,"PAIRENTROPY")

void PairEntropy::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","MAXR","1","Maximum distance for the radial distribution function ");
  keys.add("compulsory","NHIST","300","Number of bins in the rdf ");
  keys.add("compulsory","SIGMA","0.1","Width of gaussians ");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); 
}

PairEntropy::PairEntropy(const ActionOptions&ao):
Action(ao),
MultiColvarBase(ao)
{
  parse("MAXR",maxr);
  parse("NHIST",nhist);
  parse("SIGMA",sigma);

  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); checkRead();

  // Define heavily used constants
  double sqrt2piSigma = std::sqrt(2*pi)*sigma;
  invSqrt2piSigma = 1./sqrt2piSigma;
  sigmaSqr2 = 2.*sigma*sigma;
  deltar=maxr/nhist;
  deltaBin = std::floor(3*sigma/deltar); //3*sigma is 99.7 %

  // Set the link cell cutoff
  setLinkCellCutoff( maxr + 3*sigma );
  rcut2 = (maxr + 3*sigma)*(maxr + 3*sigma);
  log.printf("Setting cut off to %f \n ", maxr + 3*sigma );
}

double PairEntropy::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
   // Calculate the coordination number
   double dfunc, d2;
   // Construct g(r)
   vector<double> gofr(nhist);
   for(unsigned i=1;i<myatoms.getNumberOfAtoms();++i){
      Vector& distance=myatoms.getPosition(i);  
      if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
           double distanceModulo=std::sqrt(d2);
           unsigned bin=std::floor(distanceModulo/deltar);
           int minBin, maxBin;
           // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
           minBin=bin - deltaBin;
           if (minBin < 0) minBin=0;
           if (minBin > (nhist-1)) minBin=nhist-1;
           maxBin=bin +  deltaBin;
           if (maxBin > (nhist-1)) maxBin=nhist-1;
           for(unsigned int j=minBin;j<maxBin+1;j+=1) {   
             double x=deltar*(j+0.5);
             gofr[j] += kernel(x-distanceModulo, dfunc);
             //log.printf("Distance x gofr maxBin minBin %f %f %f %u %u \n ", distanceModulo, x, kernel(x-distanceModulo, dfunc), maxBin, minBin);
	   } 
      }
   }
   // Normalize g(r)
   double volume=getBox().determinant(); 
   double density=getNumberOfAtoms()/volume;
   for(unsigned i=0;i<nhist;++i){
     double x=deltar*(i+0.5);
     gofr[i] /= 4*pi*density*x*x;
     //log.printf("x gofr %f %f \n ", x, gofr[i]);
   }
   // Construct integrand
   vector<double> integrand(gofr.size());
   for(unsigned i=0;i<gofr.size();++i){
     double x=deltar*(i+0.5);
     if (gofr[i]<1.e-10) {
       integrand[i] = x*x;
     } else {
       integrand[i] = (gofr[i]*std::log(gofr[i])-gofr[i]+1)*x*x;
     }
   }
   // Integrate to obtain pair entropy;
   double entropy = -2*pi*density*integrate(integrand,deltar); 
   return entropy;
}

double PairEntropy::kernel(double distance,double&dfunc)const{
  // Gaussian function and derivative
  double result = invSqrt2piSigma*std::exp(-distance*distance/sigmaSqr2) ;
  dfunc = 0.;
  return result;
}

double PairEntropy::integrate(vector<double> integrand, double delta)const{
  // Trapezoid rule
  double result = 0.;
  for(unsigned i=1;i<(integrand.size()-1);++i){
    result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()];
  result *= delta;
  return result;
}


}
}

