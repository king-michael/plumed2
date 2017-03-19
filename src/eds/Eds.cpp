/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "bias/Bias.h"
#include "core/ActionRegister.h"
#include "tools/Random.h"
#include "tools/File.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <iostream>


using namespace PLMD;
using namespace bias;


namespace EDS{

//+PLUMEDOC BIAS EDS
/*
Add a linear bias on a set of observables.

This force is the same as the linear part of the bias in \ref RESTRAINT, 
but this bias has the ability to compute prefactors
adaptively using the scheme of White and Voth (JCTC 2014) in order to match 
target observable values for a set of CVs.

The addition to the potential is of the form 
\f[
  \sum_i \frac{\alpha_i}{s_i} x_i
\f]

where for CV \f$x_i\f$, a coupling \f${\alpha}_i\f$ is determined
adaptively or set by the user to match a target value for
\f$x_i\f$. \f$s_i\f$ is a scale parameter, which by default is set to
the target value. It may also be set separately. 

\warning
It is not possible to set the target value of the observable to zero with the default value of \f$s_i\f$ as this will cause a divide-by-zero error. Instead, set \f$s_i=1\f$ or modify the CV so the desired target value is no longer zero.

\par Examples

The following input for a harmonic oscillator of two beads will adaptively find a linear bias to change the mean and variance to the target values. The PRINT line shows how to access the value of the coupling constants. 

\verbatim 
dist: DISTANCE ATOMS=1,2
# this is the squared of the distance
dist2: COMBINE ARG=dist POWERS=2 PERIODIC=NO

#bias mean and variance
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 
PRINT ARG=dist,dist2,eds.dist_coupling,eds.dist2_coupling,eds.bias,eds.force2 FILE=colvars.dat STRIDE=100
\endverbatim

Rather than trying to find the coupling constants adaptively, can ramp up to a constant value.
\verbatim
#ramp couplings from 0,0 to -1,1 over 50000 steps
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 FIXED=-1,1 RAMP PERIOD=50000 TEMP=1.0

#same as above, except starting at -0.5,0.5 rather than default of 0,0
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 FIXED=-1,1 INIT=-0.5,0.5 RAMP PERIOD=50000 TEMP=1.0
\endverbatim

A restart file can be added to dump information needed to restart/continue simulation using these parameters every STRIDE.
\verbatim 
#add the option to write to a restart file
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 OUT_RESTART=restart.dat

#add the option to read in a previous restart file
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 IN_RESTART=restart.dat EDSRESTART

#add the option to read in a previous restart file and freeze the bias at the final level from the previous simulation
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 IN_RESTART=restart.dat EDSRESTART FREEZE

#add the option to read in a previous restart file and freeze the bias at the mean from the previous simulation
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 IN_RESTART=restart.dat EDSRESTART FREEZE MEAN

#add the option to read in a previous restart file and continue the bias, but use the mean from the previous run as the starting point
eds: EDS ARG=dist,dist2 CENTER=2.0,1.0 PERIOD=50000 TEMP=1.0 IN_RESTART=restart.dat EDSRESTART MEAN
\endverbatim


*/
//+ENDPLUMEDOC

class EDS : public Bias{

private:
  std::vector<double> center_;
  std::vector<double> scale_;
  std::vector<double> current_coupling_;
  std::vector<double> set_coupling_;
  std::vector<double> target_coupling_;
  std::vector<double> max_coupling_range_;
  std::vector<double> max_coupling_grad_;
  std::vector<double> coupling_rate_;
  std::vector<double> coupling_accum_;
  std::vector<double> means_;
  std::vector<double> ssds_;
  std::vector<Value*> out_coupling_;
  std::string in_restart_name_;
  std::string out_restart_name_;
  OFile out_restart_;
  IFile in_restart_;
  bool b_adaptive_;
  bool b_freeze_;
  bool b_equil_;
  bool b_ramp_;
  bool b_restart_;
  bool b_write_restart_;
  bool b_hard_c_range_;
  int seed_;
  int update_period_;
  int avg_coupling_count_;
  int update_calls_;
  double kbt_;
  double c_range_increase_f_;
  Random rand_;
  Value* value_force2_;

  /*read input restart. b_mean sets if we use mean or final value for freeze*/
  void readInRestart_(const bool b_mean);
  /*setup output restart*/ 
  void setupOutRestart_();
  /*write output restart*/
  void writeOutRestart_();


public:
  explicit EDS(const ActionOptions&);
  void calculate();
  void update();
  void turnOnDerivatives();
  static void registerKeywords(Keywords& keys);
  ~EDS();
};

PLUMED_REGISTER_ACTION(EDS,"EDS")

void EDS::registerKeywords(Keywords& keys){
   Bias::registerKeywords(keys);
   keys.use("ARG");
   keys.add("compulsory","CENTER","The desired centers (equilibrium values) which will be sought during the adaptive linear biasing.");
   keys.add("compulsory","PERIOD","Steps over which to adjust bias");

   keys.add("compulsory","RANGE","3.0","The largest magnitude of the force constant which one expects (in kBT) for each CV based");
   keys.add("compulsory","SEED","0","Seed for random order of changing bias");
   keys.add("compulsory","INIT","0","Starting value for coupling coefficients");
   keys.add("compulsory","FIXED","0","Fixed target values for bias factors (not adaptive)");
   keys.add("optional","BIAS_SCALE","A divisor to set the units of the bias. If not set, this will be the experimental value by default (as is done in White and Voth 2014).");
   keys.add("optional","TEMP","The system temperature. If not provided will be taken from MD code (if available)");

   keys.add("optional","OUT_RESTART","Output file for all information needed to continue EDS simulation");
   keys.add("optional","IN_RESTART","Read this file to continue an EDS simulation (if same as above, will be overwritten)");

   keys.addFlag("RAMP",false,"Slowly increase bias constant to a fixed value");
   keys.addFlag("FREEZE",false,"Fix bias at current level (only used for restarting). Can also set PERIOD=0 if not using EDSRESTART.");
   keys.addFlag("MEAN",false,"Instead of using final bias level from restart, use average");

   keys.addOutputComponent("force2","default","squared value of force from the bias");
   keys.addOutputComponent("_coupling","default","For each named CV biased, there will be a corresponding output CV_coupling storing the current linear bias prefactor.");
}

EDS::EDS(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
center_(getNumberOfArguments(),1.0),
current_coupling_(getNumberOfArguments(),0.0),
set_coupling_(getNumberOfArguments(),0.0),
target_coupling_(getNumberOfArguments(),0.0),
max_coupling_range_(getNumberOfArguments(),3.0),
max_coupling_grad_(getNumberOfArguments(),0.0),
coupling_rate_(getNumberOfArguments(),1.0),
coupling_accum_(getNumberOfArguments(),1.0),
means_(getNumberOfArguments(),0.0),
ssds_(getNumberOfArguments(),0.0),
out_coupling_(getNumberOfArguments(),NULL),
in_restart_name_(""),
out_restart_name_(""),
b_adaptive_(true),
b_freeze_(false),
b_equil_(true),
b_ramp_(false),
b_restart_(false),
b_write_restart_(false),
b_hard_c_range_(false),
seed_(0),
update_period_(0),
avg_coupling_count_(1),
update_calls_(0),
kbt_(0.0),
c_range_increase_f_(1.25),
value_force2_(NULL)
{
  double temp=-1.0;
  bool b_mean=false;
  
  addComponent("force2");
  componentIsNotPeriodic("force2");
  value_force2_=getPntrToComponent("force2");
  
  for(unsigned i=0;i<getNumberOfArguments();i++){
    std::string comp=getPntrToArgument(i)->getName()+"_coupling";
    addComponent(comp);
    componentIsNotPeriodic(comp);
    out_coupling_[i]=getPntrToComponent(comp);
  }
  
  parseVector("CENTER",center_);
  parseVector("BIAS_SCALE", scale_);
  parseVector("RANGE",max_coupling_range_);
  parseVector("FIXED",target_coupling_);
  parseVector("INIT",set_coupling_);
  parse("PERIOD",update_period_);
  parse("TEMP",temp);
  parse("SEED",seed_);
  parse("OUT_RESTART",out_restart_name_);
  parseFlag("RAMP",b_ramp_);
  parseFlag("FREEZE",b_freeze_);
  parseFlag("MEAN",b_mean);
  parse("IN_RESTART",in_restart_name_);
  parse("OUT_RESTART",out_restart_name_);
  checkRead();  

  log.printf("Setting scaling:");
  if(scale_.size() > 0  && scale_.size() < getNumberOfArguments()) {
    error("the number of BIAS_SCALE values be the same as number of CVs");
  } else if(scale_.size() == 0) {
    log.printf(" (default) ");
    
    scale_.resize(center_.size());
    for(unsigned int i = 0; i < scale_.size(); i++) {
      if(center_[i]==0) 
        error("BIAS_SCALE parameter has been set to CENTER value of 0 (as is default). This will divide by 0, so giving up. See doc for EDS bias");
      scale_[i] = center_[i];
    }
  } else {
    for(unsigned int i = 0; i < scale_.size(); i++)
      log.printf(" %f",scale_[i]);
  }
  log.printf("\n");

  if(in_restart_name_ != ""){
    b_restart_ = true;
    log.printf("  reading simulation information from file: %s\n",in_restart_name_.c_str());
    readInRestart_(b_mean);
  }
  else{

    
    if(temp>=0.0) kbt_=plumed.getAtoms().getKBoltzmann()*temp;
    else kbt_ = plumed.getAtoms().getKbT();

    //in driver, this results in kbt of 0
    if(kbt_ == 0){ 
      error("  Unable to determine valid kBT. Probably because you are runnning from driver.\n  Consider setting temperature manually.");
      kbt_ = 1;
    }
    
    log.printf("  with kBT = %f\n",kbt_);
    log.printf("  Updating every %i steps\n",update_period_);
    
    log.printf("  with centers:");
    for(unsigned i=0;i<center_.size();i++){
      log.printf(" %f",center_[i]);
    }
    
    for(unsigned int i = 0; i < scale_.size(); i++) 
      log.printf(" %f",center_[i]);
    
    log.printf("\n  with initial ranges / rates:\n");
    for(unsigned i=0;i<max_coupling_range_.size();i++) {
      //this is just an empirical guess. Bigger range, bigger grads. Less frequent updates, bigger changes
      max_coupling_range_[i]*=kbt_;
      max_coupling_grad_[i] = max_coupling_range_[i]*update_period_/100.;
      log.printf("    %f / %f\n",max_coupling_range_[i],max_coupling_grad_[i]);
    }
    
    if(seed_>0){
      log.printf("  setting random seed = %i",seed_);
      rand_.setSeed(seed_);
    }
    
    for(unsigned i=0;i<getNumberOfArguments();++i) if(target_coupling_[i]!=0.0) b_adaptive_=false;
    
    if(!b_adaptive_){
      if(b_ramp_>0) {
	log.printf("  ramping up coupling constants over %i steps\n",update_period_);
      }
      
      log.printf("  with starting coupling constants");
      for(unsigned i=0;i<set_coupling_.size();i++) log.printf(" %f",set_coupling_[i]);
      log.printf("\n");
      log.printf("  and final coupling constants");
      for(unsigned i=0;i<target_coupling_.size();i++) log.printf(" %f",target_coupling_[i]);
      log.printf("\n");
    }
    
    //now do setup
    if(b_ramp_){
      update_period_*=-1;
    }
    
    for(unsigned i=0;i<set_coupling_.size();i++) current_coupling_[i] = set_coupling_[i];
    
    // if b_adaptive_, then first half will be used for equilibrating and second half for statistics
    if(update_period_>0){
      update_period_ /= 2;
    }
    
    
  }
  
  if(b_freeze_){
    b_adaptive_=false;
    log.printf("  freezing bias at current level\n");
  }
  
  if(out_restart_name_.length()>0) {
    log.printf("  writing restart information every %i steps to file: %s\n",abs(update_period_),out_restart_name_.c_str());
    b_write_restart_ = true;
    setupOutRestart_();
  }
  
  log<<"  Bibliography "<<plumed.cite("White and Voth, J. Chem. Theory Comput. 10 (8), 3023-3030 (2014)")<<"\n";
}

void EDS::readInRestart_(const bool b_mean){
  int adaptive_i;
  
  
  in_restart_.open(in_restart_name_);
  
  //some sample code to get the field names:
  /*
    std::vector<std::string> fields;
    in_restart_.scanFieldList(fields);
    log.printf("field");
    for(unsigned i=0;i<fields.size();i++) log.printf(" %s",fields[i].c_str());
    log.printf("\n");
  */
  
  if(in_restart_.FieldExist("kbt")){
    in_restart_.scanField("kbt",kbt_);
  }else{ error("No field 'kbt' in restart file"); }
  log.printf("  with kBT = %f\n",kbt_);
  
  if(in_restart_.FieldExist("update_period")){
    in_restart_.scanField("update_period",update_period_);
  }else{ error("No field 'update_period' in restart file"); }
  log.printf("  Updating every %i steps\n",update_period_);
  
  if(in_restart_.FieldExist("adaptive")){
    //note, no version of scanField for boolean
    in_restart_.scanField("adaptive",adaptive_i);
  }else{ error("No field 'adaptive' in restart file"); }
  b_adaptive_ = bool(adaptive_i);
  
  if(in_restart_.FieldExist("seed")){
    in_restart_.scanField("seed",seed_);
  }else{ error("No field 'seed' in restart file"); }
  if(seed_>0){
    log.printf("  setting random seed = %i",seed_);
    rand_.setSeed(seed_);
  }
  
  double time;
  std::vector<double> avg_bias = std::vector<double>(center_.size());
  unsigned int N = 0;
  std::string cv_name;
  
  while(in_restart_.scanField("time",time)){
    
    for(unsigned i=0;i<getNumberOfArguments();++i) {
      cv_name = getPntrToArgument(i)->getName();
      in_restart_.scanField(cv_name + +"_center",center_[i]);
      in_restart_.scanField(cv_name + +"_scale", scale_[i]);
      in_restart_.scanField(cv_name + "_init", set_coupling_[i]);
      in_restart_.scanField(cv_name + "_target",target_coupling_[i]);
      in_restart_.scanField(cv_name + "_coupling",current_coupling_[i]);
      in_restart_.scanField(cv_name + "_maxrange",max_coupling_range_[i]);
      in_restart_.scanField(cv_name + "_maxgrad",max_coupling_grad_[i]);

      avg_bias[i] += current_coupling_[i];
      N++;
    }
    
    in_restart_.scanField();
  }

  
  log.printf("  with centers:");
  for(unsigned i=0;i<center_.size();i++) {
    log.printf(" %f",center_[i]);
  }
  log.printf("\n  and scaling:");
  for(unsigned i=0;i<scale_.size();i++) {
    log.printf(" %f",scale_[i]);
  }
  
  log.printf("\n  with initial ranges / rates:\n");
  for(unsigned i=0;i<max_coupling_range_.size();i++) {
    log.printf("    %f / %f\n",max_coupling_range_[i],max_coupling_grad_[i]);
  }
  
  if(!b_adaptive_ && update_period_<0){
    log.printf("  ramping up coupling constants over %i steps\n",-update_period_);
  }

  if(b_mean) {
    log.printf("Loaded in averages for coupling constants...\n");
    for(unsigned i=0;i<current_coupling_.size();i++) current_coupling_[i] = avg_bias[i] / N;
  }
  
  log.printf("  with current coupling constants:\n    ");
  for(unsigned i=0;i<current_coupling_.size();i++) log.printf(" %f",current_coupling_[i]);
  log.printf("\n");
  log.printf("  with initial coupling constants:\n    ");
  for(unsigned i=0;i<set_coupling_.size();i++) log.printf(" %f",set_coupling_[i]);
  log.printf("\n");
  log.printf("  and final coupling constants:\n    ");
  for(unsigned i=0;i<target_coupling_.size();i++) log.printf(" %f",target_coupling_[i]);
  log.printf("\n");
  
  in_restart_.close();
}
  
void EDS::setupOutRestart_(){
  out_restart_.link(*this);
  out_restart_.open(out_restart_name_);
  out_restart_.setHeavyFlush();
  
  out_restart_.addConstantField("adaptive").printField("adaptive",b_adaptive_);
  out_restart_.addConstantField("update_period").printField("update_period",update_period_);
  out_restart_.addConstantField("seed").printField("seed",seed_);
  out_restart_.addConstantField("kbt").printField("kbt",kbt_);
  
  writeOutRestart_();
}
  
void EDS::writeOutRestart_(){
  std::string cv_name;
  out_restart_.printField("time",getTimeStep()*getStep());
  
  for(unsigned i=0;i<getNumberOfArguments();++i) {
    cv_name = getPntrToArgument(i)->getName();
    
    out_restart_.printField(cv_name + "_center",center_[i]);
    out_restart_.printField(cv_name + "_scale",scale_[i]);
    out_restart_.printField(cv_name + "_init",set_coupling_[i]);
    out_restart_.printField(cv_name + "_target",target_coupling_[i]);
    out_restart_.printField(cv_name + "_coupling",current_coupling_[i]);
    out_restart_.printField(cv_name + "_maxrange",max_coupling_range_[i]);
    out_restart_.printField(cv_name + "_maxgrad",max_coupling_grad_[i]);
  }
  out_restart_.printField();
}
  
void EDS::calculate(){
  unsigned int ncvs = getNumberOfArguments();
  
  //Compute linear force as in "restraint"
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<ncvs;++i){
    const double cv=difference(i,center_[i],getArgument(i));
    const double m=current_coupling_[i];
    const double f=-m;
    ene+=m*cv;
    setOutputForce(i,f);
    totf2+=f*f;
  };
  setBias(ene);
  value_force2_->set(totf2);
  
  //adjust parameters according to EDS recipe
  update_calls_++;

  //check if we're ramping or doing normal updates and then restart if needed. The ramping check
  //is complicated because we could be frozen, finished ramping or not ramping.
  //The + 2 is so we have an extra line showing that the bias isn't changing (for my sanity and yours)
  if( b_write_restart_ &&
      ( (update_period_ < 0 && !b_freeze_ && update_calls_ <= fabs(update_period_) + 2) ||
	(update_period_ > 0 && update_calls_ % update_period_ == 0 ) ) ) {
    writeOutRestart_();
  }
  
  int b_finished_equil_flag = 1;
  double delta;
  
  //assume forces already applied and saved
  
  for(unsigned i=0;i<ncvs;++i){
    //are we ramping to a constant value and not done equilibrating
    if(update_period_ < 0){
      if(update_calls_ <= fabs(update_period_) && !b_freeze_){
	current_coupling_[i] += (target_coupling_[i]-set_coupling_[i])/fabs(update_period_);
      }
      //make sure we don't reset update calls
      b_finished_equil_flag = 0;
      continue;
    } 
    //not updating
    else if(update_period_==0){
      continue;
    }
    
    //if we aren't wating for the bias to equilibrate, collect data
    if(!b_equil_){
      //Welford, West, and Hanso online variance method
      delta = difference(i,means_[i],getArgument(i));
      means_[i]+=delta/update_calls_;
      ssds_[i]+=delta*difference(i,means_[i],getArgument(i));
    }
    // equilibrating
    else {
      //check if we've reached the setpoint
      if(coupling_rate_[i]==0 || pow(current_coupling_[i]-set_coupling_[i],2)<pow(coupling_rate_[i],2)) {
	b_finished_equil_flag &= 1;
      }
      else{
	current_coupling_[i]+=coupling_rate_[i];
	b_finished_equil_flag = 0;
      }
    }
    //Update max coupling range if allowed
    if(!b_hard_c_range_ && fabs(current_coupling_[i])>max_coupling_range_[i]) {
      max_coupling_range_[i]*=c_range_increase_f_; 
      max_coupling_grad_[i]*=c_range_increase_f_; 
    }
  }
  
  //reduce all the flags 
  if(b_equil_ && b_finished_equil_flag) {
    b_equil_ = false;
    update_calls_ = 0;
  }
  
  //Now we update coupling constant, if necessary
  if(!b_equil_ && update_period_ > 0 && update_calls_ == update_period_ && !b_freeze_) {
    double step_size = 0;
    double tmp;
    for(unsigned i=0;i<ncvs;++i){
      //calulcate step size
      //uses scale here, which by default is center
      tmp = 2. * (means_[i]/scale_[i] - 1) * ssds_[i] / (update_calls_ - 1);
      step_size = tmp / kbt_;
      
      //check if the step_size exceeds maximum possible gradient
      step_size = copysign(fmin(fabs(step_size), max_coupling_grad_[i]), step_size);
      
      //reset means/vars
      means_[i] = 0;
      ssds_[i] = 0;
      
      //multidimesional stochastic step
      if(ncvs == 1 || (rand_.RandU01() < (1. / ncvs) ) ) {
	coupling_accum_[i] += step_size * step_size;

	//equation 5 in White and Voth, JCTC 2014
	//no negative sign because it's in step_size
	set_coupling_[i] += max_coupling_range_[i]/sqrt(coupling_accum_[i])*step_size;
	coupling_rate_[i] = (set_coupling_[i]-current_coupling_[i])/update_period_;

      } else {
	//do not change the bias
	coupling_rate_[i] = 0;
      }
      
      
    }
    
    update_calls_ = 0;
    avg_coupling_count_++;
    b_equil_ = true; //back to equilibration now
  } //close update if
  
    //pass couplings out so they are accessible
  for(unsigned i=0;i<ncvs;++i){
    out_coupling_[i]->set(current_coupling_[i]);
  }
}
  
void EDS::update(){
  //pass
}
  
EDS::~EDS(){
  out_restart_.close();
}
  
void EDS::turnOnDerivatives(){
  // do nothing
  // this is to avoid errors triggered when a bias is used as a CV
  // (This is done in ExtendedLagrangian.cpp)
}
  
  
}
