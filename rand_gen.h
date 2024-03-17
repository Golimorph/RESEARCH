/***********************************************************************
Patrik Henelius March 20 2000
Mikael Twengström, Extended and modified, October 2014 
Own randin # class
************************************************************************
The following functions are provided for optimized speed/convenience:

rand_class rnd;            class constructor
rnd.set_seed(i,j,k,l) 
                      sets the seed (type unsigned int) to i,j,k,l.
rnd.get_seed(i,j,k,l) 
                      returns seed as  i,j,k,l.
rnd.read_seed()       
                      reads seed from file, overrites with new seeds.
rnd.read_seed_local()       
                      reads seed from local file, overrites with new seeds.
rnd.write_seed_local()       
                      writes with new seeds to local file, generates 
                      different seeds than read_seed_local.
unint rnd.r_unint()   
                      returns unsigned int between 0 and 4294967295.
unint rnd.r_unint0()   
                      returns unsigned int between 1 and 4294967295.
unint rnd.r_unint1()   
                      returns unsigned int between 0 and 4294967294.
unint rnd.r_unint01()   
                      returns unsigned int between 1 and 4294967294.
unint r_pow2(unint i) 
                      returns (unsigned) integer between 0 and i-1, works
                      only if i is power of 2.
int find_shift(unint i)
                      returns the number of right shifts needed to 
                      return an unint between 0 and i-1. Assumes i is a 
		      power of 2, otherwise returns -1.
unint r_shift(int nshift) 
                      returns unsigned int between 0 and 4294967295,
                      shifted nshift steps to the right.

unint r_int(unint i)  
                      returns (unsigned) integer between 0 and i-1.
                      works for i=1..4294967295. 
find_limits(unint i, unint bin, unint r_max) 
                      returns the bin size and maximum approved random 
                      number for generator to make evenly distributed
                      unints from 0 to i-1. 
unint r_int(unint bin, unint r_max) 
                      returns (unsigned) integer determined by bin and r_max.
                      works for i=2..4294967295.

unint r_2() ..r_1024() 
                      returns unint between 0 and i-1, with i=2..1024.
double r_double01()     
                      returns double value (0..1) (excludes both endpoints)
double r_double0()     
                      returns double value (0..1] (excludes 0)
double r_double1()     
                      returns double value [0..1) (excludes 1)
double r_double()     
                      returns double value [0..1] (includes endpoints)

*************************************************************************
It is faster to make int randomnumbers when ii is a power of 2:
The fastest way to generate int randomnumbers between 0 and ii-1, when
ii is a power of  2 is:
int nshifts=rnd.find_shifts(ii);
for(i=1;i<=10000;++i){ sum+=rnd.r_shift(nshifts)} 

A more convenient form that is slower (used when making single calls):
sum=rnd.r_pow2(ii)} 
**************************************************************************
The fastest way to generate many int randomnumbers between 0 and ii-1 is:
unsigned int bin, r_max
rnd.find_limits(ii,&bin,&r_max);
for(i=1;i<=10000;++i){ sum+=rnd.r_int(bin,r_max)} 

An more convenient way that is slower (used with single calls):
sum=rnd.r_int(ii)
**************************************************************************
The fastest way to return 0 or 1 is: rnd.r_2().

The fastest routine returning doubles is rnd.double(), 
which includes both endpoints.

If an event is accepted with probability p, use:
if(rnd.r_double0<=p) accept;
this way the event cannot be accepted if p=0.

************************************************************************/
#ifndef RAND_OWN_H
#define RAND_OWN_H
#define __USE_STD_IOSTREAM
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
//#include "mpi.h"
//using namespace std;
typedef unsigned int unint;
/***********************************************************************/
class rand_class {
/***********************************************************************/
  friend std::ostream &operator<<(std::ostream&,const rand_class&);
  friend std::istream &operator>>(std::istream&,rand_class&);
 private:
  unint x,y,z,c,n,max_unint;
  unint x_init,y_init,z_init,n_init;
  double scale;
 public:
  rand_class();
  void read_seed(); 
  void read_seed_local(std::string);
  void write_seed_local(); 
  void set_seed(const unint i,const unint j,const unint k,const unint l);
  void get_seed(unint i,unint j,unint k,unint l);
  int find_shift(unint i);
  void find_limits(unint i,unint *bin, unint *r_max);
  std::vector<unint> get_seed_vector();
  std::vector<unint> get_current_seed_vector();
  void print_seed(int,int);

  void rand_step(){
    int s;
    if (y>x+c) {s=y-(x+c);c=0;}
    else {s=y-(x+c)-18; c=1;}
    x=y;y=z;z=s;n=69069*n+1013904243;
  }

  unint r_unint(){
    rand_step();
    return (z+n);
  }

  unint r_unint0(){
    while(1){
      rand_step();
      if(z+n) return (z+n);
    }
  }

  unint r_unint1(){
    while(1){
      rand_step();
      if ((z+n)<max_unint) return (z+n);
    }
  }

  unint r_unint01(){
    while(1){
      rand_step();
      if ((z+n)!=0&&(z+n)!=max_unint) return (z+n);
    }
  }

  unint r_shift(int nshift){
    rand_step();
    return ((z+n)>>nshift);
  }

  unint r_pow2(unint i){
    int ii=find_shift(i);
    rand_step();
    return ((z+n)>>ii);
  }

  unint r_int(unint bin, unint r_max){
    while (1) {
      rand_step();
      if ((z+n)<r_max) return ((z+n)/bin);
    }
  }

  unint r_int(unint nr){
    if (nr==1) return(0);
    unint r_max=max_unint-(max_unint-nr+1)%nr;
    unint bin=(r_max-nr+1)/nr+1;
    while (1) {
      rand_step();
      if ((z+n)<r_max) return ((z+n)/bin);
    }
  }
  
  double r_double(){
    rand_step();
    return ((static_cast<double>(z+n))*scale);
  }

  double r_double1(){
    while (1) {
      rand_step();
      if ((z+n)<max_unint){
	return ((static_cast<double>(z+n))*scale);
      }
    }
  }

  double r_double0(){
    while (1) {
      rand_step();
      if (z+n){
	return ((static_cast<double>(z+n))*scale);
      }
    }
  }

  double r_double01(){
    while (1) {
      rand_step();
      if ((z+n)!=0&&(z+n)!=max_unint){
	return ((static_cast<double>(z+n))*scale);
      }
    }
  }

  unint r_2(){
    rand_step();
    return ((z+n)>>31);
  }
  unint r_4(){
    rand_step();
    return ((z+n)>>30);
  }
  unint r_8(){
    rand_step();
    return ((z+n)>>29);
  }
  unint r_16(){
    rand_step();
    return ((z+n)>>28);
  }
  unint r_32(){
    rand_step();
    return ((z+n)>>27);
  }
  unint r_64(){
    rand_step();
    return ((z+n)>>26);
  }
  unint r_128(){
    rand_step();
    return ((z+n)>>25);
  }
  unint r_256(){
    rand_step();
    return ((z+n)>>24);
  }
  unint r_512(){
    rand_step();
    return ((z+n)>>23);
  }
  unint r_1024(){
    rand_step();
    return ((z+n)>>22);
  }

};
#endif
/***********************************************************************/
