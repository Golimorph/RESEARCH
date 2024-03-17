#include "rand_gen.h"
/***********************************************************************/
rand_class::rand_class(){
/***********************************************************************/
  x=521288629;
  y=362436069;
  z=16163801; 
  c=1; 
  n=1131199209;
  max_unint=4294967295;
  scale=1.0/static_cast<double>(max_unint);
  x_init = x;
  y_init = y;
  z_init = z;
  n_init = n;
}
/***********************************************************************/
void rand_class::set_seed(const unint i,const unint j,const unint k,
			  const unint l){
/***********************************************************************/
  x=i;y=j;z=k;n=l;
    
  c=y>z;
  for (int p=1;p<=4;p++){
      for (int q=1;q<=10;q++){r_unint();}
  }

}  
/***********************************************************************/
void rand_class::get_seed(unint i,unint j,unint k,unint l){
/***********************************************************************/
  i=x;j=y;k=z;l=n;
}  
/***********************************************************************/
void rand_class::find_limits(unint nr,unint *bin, unint *r_max){
/***********************************************************************/
  *r_max=max_unint-(max_unint-nr+1)%nr;
  *bin=(*r_max-nr+1)/nr+1; 
}
/***********************************************************************/
//int rand_class::find_shift(unint i){
/***********************************************************************/
//return (static_cast<int>(32-log(static_cast<double>(i))/log(2)+0.01));
//}  
/***********************************************************************/
int rand_class::find_shift(unint i){
/***********************************************************************/
  if (i==2) return(31);if (i==4) return(30);if (i==8) return(29);
  if (i==16) return(28);if (i==32) return(27);if (i==64) return(26);
  if (i==128) return(25);if (i==256) return(24);if (i==512) return(23);
  if (i==1024) return(22);if (i==2048) return(21);if (i==4096) return(20);
  if (i==8192) return(19);if (i==16384) return(18);if (i==32768) return(17);
  if (i==65536) return(16);if (i==131072) return(15);if (i==262144) return(14);
  if (i==524288) return(13);if (i==1048576) return(12);if (i==2097152) return(11);
  if (i==4194304) return(10);if (i==8388608) return(9);
  if (i==16777216) return(8);if (i==33554432) return(7);
  if (i==67108864) return(6);if (i==134217728) return(5);
  if (i==268435456) return(4);if (i==536870912) return(3);
  if (i==1073741824) return(2);if (i==2147483648) return(1);
  //if (i==4294967296) return(0);
  return(-1);
}  
/***********************************************************************/
void rand_class::read_seed(){
/***********************************************************************/
  std::fstream data_file;
  data_file.open("Some user path",std::ios::in);
  data_file >> x>>y>>z>>n;
  c=y>z;
  data_file.close();
    data_file.open("Some user path",std::ios::out);
  for (int i=1;i<=4;i++){
    unint r=0;
    for (int j=1;j<=10;j++){r=r_unint();}
    data_file.width(10);
    data_file << r <<"\n";
  }
  data_file.close();
}
/***********************************************************************/
void rand_class::read_seed_local(std::string path){
/***********************************************************************/
  std::fstream data_file;
  data_file.open(path.c_str(),std::ios::in);
    
/* for MPI use
  if (!data_file.is_open()){
    printf("RNG seedfile not found, exiting\n");
    MPI_Abort(MPI_COMM_WORLD,17);
  }
 
*/
    
  data_file >> x>>y>>z>>n;
  x_init = x; 
  y_init = y; 
  z_init = z; 
  n_init = n; 
  
  c=y>z;
  
  data_file.close();
  data_file.open(path.c_str(),std::ios::out);
  for (int i=1;i<=4;i++){
    unint r=0;
    for (int j=1;j<=10;j++){r=r_unint();}
    data_file.width(10);
    data_file << r <<"\n";
  }

  data_file.close();
    
}
/***********************************************************************/
void rand_class::write_seed_local(){
/***********************************************************************/
  std::fstream data_file;
  data_file.open("InputData/inputSeed.txt",std::ios::out);
  for (int i=1;i<=4;i++){
    unint r=0;
    for (int j=1;j<=i*211;j++){r=r_unint();}
    data_file.width(10);
    data_file << r <<"\n";
  }
  data_file.close();
}
/***********************************************************************/
std::vector<unint> rand_class::get_current_seed_vector(){
/***********************************************************************/
    //Can only be called after calling read_seed_local()
    std::vector<unint> seed_vector;
    seed_vector.resize(0);
    seed_vector.push_back(x_init);
    seed_vector.push_back(y_init);
    seed_vector.push_back(z_init);
    seed_vector.push_back(n_init);
    
    return seed_vector;
}
/***********************************************************************/
std::vector<unint> rand_class::get_seed_vector(){
/***********************************************************************/
    std::vector<unint> seed_vector;
    for (int i=1;i<=4;i++){
        unint r = 0;
        for (int j=1;j<=i*211;j++){r=r_unint();}
        seed_vector.push_back(r);
    }
    return seed_vector;
}
/***********************************************************************/
void rand_class::print_seed(int myrank,int pos){
/***********************************************************************/
    printf("%u %u %u %u %i %i\n",x,y,z,n,myrank,pos);
}
/***********************************************************************/
std::ostream  &operator<<(std::ostream& file,const rand_class& rnd){
/***********************************************************************/
  file<<rnd.x<<" "<<rnd.y<<" "<<rnd.z<<" "<<rnd.c<<" "<<rnd.n<<"\n";
  return(file);
}
/***********************************************************************/
std::istream  &operator>>(std::istream& file,rand_class& rnd){
/***********************************************************************/
  file>>rnd.x>>rnd.y>>rnd.z>>rnd.c>>rnd.n;
  return(file);
}
