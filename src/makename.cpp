#include <iostream>
#include <sstream>
#include <cstring>

extern "C"
{
  void fctmakename_(char *resu, char *filename, float &lambda, int &nst, float &lscale);
  void newwritechar_(char *char1);
  void newwrite2char_(char *char1, char *char2);
  void newwritecharint_(char *char1, int &int1);
  void newwritechar2int_(char *char1, int &int1, int &int2);
  void newwritechar2ints_(char *char1, int &int1, int &int2);
}

std::string remove_spaces(std::string &s)
{
  int last = s.size() - 1;
  while (last >= 0 && s[last] == ' ')
    --last;
  return s.substr(0, last + 1);
}


void fctmakename_(char *resu, char *filename, float &lambda, int &nst, float &lscale)
{
  //std::cout<<filename<<std::endl;
  
  std::ostringstream oss;
  oss<<filename<<"_"<<lambda<<"_"<<nst<<"_"<<lscale;

  strcpy(resu,oss.str().data());
  //resu=const_cast< char* > ( oss.str().data() );
}

void newwritechar_(char *char1)
{
  std::string mastring=char1;
  std::cout<<" "<<remove_spaces(mastring)<<std::endl;
}

void newwrite2char_(char *char1, char *char2)
{
  std::string mastring1=char1;
  std::string mastring2=char2;
  std::cout<<" "<<remove_spaces(mastring1)<<" "<<remove_spaces(mastring2)<<std::endl;
}

void newwritecharint_(char *char1, int &int1)
{
  std::string mastring=char1;
  std::cout<<" "<<remove_spaces(mastring)<<" "<<int1<<std::endl;
}

void newwritechar2int_(char *char1, int &int1, int &int2)
{
  std::string mastring=char1;
  std::cout<<" "<<remove_spaces(mastring)<<" "<<int1<<" "<<int2<<std::endl;
}
