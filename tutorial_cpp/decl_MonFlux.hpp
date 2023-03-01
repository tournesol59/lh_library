#include <stdio.h>  // non decisive, but we use many of strncpy
#include <iostream>
#include <fstream>    // std::ifstream, std::ofstream 
#include <sstream>
#include <cstring>

/**
 *  Stream Class
 */
class MyStreamFile {
   // simply stores the name of an asciii -file
   // and function to display its content line-by-line
   // on the screen. Used to test the C++ exception dealing with
	public:
    MyStreamFile(const char *fname, int taille, char *sdata);
    ~MyStreamFile();
    int ReadFromFile() throw (int);
    int WriteData();
    // int EndStreamFile();

	protected:        
    char *szFileName;
    int taille;
    char *szData;
};

/**
 * Exception definition
 */

class MyException {
   public:
      MyException(int c); // constructor
      MyException(const MyException &);  // copy constructor, necessary to pass Exception to caller
      ~MyException(void); // destructor
      int getcode(void);
      void printex(void);
   protected:
      int code;
};

class MyFileNameExcept : public MyException {
   public:
      MyFileNameExcept(int c);
      MyFileNameExcept(const MyFileNameExcept &);
      ~MyFileNameExcept(void);
  //    int getcode(void);
      void printex(void);
};

class MyFullDirExcept : public MyException {
   public:
      MyFullDirExcept(int c);
      MyFullDirExcept(const MyFullDirExcept &);
      ~MyFullDirExcept(void);
  //    int getcode(void);
      void printex(void);
};

