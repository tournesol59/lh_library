#include "decl_MonFlux.hpp"

#define LH_EX_READ    1
#define LH_EX_WRITE   2

MyStreamFile::MyStreamFile(const char* fname, int dtaille, char *sdata) :
	taille(dtaille)
{
   strncpy(szFileName, fname, 14);  // fname must be of size <= 14 
   szData = new char[256];
   strncpy(szData, (const char*) sdata, dtaille);
   
}

MyStreamFile::~MyStreamFile()
{
}

int MyStreamFile::ReadFromFile() throw (int)
{
        /**
	 * print file as a whole version
	 */
   try {
      std::ifstream is((const char*) szFileName, std::ifstream::binary);
      // get length of file
      is.seekg(0, is.end);
      int length = is.tellg();
      is.seekg(0, is.beg);

      char * buffer = new char[length];

      std::cout << "Reading " << length<< " characters..." << std::endl;
      // read data as a block
      is.read(buffer, length);
      buffer[length]='\0'; // add a null characte at the end of the buffer
      if (is) {
	 std::cout << "all characters read sucessfully." << std::endl;
         std::cout << buffer << std::endl;
      }
      else
	 std::cout << "error: only ' << is.gcount() << ' could be read" << std::endl;
      is.close();
      // .. buffer contains the entire file...
      delete[] buffer;

   } catch (int e) {
      std::cout << "An error has occurred at the opening of the file for reading" << std::endl;
      throw LH_EX_READ;  // 1 this will be our code of bad read  
   }
	/**
	 * line by line version
   std::fstream :: fs.open((const char*) szFileName, std::fstream::in );
   string szline;
   int flag=std::getline(fs, szline); 
   while (flag != eof()) {
      std::stringstream sts(szline);  // decompose, separated by \t
      for (int i=0; i<3; i++) {
          std::cout <<  sts << "\t";
      }
      std::cout << std:endl;
      flag=std::getline(fs, szline); 
    }
   fs.close();  // close file
        */

   return 0;
}

int MyStreamFile::WriteData()
{
   //write sdata array of char into file at the end

   try {
      std::ifstream infile ((const char*) szFileName, std::ifstream::binary);
      std::ofstream outfile ((const char*) szFileName, std::ofstream::binary);
//      fs.open((const char*) szFileName, std::fstream::in | std::fstream::out);
      // get size of the file
      infile.seekg(0, infile.end);
      long size = infile.tellg();
      infile.seekg(0);

      // allocate memory for file content + preview for data
      char* buffer = new char[size+taille]; 
      // add data, for this target the buffer at position size
      char* temp = &buffer[size];
      strncpy(temp, (const char*) szData, taille);
      // write to outfile
      outfile.write(buffer,size+taille);
      // release dynamically-allocated memory
      delete[] buffer;
      
      outfile.close();
      infile.close();  // close file
   }
   catch (int e) {
      std::cout << "A problem in writing at the end of the file has occurred " << std::endl;
//      throw LH_EX_WRITE;
   }
   return 0;
}

/******* Class MyException *********/

MyException :: MyException(int c) : // constructor
   code(c)
{

}

MyException :: MyException(const MyException &source)  // copy constructor, necessary to pass Exception to caller
{
   code=source.code;
}

MyException :: ~MyException(void) 
{ // destructor

}

int MyException :: getcode(void) 
{
   return code;
}

void MyException :: printex(void) 
{
   std::cout << "An exception was thrown with code " << code << std::endl;
}

/* special exception "Invalid File name for opening" */

MyFileNameExcept :: MyFileNameExcept(int c) : MyException(c)
{

}

MyFileNameExcept :: MyFileNameExcept(const MyFileNameExcept &source) : MyException(source)
{

}

MyFileNameExcept :: ~MyFileNameExcept(void)
{

}

void MyFileNameExcept :: printex(void) 
{
   std::cout << "An invalid file name exception was thrown with code " << code << std::endl;
}

/*  other exception reserved to error in file writing */
MyFullDirExcept :: MyFullDirExcept(int c) :  MyException(c)
{

}

MyFullDirExcept :: MyFullDirExcept(const MyFullDirExcept &source) : MyException(source)
{
	
}

MyFullDirExcept :: ~MyFullDirExcept(void) 
{

}

void MyFullDirExcept :: printex(void)
{
   std::cout << "An invalid writing file exception was thrown with code " << code << std::endl;
}

