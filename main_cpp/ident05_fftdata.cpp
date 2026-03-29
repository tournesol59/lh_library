#include "../include/ident05_fftdata.hpp" 
#include <iterator>
#include <numeric>
#include <algorithm>
#include <cstring>
#include <string>

// constructor
 IDENT05_IODATA::IDENT05_IODATA(Index size, Number fTs, const char* giveFileName) : 
   sizeevalsol(size),
   Ts(fTs)
 { 
  sizefft=2^(8);
  samplefft=2*3.1415/Ts;
  //strncpy(FileName, "import.dat", 10);
  strncpy(FileName, giveFileName, 15);
  }
 
//cpy constructor
 IDENT05_IODATA::IDENT05_IODATA(const IDENT05_IODATA & idsrc) :
   sizeevalsol(idsrc.sizeevalsol),
   Ts(idsrc.Ts)
{
  sizefft=2^(8);
  samplefft=2*3.1415/Ts;
  //strncpy(FileName, "import.dat", 9);
  strncpy(FileName, (const char*) idsrc.FileName, 15);
}

// affect constructor
 IDENT05_IODATA &IDENT05_IODATA::operator=(const IDENT05_IODATA & idsrc) // :
{
   sizeevalsol = idsrc.sizeevalsol;
   Ts = idsrc.Ts;
	
  sizefft=2^(8);
  samplefft=2*3.1415/Ts;
  //strncpy(FileName, "import.dat", 10);
  strncpy(FileName, (const char*) idsrc.FileName, 15);
}

// destructor
 IDENT05_IODATA::~IDENT05_IODATA()
 {
	 
 }
 
 bool IDENT05_IODATA::read_extern_output(Index dim, li_doubles &li_reals)
 {
   std::pair<double,double> singleton; // TBRD
   //std::cout << FileName << "\n";  // DEBUG
   std::ifstream Inpfile(FileName);
   std::istream_iterator<double> my_it(Inpfile);
   std::cout << "it begins at " << *my_it << "\n";
   for (; my_it != std::istream_iterator<double>(); my_it++)
   {
	   singleton.first=(*my_it);
	   my_it++;
	   singleton.second=(*my_it);
	   li_reals.push_back( singleton );
           std::cout << singleton.first << " : " << singleton.second << "\n"; // DEBUG
	   
   }
	
	/* cpp reference
	std::istringstream str("0.1 0.2 0.3 0.4");
    std::partial_sum(std::istream_iterator<double>(str),
                     std::istream_iterator<double>(),
                     std::ostream_iterator<double>(std::cout, " "));
 
    std::istringstream str2("1 3 5 7 8 9 10");
    auto it = std::find_if(std::istream_iterator<int>(str2),
                      std::istream_iterator<int>(),
                      [](int i){return i%2 == 0;});
    if (it != std::istream_iterator<int>())
        std::cout << "\nThe first even number is " << *it << ".\n";
    //" 9 10" left in the stream
}
	*/
	
	// std::fclose(Inpfile);// is it needed
	return 0;
 }
 
 bool IDENT05_IODATA :: exportToDisk(li_doubles &li_reals) 
 {
	 FILE *fp = fopen(FileName, "wb");
	 int k;  // r1
	 k=0;
	 char* bufnum=NULL;

	 for (auto it=li_reals.begin(); it != li_reals.end(); it++) 
		 // li_doubles is typedef vector<dpair> hence an iterator can be used
	 {
		std::pair<double, double> singleton = (*it);  // TBRD
		strcpy((char *) bufnum, (std::to_string(singleton.first)).c_str());
	       fwrite(bufnum, sizeof(bufnum), 1, fp);
	       fwrite("\t", sizeof(char), 1, fp);
		strcpy((char *) bufnum, (std::to_string(singleton.second)).c_str());
	       fwrite(bufnum, sizeof(bufnum), 1, fp);
	    k++;
	 }
	 std::cout << "wrote " << k << " element lines of pair of data ";
	 std::fclose(fp);
	 return 0;
 }
