inline double lh_atof(const char *str)
{
      while ((*str) && isspace(*str))
       ++str;
      return strcmp(str, "nan") ? atof(str) : NAN;
}

// Exercice: do a little bit differently for the next function
bool IDENT05_COLL::read_parse_specfile() {
	
   //std::ifstream lh_fspec;
   char head[256];
   std::string headstr;
   char flag[256];
   std::string flagstr;
   char fl[3];
   char value[256];
   std::string valuestr;

   std::cout << "First line " << headstr << "\n";
   strcpy(head, headstr.c_str());
   if (!strncmp(head, "INT:RNGORD0", 11)) {
      std::cout << "First line was read: \n";
      std::cout << head << "\n";
  //    lh_fspec.getline(flagstr >> std::ws, 256, '\n');
      getline(lh_fspec, flagstr);
  //    // SO: the std::ws IO Manipulator can be used to discard the leading whitespace
      strcpy(flag, flagstr.c_str());
      std::cout << "Second line was read: \n";
      std::cout << flagstr << "\n";
      strncpy(fl, (const char*) flag, 3);
      
      if (!strncmp(flag, "DOU", 3)) {
         std::cout << "A third line will be read \n";
      //   lh_fspec.getline(valuestr, 256, '\n');
         getline(lh_fspec, valuestr);
         strcpy(value, (const char*) valuestr.c_str());
         valued1 = lh_atof(value);

	 // valued1 << valuestr; this causes an error of compiling
         std::cout << "Parsed number read: " << valued1 << "\n" ;
      }
      else {
         std::cout << "The type of args does not match (DOU) \n";
      }
   }  
   lh_fspec.close();
   return 0;
}
