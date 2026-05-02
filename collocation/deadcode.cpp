   if (!(lh_code.eof()) && (row<1)) {
      std::string line;
      std::getline(lh_code, line);
      std::stringstream ss(line);  
      while(ss >> intarray[0][col]) col++;
      type_ovp=intarray[0][0];
      type_eqn=intarray[0][1];      
   }
   row++;col=0;

   if (!(lh_code.eof()) && (row<2)) {
      std::string line;
      std::getline(lh_code, line);
      std::stringstream ss(line);  
      while(ss >> intarray[0][col]) col++;
      type_predict=intarray[0][0];
      repeat_predict=intarray[0][1];      
   }
   row++;col=0;

   if (!(lh_code.eof()) && (row<3)) {
      std::string line;
      std::getline(lh_code, line);
      std::stringstream ss(line);  
      while(ss >> doublearray[0][col]) col++;
      boundry[0]=doublearray[0][0];
      boundry[1]=doublearray[0][1];      
   }   
   row++;col=0;

   if (!(lh_code.eof()) && (row<4)) {
      std::string line;
      std::getline(lh_code, line);
      std::stringstream ss(line);  
      while(ss >> doublearray[0][col]) col++;
      tinit=doublearray[0][0];
      tend=doublearray[0][1];      
   }
   row++;col=0;

   if (!(lh_code.eof()) && (row<5)) {
      std::string line;
      std::getline(lh_code, line);
      std::stringstream ss(line);  
      while(ss >> intarray[0][col]) col++;
      num_ranges=intarray[0][0];
      num_points=intarray[0][1];      
   }
   row++;col=0;

   if (!(lh_code.eof()) && (row<6)) {
      std::string line;
      std::getline(lh_code, line);
      std::stringstream ss(line);  
      while(ss >> doublearray[0][col]) col++;
      predictparams[0]=doublearray[0][0];
      predictparams[1]=doublearray[0][1];      
   }
   row++;col=0;


   /* read_parse_file() notdeadbut moved to another file */

bool IDENT05_COLL::read_parse_file() 
{
   Index row,col;

// Open the input file, which is a text file with 6 columns
  std::ifstream lh_file;
  lh_file.open("TheFile", std::ifstream::in);  // test w/o variable Name
  
  row=0;
  while ((!lh_file.eof()) && (row<41)) {
//  while (lh_file.good()) {
      std::string line; 
      std::getline(lh_file, line);    
    
      std::stringstream ss(line);
      col = 0;
      while(ss >> dataarray[row][col]) col++;
#ifdef __TEST_COLL_ONLY__
      std::cout << "tf= "  << dataarray[row][0] << " ";
      std::cout << "utf= " << dataarray[row][1] << " ";
      std::cout << "ytf= " << dataarray[row][2] << " ";
      std::cout << "c2f= " << dataarray[row][3] << " ";
      std::cout << "c1f= " << dataarray[row][4] << " ";
      std::cout << "c0f= " << dataarray[row][5] << "\n";
#endif
      row++;
  }
   num_rows=row;
   lh_file.close();

   return 0;
}


