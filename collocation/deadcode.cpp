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
