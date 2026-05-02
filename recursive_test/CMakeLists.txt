cmake_minimum_required(VERSION 3.23)

project(lh_library_relsq)

add_executable(test_relsq)

target_sources(test_relsq
   PRIVATE
      ../recursive_test/test_relsq.cpp
)

add_library(LHRecursive)

target_sources(LHRecursive
   PRIVATE
      ../recursive_LSQ/impl_relsq.cpp

   PUBLIC
      FILE_SET LHRecHeaders
      TYPE HEADERS
      BASE_DIRS
         ../include
         ../recursive_LSQ/ 
      FILES
         ../include/ident05_colloc.hpp 
         ../recursive_LSQ/decl_relsq.h
)

target_link_libraries(test_relsq.exe
   PRIVATE
      LHRecursive
)


