/*
 * testlibpq.c
 *
 *      Test the C version of libpq, the PostgreSQL frontend library.
 */
#include <stdio.h>
#include <stdlib.h>
#include <libpq-fe.h>

//#include "lh_solver_postgresql.h"
// include itself header file for the libpq C  library


// the main program here only for developppement
#ifndef __LH_POSTGRESQL_TEST__
#define __LH_POSTGRESQL_TEST__
#endif

#ifdef __LH_POSTGRES_TEST__
int main() {
PGconn * this_conn;
PostgresPollingStatusType status;
PQconninfoOption * info;
int i; 

const char * db_cmd;
PGresult * q_result;
char * columnName1, columnName2;
char * column1Value;
long int * column2Value;

	/*  establish a connection to Variables.pgsql , host/hostaddr is omitted s.t. attempst to connect to localhost*/
   this_conn *PQconnectStartParams("port", "dbname", "user", "passwd",
                             "5432", "postgres", "postgres", "lpf6lpgres",
                             );
   if this_conn == NULL {
      printf("Database Connection Start failed");
      exit(1);
   }
   
   status = PQconnectPoll(this_conn);
   if status == CONNECTION_MADE {
      printf("Connected to Server..");
   }
   /* get the infos used by this live connection */
   info = PQconninfo(this_conn);
   
   db_cmd="SELECT * FROM modelparams"
   q_result=PQexec(this_conn, db_cmd);
   if q_result == NULL {
      printf("query has failed");
      exit(1);
   }
   /* retrieve the results data from the last query and print the First two columns */
   columnName1=PQfname(q_result, 0);
   printf(columnName1);
   columnName2=PQfname(q_result, 1);
   printf(columnName2);
   for (i=0;i<10;i++) {
      column1Value=PQgetvalue(q_result,i,0);
      if column1Value == NULL {
         printf("Col 1 Row %i is not good retrieved", i);
	 exit(1);
      }
      printf("%s", column1Value); //assumption string 
      column2Value=PQgetvalue(q_result,i,1);
      if column2Value == NULL {
         printf("Col 2 Row %i is not good retrieved", i);
	 exit(1);
      printf("%d", column2Value); 
      printf("\n");
   }
   /* free memory */    
   PQclear();  // do not know if it has an argument
   
   /* closes the connection to the db */
   PQfinish(this_conn);
}

#endif
