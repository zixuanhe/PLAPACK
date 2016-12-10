/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static double read_time=0.0;
static double write_time=0.0;

static int read_size=0;
static int write_size=0;

static int collect_stats=FALSE;

int PLA_File_stats_start()
{

  read_time = 0.0;
  write_time = 0.0;
  read_size = 0;
  write_size = 0;
  collect_stats = TRUE;

  return PLA_SUCCESS;
}
int PLA_File_stats_read(int size, double time)
{

  if(collect_stats ) {
    read_size += size;
    read_time += time;
  }
  return PLA_SUCCESS;
}

int PLA_File_stats_write(int size, double time)
{

  if(collect_stats) {
    write_size += size;
    write_time += time;
  }
  return PLA_SUCCESS;

}

int PLA_File_stats_stop()
{
  
  
  int
    me,
    total_read_size,
    total_write_size;

  double
    total_read_time,
    total_write_time;
 
  collect_stats = FALSE;

  /* Get the total size */
  MPI_Reduce(&read_size, &total_read_size, 1, MPI_INT, MPI_SUM,
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&write_size, &total_write_size, 1, MPI_INT, MPI_SUM,
	     0, MPI_COMM_WORLD);

  /* Get the largest time */
  MPI_Reduce(&read_time, &total_read_time, 1, MPI_DOUBLE, MPI_MAX,
	     0, MPI_COMM_WORLD);
  MPI_Reduce(&write_time, &total_write_time, 1, MPI_DOUBLE, MPI_MAX,
	     0, MPI_COMM_WORLD);


  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  if(me == 0) {
    if(total_read_time > 0.0)
      printf("READ : size(MB)=%4.2lf  time(s)=%4.2lf MB/s=%4.2lf\n", 
	     (double)total_read_size/1000000.0, 
	     total_read_time, 
	     ((double)total_read_size/1000000.0)/total_read_time);
    
    if(total_write_time > 0.0)
      printf("WRITE: size(MB)=%4.2lf time(s)=%4.2lf MB/s=%4.2lf\n", 
	     (double)total_write_size/1000000.0, 
	     total_write_time,
	     ((double)total_write_size/1000000.0)/total_write_time);
    
    if(total_read_time+total_write_time > 0.0)
      printf("BOTH : size(MB)=%4.2lf time(s)=%4.2lf MB/s=%4.2lf\n", 
	     (double)(total_read_size+total_write_size)/1000000.0,
	     total_read_time+total_write_time,
	     ((double)(total_read_size+total_write_size)/1000000.0)/(total_read_time+total_write_time));
  }
}

int PLA_File_stats_local_print()
{
  
  if(read_time > 0.0)
    printf("READ : size(MB)=%4.2lf time(s)=%4.2lf MB/s=%4.2lf\n", 
	   (double)read_size/1000000.0, read_time, 
	   ((double)read_size/1000000.0)/read_time);
  
  if(write_time > 0.0)
    printf("WRITE: size(MB)=%4.2lf time(s)=%4.2lf MB/s=%4.2lf\n", 
	   (double)write_size/1000000.0, write_time,
	   ((double)write_size/1000000.0)/write_time);
  
  if(read_time+write_time > 0.0)
    printf("BOTH : size(MB)=%4.2lf time(s)=%4.2lf MB/s=%4.2lf\n", 
	   (double)(read_size+write_size)/1000000.0,
	   read_time+write_time,
	   ((double)(read_size+write_size)/1000000.0)/(read_time+write_time));
  

}
