#ifndef AuxCortexNematic2D
#define AuxCortexNematic2D

#include <iostream>
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <vector>
using namespace std;



void splittingCommunicators()
{


  int color, size, rank, subRank;

  MPI_Comm new_comm;

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank%2==0)
      color = 1;
  else
      color = 2;


  MPI_Comm_split(MPI_COMM_WORLD,color,rank,&new_comm);


  MPI_Comm_rank(new_comm,&subRank);


  printf("I am rank in %d in MPI_COMM_WORLD and %d in new_comm. \r\n ", rank,subRank);


}                                                                                             



struct Partstruct
{
    int c;
    double d[6];
    char b[7];
};

void nonBlockingCommunication()
{

    MPI_Request request;
    MPI_Status  status;
    int rank, size;
    int buffer_count=100;
    int buffer[buffer_count];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double deltat_p0;
    if (rank == 0)
    {

        double start = MPI_Wtime();
        for (int i=0; i<buffer_count ; i++)
            buffer[i] = i;


        MPI_Isend(buffer, buffer_count, MPI_INT, 1, 666,MPI_COMM_WORLD,&request);
    
        // Performed 5 second long computations by processor 0
        sleep(5);
        deltat_p0 = MPI_Wtime()- start;

        cout<<"NON-BLOCKING SETTING::: Time take by the processor 0 to send buffer and perform computation is "<<deltat_p0<<"seconds"<<endl;

    }
    else if (rank == 1)
    {

        //Processor 1  is going to performs 7 second of computation before its receive entires of the buffer
        sleep(7);

        MPI_Irecv(buffer, buffer_count, MPI_INT, 0, 666,MPI_COMM_WORLD, &request);
        
        MPI_Wait(&request, &status);

    }

    MPI_Barrier(MPI_COMM_WORLD); 
}


void BlockingCommunication()
{

    MPI_Request request;
    MPI_Status  status;
    int rank, size;
    int buffer_count=100;
    int buffer[buffer_count];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double deltat_p0;
    if (rank == 0)
    {

        double start = MPI_Wtime();
        for (int i=0; i<buffer_count ; i++)
            buffer[i] = i;

        MPI_Send(buffer, buffer_count, MPI_INT, 1, 666,MPI_COMM_WORLD);

        // Performed 5 second long computations by processor 0
        sleep(5);
        deltat_p0 = MPI_Wtime()- start;

        cout<<"BLOCKING SETTING::: Time take by the processor 0 to send buffer and perform computation is "<<deltat_p0<<"seconds"<<endl;

    }
    else if (rank == 1)
    {
        //Processor 1  is going to performs 7 second of computation before its receive entires of the buffer
        sleep(7);

        MPI_Recv(buffer, buffer_count, MPI_INT, 0, 666,MPI_COMM_WORLD, &status);
        MPI_Wait(&request, &status);        
    }

     MPI_Barrier(MPI_COMM_WORLD); 

}


void nonBlockingCommunicationWthRaceCondition()
{

    MPI_Request request;
    MPI_Status  status;
    int rank, size;
    int buffer_count=10;
    int buffer[buffer_count];
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double deltat_p0;
    if (rank == 0)
    {


        for (int i=0; i<buffer_count ; i++)
            buffer[i] = i;

        MPI_Isend(buffer, buffer_count, MPI_INT, 1, 666,MPI_COMM_WORLD,&request);

        // The array which was supposed to be sent
        printf("The sent entries of the array buffer by precoessor 0 are \r\n");
        for (int i = 0; i < buffer_count; i++)
            cout << buffer[i]<<"  ";
        printf("\r\n");

        // Processor 0 does rest of the computation
        sleep(5);

    }
    else if (rank == 1)
    {

        //Initializing the buffer array for processor 1
        for (int i=0; i<buffer_count ; i++)
            buffer[i] = 0;

        //Processor 1  is going to perform 3 second of computation before its receive entires of the buffer
        sleep(3);

        MPI_Irecv(buffer, buffer_count, MPI_INT, 0, 666,MPI_COMM_WORLD, &request);


        // The buffer array received by the processor 1
        printf("The received entries of the array buffer by precoessor 1 are \r\n");
        for (int i = 0; i < buffer_count; i++)
            cout << buffer[i]<<"  ";
        printf("\r\n");

        printf("Use MPI_WAIT( ) right after MPI_Irecv( ) to avoid the race condition \r\n"); 
        MPI_Wait(&request, &status);

    }

}

void CommunicationWithMPIProbe()
{

    MPI_Request request;
    MPI_Status  status;
    int rank, size;
    int buffer_count=99;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double deltat_p0;
    if (rank == 0)
    {
        int tag = 666;
        int buffer[buffer_count];
        for (int i=0; i<buffer_count ; i++)
            buffer[i] = i;

        MPI_Send(buffer, buffer_count, MPI_INT, 1, tag,MPI_COMM_WORLD);

    }
    else if (rank == 1)
    {

        int buffer_count_received;
        MPI_Status status;
        int tag = 666;
        MPI_Probe(0, tag, MPI_COMM_WORLD, &status);

        // From the probed status we get the number of elements to receive
        MPI_Get_count(&status, MPI_INT, &buffer_count_received);
        int buffer[buffer_count_received];
        printf("The received size of the buffer is %d",buffer_count_received);
        MPI_Recv(buffer, buffer_count_received, MPI_INT, 0, 666,MPI_COMM_WORLD, &status);

    }

}



void communicateMultiDimensionalArray()
{

    MPI_Request request;
    MPI_Status  status;
    int rank, size, tag=666;
    int buffer_count_columns=10;
    int buffer_count_rows=5;
    int buffer[buffer_count_rows][buffer_count_columns];

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double deltat_p0;
    if (rank == 0)
    {

        for (int i=0; i<buffer_count_rows ; i++)
        {
            for (int j=0; j<buffer_count_columns ; j++)
            {
                buffer[i][j] = i*j;
            }
        }


        MPI_Send(buffer, buffer_count_rows*buffer_count_columns, MPI_INT, 1, tag,MPI_COMM_WORLD);




    }
    else if (rank == 1)
    {

        sleep(5);

        MPI_Recv(buffer, buffer_count_rows*buffer_count_columns, MPI_INT, 0, tag,MPI_COMM_WORLD, &status);

    }


}



void serialTrapazoidalRule(int index_init, int index_final, double h, double a, int rank, double *sum_local)
{

    for (int i=index_init; i<index_final; i++)
    {
        *sum_local = *sum_local + h * sin(a + i * h); // sin(x) is the function to be integrated here
    }
    if (rank==0)
    {
        *sum_local = *sum_local + h * (sin(90. * M_PI / 180.) + sin(180. * M_PI / 180.)) / 2.0;
         printf("The local sum is %lf \r\n",*sum_local);

    }
}


void parallelTrapazoidalRule(int rank, int size)
{

//    int rank, size;
//    MPI_Init(NULL, NULL);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &size);


    //total number of points. Here each processor will obtain 100 sample points
    int  n = 1E9; //1E10 number of sample points
    // Initial index of the processor "rank", final index of the processor "rank"
    int index_init, index_final;
    //Local contribution of of calculating the integral by each processor
    double sum_local{};
    //  Limits of integration [90,180] degrees
    double a = 90.*M_PI/180.,b=180.0*M_PI/180.;

    // Grid size
    double h  = (b-a)/n;

    // Work distrbution between the processors start from here
    int n_p = floor( ( (n-1)/size ));
    if (rank< (n-1)%size)
        n_p = n_p+1;
   
     
    index_init= 1+ rank*floor( (double) (n-1)/ (double) size ) + min(rank,  (n-1)%size ) ;

    index_final =  index_init + n_p - 1;                                                                

    // Each processor performs as serial computation
    serialTrapazoidalRule(index_init, index_final, h, a, rank, &sum_local);

    // Now send the answer of the numerical integragtion to the root node 0
    double global_sum{};

    MPI_Reduce(&sum_local, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

    
 //   if (rank==0)
 //       cout<<global_sum<<endl;

    return ;


}

void communicateBCast(int argc, char **argv, bool diagnosticFlag)
{

   MPI_Init(&argc, &argv);
   MPI_Request request;
   MPI_Status  status;    
   int size, rank;  
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   int   n_elements = 100;  
   double buffer[n_elements];
  
   if (rank==0)
   {
 
      for (int i = 0 ; i< n_elements; i++)
      { 
        buffer[i] =  i;    
        
      }                        
        
   }  
   

   MPI_Bcast(&buffer,100,MPI_DOUBLE,0,MPI_COMM_WORLD);

   //cout<<buffer[50]<<endl;

/*
   //cout<<n_elements<<endl;  
   if (rank!=0)
   {
       buffer = new double[n_elements];
   }

   for (int i = 0 ; i< n_elements; i++)
          buffer[i] = i;                     
   

   for (int i = 0 ; i< n_elements; i++)
       cout<<buffer[i]<<endl; 

   
   MPI_Bcast(&buffer,100,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   cout<<n_elements<<endl;
    

 //  MPI_Barrier(MPI_COMM_WORLD); 

   if (rank==0)
   { 
       cout<<"proloterian"<<buffer[99]<<endl;                    
   } 


*/

    

}


void communicateAll_Reduce(double globalCoords[16][2])
{

   MPI_Request request;
   MPI_Status  status;    
   int size, rank;  
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
   int n = 16; // size of the global coordinates   
   int index_init, index_final;
   int n_p = floor(  n/size );
   if (rank< n%size)
       n_p = n_p+1;
   
   
   index_init= rank*floor(n/size) + min(rank,n%size);

   index_final =  index_init + n_p - 1;                                                        

   double localCoords[n_p][2];
  
   int local_idx = 0 ;

   double sumLocalCoords[2]{}; 
   double sumGlobalCoords[2]{};
   
   cout<<index_init<<" "<<index_final<<" "<<rank<<endl; 


   for (int global_idx = index_init ; global_idx<index_final+1; global_idx++)
   { 
       localCoords[local_idx][0]=   globalCoords[global_idx][0]; 
       localCoords[local_idx][1]=   globalCoords[global_idx][1];
       sumLocalCoords[0] += localCoords[local_idx][0];
       sumLocalCoords[1] += localCoords[local_idx][1];    
       local_idx++;
   }

   
   MPI_Barrier(MPI_COMM_WORLD);
   // cout<<"hello shafooo"<<endl;   
   MPI_Allreduce(&sumLocalCoords, &sumGlobalCoords, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
   // Coordinates of the barycenter are:
   sumGlobalCoords[0] = sumGlobalCoords[0]/n;
   sumGlobalCoords[1] = sumGlobalCoords[1]/n;
    
   
}


void communicateMPIScatterGather()
{

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int globaldata[size];/*wants to declare array this way*/
    int localdata;/*without using pointers*/

    int i;
    if (rank == 0) {

        for (i=0; i<size; i++)
            globaldata[i] = i;

        printf("1. Processor %d has data: ", rank);
        for (i=0; i<size; i++)
            printf("%d ", globaldata[i]);
        printf("\n");
    }

    MPI_Scatter(globaldata, 1, MPI_INT, &localdata, 1, MPI_INT, 0, MPI_COMM_WORLD);

    printf("2. Processor %d has data %d\n", rank, localdata);
    localdata= 5;
    printf("3. Processor %d now has %d\n", rank, localdata);

    MPI_Gather(&localdata, 1, MPI_INT, globaldata, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("4. Processor %d has data: ", rank);
        for (i=0; i<size; i++)
            printf("%d ", globaldata[i]);
        printf("\n");
    }


   // MPI_Finalize();
    return;



}

 

void communicateMultiB_SEND(int argc, char **argv, bool diagnosticFlag)
{

}

void communicateMultiR_SEND(int argc, char **argv, bool diagnosticFlag)
{

}

void communicateMultiS_SEND(int argc, char **argv, bool diagnosticFlag)
{

}



void communicateRMA_(int argc, char** argv, bool diagnosticFlag)
{

    int NUM_ELEMENT=4;

    int i, id, num_procs, len, localbuffer[NUM_ELEMENT], sharedbuffer[NUM_ELEMENT];

    char name[MPI_MAX_PROCESSOR_NAME];

    MPI_Win win;

	 

    MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Get_processor_name(name, &len);

    printf("Rank %d running on %s\n", id, name);
    MPI_Win_create(sharedbuffer, NUM_ELEMENT, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

	 

    for (i = 0; i < NUM_ELEMENT; i++)
    {
       sharedbuffer[i] = 10*id + i;
       localbuffer[i] = 0;
    }

	 

	   printf("Rank %d sets data in the shared memory:", id);

	 

	   for (i = 0; i < NUM_ELEMENT; i++)
	      printf(" %02d", sharedbuffer[i]);

	 

	   printf("\n");

	  

	   MPI_Win_fence(0, win);

	 
//	   if (id != 0)
//	      MPI_Get(&localbuffer[0], NUM_ELEMENT, MPI_INT, id-1, 0, NUM_ELEMENT, MPI_INT, win);
//	   else
	   MPI_Get(&localbuffer[0], NUM_ELEMENT, MPI_INT, num_procs-1, 0, NUM_ELEMENT, MPI_INT, win);

	   

	   MPI_Win_fence(0, win);

	 

	   printf("Rank %d gets data from the shared memory:", id);

	 

	   for (i = 0; i < NUM_ELEMENT; i++)
	      printf(" %02d", localbuffer[i]);

	 

	   printf("\n");

	 

	   MPI_Win_fence(0, win);

	 

	   if (id < num_procs-1)
	      MPI_Put(&localbuffer[0], NUM_ELEMENT, MPI_INT, id+1, 0, NUM_ELEMENT, MPI_INT, win);
	   else
	      MPI_Put(&localbuffer[0], NUM_ELEMENT, MPI_INT, 0, 0, NUM_ELEMENT, MPI_INT, win);

	 

	   MPI_Win_fence(0, win);

	 

	   printf("Rank %d has new data in the shared memory:", id);

	  

	   for (i = 0; i < NUM_ELEMENT; i++)
	      printf(" %02d", sharedbuffer[i]);

	 

	   printf("\n");

	 

	   MPI_Win_free(&win);

	   return;

}



void communicateRMA(int argc, char** argv, bool diagnosticFlag)
{

    int  NROWS =  100;
    int  NCOLS =  100;
    int rank, nprocs, A[NROWS][NCOLS], i, j;
    MPI_Win win;
    MPI_Datatype column, xpose;
    int errs = 0;
    
    


    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (nprocs != 2) {
        printf("Run this program with 2 processes\n");fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
 
    if (rank == 0)
    {
        for (i=0; i<NROWS; i++)
            for (j=0; j<NCOLS; j++)
                A[i][j] = -1;
 
        /* create datatype for one column */
        MPI_Type_vector(NROWS, 1, NCOLS, MPI_INT, &column);
        /* create datatype for matrix in column-major order */
        MPI_Type_hvector(NCOLS, 1, sizeof(int), column, &xpose);
        MPI_Type_commit(&xpose);

        MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
 
        MPI_Win_fence(0, win);
 
        MPI_Get(A, NROWS*NCOLS, MPI_INT, 1, 0, 1, xpose, win);

        MPI_Type_free(&column);
        MPI_Type_free(&xpose);
 
        MPI_Win_fence(0, win);
 
        for (j=0; j<NCOLS; j++)
        {
            for (i=0; i<NROWS; i++)
            {
                if (A[j][i] != i*NCOLS + j)
                {
                    if (errs < 50)
                    {
                        printf("Error: A[%d][%d]=%d should be %d\n", j, i, A[j][i], i*NCOLS + j);fflush(stdout);
                    }
                    errs++;
                }
            }
        }
        if (errs >= 50)
        {
            printf("Total number of errors: %d\n", errs);fflush(stdout);
        }
    }
    else
    { /* rank = 1 */
        for (i=0; i<NROWS; i++)
            for (j=0; j<NCOLS; j++)
                A[i][j] = i*NCOLS + j;
        MPI_Win_create(A, NROWS*NCOLS*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
        MPI_Win_fence(0, win);
        MPI_Win_fence(0, win);
    }
 
    MPI_Win_free(&win);
    MPI_Finalize();
}


void  mpi_custom_types(int argc, char** argv, bool diagnosticFlag)
{

    struct Partstruct particle[2];
    int i, j, myrank;
    MPI_Status status;
    MPI_Datatype Particletype;
    MPI_Datatype type[3] = { MPI_INTEGER, MPI_DOUBLE, MPI_CHAR };
    int blocklen[3] = { 1, 6, 7 };
    /* compute displacements of structure components */
    MPI_Aint disp[3] = { offsetof(Partstruct, c), offsetof(Partstruct, d), offsetof(Partstruct, b) };
    //The following is the longer way of calculating the displacements
    //  MPI_Get_address(particle, disp);
    //  MPI_Get_address(particle[0].d, disp+1);
    //  MPI_Get_address(particle[0].b, disp+2);
    //  MPI_Aint base = disp[0];
    //  for (i=0; i < 3; i++) disp[i] = MPI_Aint_diff(disp[i], base);

    
    MPI_Type_create_struct(3, blocklen, disp, type, &Particletype);
    MPI_Type_commit(&Particletype);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0)
    {
        //Populating the struct particle
        particle[0].c = 100; particle[1].c = 200;
        particle[0].d[0] = 111.1;  particle[0].d[1] = 222.2;  particle[0].d[2] = 333.3;  particle[0].d[3] = 444.4;  particle[0].d[4] = 555.5;  particle[0].d[5] = 666.6;
        particle[1].d[0] = 111.1;  particle[1].d[1] = 222.2;  particle[1].d[2] = 333.3;  particle[1].d[3] = 444.4;  particle[1].d[4] = 555.5;  particle[1].d[5] = 666.6;
        particle[0].b[0] = 'a';  particle[0].b[1] = 'b';  particle[0].b[2] = 'c';  particle[0].b[3] = 'd';  particle[0].b[4] = 'e';  particle[0].b[5] = 'f';   particle[0].b[6] = 'g';
        particle[1].b[0] = 'h';  particle[1].b[1] = 'i';  particle[1].b[2] = 'j';  particle[1].b[3] = 'k';  particle[1].b[4] = 'l';  particle[1].b[5] = 'm';   particle[0].b[6] = 'n';

        // Sending  the struct particle
        MPI_Send(particle, 2, Particletype, 1, 123, MPI_COMM_WORLD);
     }
    else if (myrank == 1)
    {
        MPI_Recv(particle, 2, Particletype, 0, 123, MPI_COMM_WORLD, &status);

        cout<<"Few entries of member array d are: "<<particle[0].d[0]<< " "<<particle[0].d[1]<< " "<< particle[0].d[2]<<endl;

        cout<<"Entries of member array c are: "<<particle[0].c<<" "<<particle[1].c<<endl;

        cout<<"Few entries of member array v are: "<<particle[1].b[0]<< " "<<particle[1].b[1]<< " "<< particle[1].b[2]<<endl;

    }
    return;


}

void communicate_graph_topology(int argc, char** argv, bool diagnosticFlag)
{


    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nnode = 4;
    int source = rank;
    int degree;
    int dest[nnode];
    int weight[nnode] = {1, 1, 1, 1};
    int recv[nnode] = {-1, -1, -1, -1};
    int send = rank;

    // set dest and degree.
    if (rank == 0)
    {
        dest[0] = 1;
        dest[1] = 3;
        degree = 2;
    }
    else if(rank == 1)
    {
        dest[0] = 0;
        degree = 1;
    }
    else if(rank == 2)
    {
        dest[0] = 3;
        dest[1] = 0;
        dest[2] = 1;
        degree = 3;
    }
    else if(rank == 3)
    {
        dest[0] = 0;
        dest[1] = 2;
        dest[2] = 1;
        degree = 3;
    }

    // create graph.
    MPI_Comm graph;
    MPI_Dist_graph_create(MPI_COMM_WORLD, 1, &source, &degree, dest, weight, MPI_INFO_NULL, 1, &graph);

    // send and gather rank to/from neighbors.
    MPI_Neighbor_allgather(&send, 1, MPI_INT, recv, 1, MPI_INT, graph);

    printf("Rank: %i, recv[0] = %i, recv[1] = %i, recv[2] = %i, recv[3] = %i\n", rank, recv[0], recv[1], recv[2], recv[3]);

    MPI_Finalize();
    return;


}
#endif
