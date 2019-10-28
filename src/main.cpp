
/// C++ headers
#include <iostream>
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
using namespace std;
#include "auxiliary.h"


int main(int argc, char **argv)
{
     
     // test case to demonstrate how to split communicators and use them for point-tp- point communication
    bool testCase_0 = false;
     // test Case to demonstrate the overhead introduced by the blocking communication
    bool testCase_1 = false; 
    // test case to demonstrate the race condition
    bool testCase_2 = false;  
    // test case to demonstrate the significance of mpi probe
    bool testCase_3 = false;  

    bool testCase_4 = false ; // test Case to demonstrate how to pass a multi-dimensional array using the blocking communication
    bool testCase_5 = false;  // test Case to demonstrate the difference between the performance of RSend( ), BSend( ) and SSend( ) routines

    // test Case to demonstrate the speed up obtain using multiple cores. The test case involves calculating nuemrical integration
    bool testCase_6 = false;   
    bool testCase_7 = false;   // test Case to demonstrate the functionality of BCast subroutine.
// test Case to demonstrate the functionality of MPI_Allreduce(). The description of the test case can be seen here: https://www.codingame.com/playgrounds/349/introduction-to-mpi/reduction---exercise-2
    bool testCase_8 = false;    

    // test Case to deminstrate the functionality of MPI_Scatter and MPI_Gather subroutines.
    bool testCase_9 = false;  
     // test Case to demonstrate the RMA model  
    bool testCase_10= false; 
     // test Case to demonstrate an example on the custom types
    bool testCase_11= false ;
     // test case to demonstrate building intra-node communication vith graph topologies 
    bool testCase_12= true;  

    int rank,size;

    MPI_Init(nullptr, nullptr);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  
    if (testCase_0)
    {
        if (size%2==0)  
        {   
            if (rank==0) printf("ERROR:: Please try with odd number of precessors \r\n");
        } 
        else 
            splittingCommunicators();      
 
    }
    if (testCase_1)
    {
        // The following two routines estimates the efficiency of non-blocking communication in terms of time in seconds.
        if (size!=2)
        {
            if (rank==0) printf("ERRRO:: Please try with 2 processors");
        }  
        else 
        {
            nonBlockingCommunication();
            BlockingCommunication();   
        }    
    }


    if (testCase_2)
    {
        // The following subroutines estimates the importance of MPI_TEST and MPI_WAIT (RACE conditions)
       if (size!=2)
       {
           if (rank==0) printf("ERRRO:: Please try with 2 processors");
       }
       else  
           nonBlockingCommunicationWthRaceCondition();
    }


    if (testCase_3)
    {

       // The following subroutine demonstrates the importance of MPI_Probe( ) in point to point nodal communication.
       if (size!=2)
       {
           if (rank==0) printf("ERRRO:: Please try with 2 processors");
       }
       else                                                                  
           CommunicationWithMPIProbe();
    }


    if (testCase_4)
    {
        // The following subroutine demonstrates  how to communicate  a multi-dimesnional array between the processor
        if (size!=2)
        {
            if (rank==0) printf("ERRRO:: Please try with 2 processors");
        }
        else                                                                  
            communicateMultiDimensionalArray();
    }

    if (testCase_5)
    {
        MPI_Init(nullptr, nullptr);
        double average_speed_up{}, total_iter=100;




        for (int index = 0; index<total_iter;index++)
        {

            double global_sum_serial{}, global_sum_parallel{}, total_time_serial{},total_time_parallel{};

            // Start time

            double start_parallel = MPI_Wtime();

            parallelTrapazoidalRule(rank,size);

            total_time_parallel = MPI_Wtime() - start_parallel;

            MPI_Barrier(MPI_COMM_WORLD);

            if (rank == 0)
            {

                double start_serial = MPI_Wtime();
                //total number of points. Here each processor will obtain 100 sample points
                int  n = 1E9; //1E9number of sample points
                // Initial index of the processor "rank", final index of the processor "rank"
                int index_init=1, index_final=1E9-1;
                //Local contribution of of calculating the integral by each processor
                double sum_local{};
                //  Limits of integration [90,180] degrees
                double a = 90.*M_PI/180.,b=180.0*M_PI/180.;

                // Grid size
                double h  = (b-a)/n;

                serialTrapazoidalRule(index_init, index_final, h, a, rank, &global_sum_serial);
                total_time_serial = MPI_Wtime() - start_serial;
                average_speed_up += total_time_serial/total_time_parallel;
            }

             if (rank == 0)
                printf("The average speed up of this test is %lf in iteration %d \r\n ", total_time_serial/total_time_parallel, index);

 
        }
        // Print the result
        if (rank == 0)
            printf("The average speed up of this test is %lf \r\n ", (double) average_speed_up/ (double) total_iter);

    }
    
    if (testCase_8)
    {   
        // In this test case we calculate the barycenter of the following coordinates and then communicate this information to all processors 
        double  coords[16][2] = { {0.,0.}, {0.1,0.}, {0.2,0.}, {0.3,0.},
                                 {0.,0.1}, {0.1,0.1}, {0.2,0.1}, {0.3,0.1}, 
                                 {0.,0.2}, {0.1,0.2}, {0.2,0.2}, {0.3,0.2},
                                 {0.,0.3}, {0.1,0.3}, {0.2,0.3}, {0.3,0.3}} ;    
         
        communicateAll_Reduce(coords);      
    }
   
    if (testCase_9)
    {
        communicateMPIScatterGather();
    }
 
    if (testCase_10)
    {     
       if (size!=2)     
       {
           if (rank==0) printf("ERROR:: Please try with 2  precessors \r\n");
       }
        else
           communicateRMA();
    }

         cout<<size<<endl;
    if (testCase_11)
    {
        if (size!=2) 
        {    
            if (rank==0) printf("ERROR:: Please try with 2  precessors \r\n"); 
        }
        else
            mpi_custom_types();
    }

   if (testCase_12)
   {

       if (size!=4)     
       {
           if (rank==0) printf("ERROR:: Please try with 4 precessors \r\n"); 
       }
       else
       {
           communicate_graph_topology();    
       }
   }



   MPI_Finalize();


    
}

