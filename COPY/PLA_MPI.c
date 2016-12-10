/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/*-----------------------------------------------------------------------
*     File: PLA_MPI.c
*   Author: Mike Kistler
*  Descrip: PLAPACK wrappers for MPI routines
*
*  Note: Routines PLA_MPI_Reduce and PLA_MPI_Allreduce are defined in
*        $PLAROOT/REDUCE/PLA_Reduce.c
*----------------------------------------------------------------------*/

#include "PLA.h"

int PLA_MPI_Bcast(
	void *		buffer, 
	int 		count, 
	MPI_Datatype	datatype,
	int		root, 
	MPI_Comm	comm)
{
	int value;
	double Start_time;

	if ( PLA_TIMINGS ){
		Start_time = MPI_Wtime();
	}

    value = MPI_Bcast ( buffer, count, datatype, root, comm );

	if ( PLA_TIMINGS ){
		PLA_TIMINGS[PLA_MPI_BCAST_TIMING] += MPI_Wtime() - Start_time;
	}

	return value;
}


int PLA_MPI_Gatherv(
	void *		sendbuf, 
	int 		scount, 
	MPI_Datatype 	stype, 
	void *		recvbuf,
	int *		rcounts,
	int *		displs,
	MPI_Datatype 	rtype, 
	int 		root, 
	MPI_Comm	comm)
{
    return (MPI_Gatherv ( sendbuf, scount, stype, recvbuf, rcounts, displs, rtype, root, comm ));
}


int PLA_MPI_Scatterv(
	void *		sendbuf, 
	int *		scounts, 
	int *		displs,
	MPI_Datatype 	stype, 
	void *		recvbuf,
	int 		rcount,
	MPI_Datatype 	rtype, 
	int 		root, 
	MPI_Comm 	comm)
{
    return (MPI_Scatterv ( sendbuf, scounts, displs, stype, recvbuf, rcount, rtype, root, comm ));
}


int PLA_MPI_Allgatherv(		/* aka Collect */
	void *		sendbuf,
	int 		scount, 
	MPI_Datatype 	stype, 
	void *		recvbuf,
	int *		rcounts,
	int *		displs,
	MPI_Datatype 	rtype, 
	MPI_Comm 	comm)
{
    return (MPI_Allgatherv ( sendbuf, scount, stype, recvbuf, rcounts, displs, rtype, comm ) );
}


int PLA_MPI_Reduce_scatter(	/* aka Distributed Reduce */
	void *		sendbuf, 
	void *		recvbuf, 
	int *		rcounts,
	MPI_Datatype 	datatype,
	MPI_Op 		op, 
	MPI_Comm 	comm)
{
    return (MPI_Reduce_scatter ( sendbuf, recvbuf, rcounts, datatype, op, comm ));
}


int PLA_MPI_Send(
	void *		buf, 
	int 		count, 
	MPI_Datatype	datatype,
	int		    dest,
	int         tag,
	MPI_Comm	comm)
{
	int value;
	double Start_time;

	if ( PLA_TIMINGS ){
		Start_time = MPI_Wtime();
	}

    value = MPI_Send ( buf, count, datatype, dest, tag, comm );

	if ( PLA_TIMINGS ){
		PLA_TIMINGS[PLA_MPI_SEND_TIMING] += MPI_Wtime() - Start_time;
	}

	return value;
}


int PLA_MPI_Recv(
	void *		buf, 
	int 		count, 
	MPI_Datatype	datatype,
	int		    source,
	int         tag,
	MPI_Comm	comm,
	MPI_Status  *status)
{
	int value;
	double Start_time;

	if ( PLA_TIMINGS ){
		Start_time = MPI_Wtime();
	}

    value = MPI_Recv ( buf, count, datatype, source, tag, comm, status );

	if ( PLA_TIMINGS ){
		PLA_TIMINGS[PLA_MPI_RECV_TIMING] += MPI_Wtime() - Start_time;
	}

	return value;
}


int PLA_MPI_Irecv(
	void *		buf, 
	int 		count, 
	MPI_Datatype	datatype,
	int		    source,
	int         tag,
	MPI_Comm	comm,
	MPI_Request *request)
{
	int value;
	double Start_time;

	if ( PLA_TIMINGS ){
		Start_time = MPI_Wtime();
	}

    value = MPI_Irecv ( buf, count, datatype, source, tag, comm, request );

	if ( PLA_TIMINGS ){
		PLA_TIMINGS[PLA_MPI_IRECV_TIMING] += MPI_Wtime() - Start_time;
	}

	return value;
}


int PLA_MPI_Wait(
	MPI_Request *request,
	MPI_Status  *status)
{
	int value;
	double Start_time;

	if ( PLA_TIMINGS ){
		Start_time = MPI_Wtime();
	}

    value = MPI_Wait ( request, status );

	if ( PLA_TIMINGS ){
		PLA_TIMINGS[PLA_MPI_WAIT_TIMING] += MPI_Wtime() - Start_time;
	}

	return value;
}
