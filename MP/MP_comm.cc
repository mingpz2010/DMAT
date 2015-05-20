/*
 * MP_comm.cc
 *
 *  Created on: 2015.5.20
 *      Author: Pingzhou
 */
#include <iostream>
#include <mpi.h>
#include "MP_comm.h"

int MP_Np = -1;
int MP_Iam = 0;

double MP_now()
{
    return MPI_Wtime();
}

void MP_get_info(int &rank, int &proc)
{
    rank = MP_Iam;
    proc = MP_Np;
}

void MP_init()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &MP_Np);
    MPI_Comm_rank(MPI_COMM_WORLD, &MP_Iam);
}

void MP_abort()
{
    MPI_Finalize();
}

