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
static int g_initflag = 0;

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
    if (g_initflag == 0) {
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &MP_Np);
        MPI_Comm_rank(MPI_COMM_WORLD, &MP_Iam);
        g_initflag = 1;
    } else {
        return;
    }
}

void MP_abort()
{
    MPI_Finalize();
}

MP_buffer_t *MP_buffer_alloc(unsigned int size, data_type_t type)
{
    MP_buffer_t *p;

    p = new MP_buffer_t;
    p->N = size;
    if (type == CHAR) {
        p->bytes = size * sizeof(char);
    } else if (type == INT) {
        p->bytes = size * sizeof(int);
    } else if (type == LONG) {
        p->bytes = size * sizeof(long);
    } else if (type == FLOAT) {
        p->bytes = size * sizeof(float);
    } else if (type == DOUBLE) {
        p->bytes = size * sizeof(double);
    } else {
        p->bytes = 0;
        free(p);
        MP_ERROR("Unknown data type, can not alloc memory!");
        return NULL;
    }

    p->buff = new char[p->bytes];

    return p;
}

void MP_buffer_free(MP_buffer_t *buff)
{
    if (buff != NULL) {
        if (buff->buff != NULL)  free(buff->buff);
        free(buff);
    }
}

bool simple_send(int dest, MP_buffer_t *buff)
{
    if (buff->buff == NULL) {
        return false;
    }

    MPI_Send(buff->buff, buff->bytes, MPI_CHAR, dest, P2PTAG, MPI_COMM_WORLD);

    return true;
}

void simple_recv(MP_buffer_t *buff)
{
    if (buff->buff == NULL) {
        return;
    }

    MPI_Recv(buff->buff, buff->bytes, MPI_CHAR, MPI_ANY_TAG, P2PTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

