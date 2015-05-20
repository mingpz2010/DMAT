/*
 * MP_comm.h
 *
 *  Created on: 2015.5.20
 *      Author: Pingzhou
 */

#ifndef MP_COMM_H_
#define MP_COMM_H_

#define MP_ERROR(x)     std::cerr<< x << std::endl

enum DATA_TYPE {
    CHAR = 1,
    INT,
    LONG,
    FLOAT,
    DOUBLE
}data_type_t;

typedef struct {
    data_type_t type;
    long N;
    long bytes;
    char *buff;
}MP_buffer_t;

#define P2PTAG      8521

extern int MP_Np;
extern int MP_Iam;

double MP_now();
void MP_get_info(int &rank, int &proc);
void MP_init();
void MP_abort();

MP_buffer_t *MP_buffer_alloc(unsigned int size, data_type_t type);
void MP_buffer_free(MP_buffer_t *buff);
bool simple_send(int dest, MP_buffer_t *buff);
void simple_recv(MP_buffer_t *buff);

void MP_bandwidth_testing();

#endif /* MP_COMM_H_ */
