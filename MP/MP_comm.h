/*
 * MP_comm.h
 *
 *  Created on: 2015.5.20
 *      Author: Pingzhou
 */

#ifndef MP_COMM_H_
#define MP_COMM_H_

extern int MP_Np;
extern int MP_Iam;

double MP_now();
void MP_get_info(int &rank, int &proc);
void MP_init();
void MP_abort();

#endif /* MP_COMM_H_ */
