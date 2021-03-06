/*
 * Copyright (c), MM Weiss
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *     1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *     2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *     3. Neither the name of the MM Weiss nor the names of its contributors
 *     may be used to endorse or promote products derived from this software without
 *     specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 * OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 *  clock_gettime_stub.c
 *  gcc -Wall -c clock_gettime_stub.c
 *  user space glue
 */

#ifdef __APPLE__

#include "mac_clock_gettime.h"

int clock_gettime(clockid_t clk_id, struct timespec *tp) {
	kern_return_t   ret;
	clock_serv_t    clk;
	mach_timespec_t tm;

	uint64_t start;
	uint64_t end;
	uint64_t delta;
	uint64_t nano;
	static mach_timebase_info_data_t __clock_gettime_inf;

	int retval = -1;
	switch (clk_id) {
		case CLOCK_REALTIME:
			if (KERN_SUCCESS == (ret = host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &clk))) {
				if (KERN_SUCCESS == (ret = clock_get_time(clk, &tm))) {
					tp->tv_sec  = tm.tv_sec;
					tp->tv_nsec = tm.tv_nsec;
					retval = 0;
				}
			}
			if (KERN_SUCCESS != ret) {
				errno = EINVAL;
				retval = -1;
			}
			break;
		case CLOCK_MONOTONIC:
			if (KERN_SUCCESS == (ret = host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &clk))) {
				if (KERN_SUCCESS == (ret = clock_get_time(clk, &tm))) {
					tp->tv_sec  = tm.tv_sec;
					tp->tv_nsec = tm.tv_nsec;
					retval = 0;
				}
			}
			if (KERN_SUCCESS != ret) {
				errno = EINVAL;
				retval = -1;
			}
			break;
		case CLOCK_PROCESS_CPUTIME_ID:
		case CLOCK_THREAD_CPUTIME_ID:
			start = mach_absolute_time();
			if (clk_id == CLOCK_PROCESS_CPUTIME_ID) {
				getpid();
			} else {
				sched_yield();
			}
			end = mach_absolute_time();
			delta = end - start;
			if (0 == __clock_gettime_inf.denom) {
				mach_timebase_info(&__clock_gettime_inf);
			}
			nano = delta * __clock_gettime_inf.numer / __clock_gettime_inf.denom;
			tp->tv_sec = nano * 1e-9;
			tp->tv_nsec = nano - (tp->tv_sec * 1e9);
			retval = 0;
			break;
		default:
			errno = EINVAL;
			retval = -1;
	}
	return retval;
}

#endif // __APPLE__

/* EOF */
