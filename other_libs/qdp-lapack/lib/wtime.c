/* $Id: wtime.c,v 1.2 2008-04-22 03:53:05 kostas Exp $ */
 /*
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
 *
 *   This file is part of PRIMME.
 *
 *   PRIMME is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   PRIMME is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <stdlib.h>
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#  include <sys/time.h>
#  include <sys/resource.h>
#endif

/* 
 * Other timers that may be of use -------------------------------------------
 */

/* Simply return the microseconds time of day */
double primme_get_wtime() {
   static struct timeval tv;

   gettimeofday(&tv, NULL);
   return ((double) tv.tv_sec) + ((double) tv.tv_usec ) / (double) 1E6;
}
