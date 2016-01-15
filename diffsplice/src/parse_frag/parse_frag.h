/*    
 *    parseFragment.cpp		
 *    DiffSplice
 *
 *    Copyright (C) 2012 University of Kentucky and
 *                       Yin Hu
 *
 *    Authors: Yin Hu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */



/************************************************************************/
/* Get Fragments From SAM File                                          */
/************************************************************************/

//Should be unique alignments

#define UNIX

#ifdef UNIX
#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib> 
#include <unistd.h>
#include <sstream>
#include <vector>
#else
#include <fstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h> 
#include <sstream>
#include <vector>
#endif

using namespace std;

//#define DEBUG

const long MAX = 2000000000;

ofstream statFile;

class record
{
public:
	long occurrence;
	long start_pos;
	long end_pos;
	int xs;

	record(long occur = 0, long startp = 0, long endp = 0, int xs_field = 0) {occurrence = occur; start_pos = startp; end_pos = endp; xs = xs_field;};
	int rec_compare(const record &ref_record, bool is_exon); //compare this record with the ref_record, return -1 if this record < ref_record, return 0 if equal, return 1 if this > ref_record
};

