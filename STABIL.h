/*	STABIL.h
	
	See STABIL.c for an explanation.
	
	- Keshav Kini <kini@member.ams.org>, 2010-10-14
*/

#define EXIT_SUCCESS 0
#define EXIT_BAD_INPUT 1
#define EXIT_ALLOC_ERROR 2
#define EXIT_OVERFLOW 3

int STABIL(unsigned long* matrix, unsigned long n, unsigned long* d);
