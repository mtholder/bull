// This file is part of BULL, a program for phylogenetic simulations
// most of the code was written by Mark T. Holder.

//	This program is for internal use by the lab of Dr. Tandy Warnow only.
//	Do not redistribute the code.  It is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//
//	Some of the code is from publically available source by Paul Lewis, Ziheng Yang, 
//	John Huelsenbeck, David Swofford , and others  (as noted in the code).
//	In fact the main structure of the program was created by modifying Paul Lewis' 
//	basiccmdline.cpp from his NCL
//
//	This code was used in Mark's dissertation, some changes were made in order to 
//	get it to compile on gcc.  It is possible that this porting introduced bugs (very
//	little debugging has been done on UNIX platforms).	I would suggest checking 
//	the simulator by generating data on trees with short branches, etc.
 

// This file is part of BULL, a program for phylogenetic simulations
// most of the code was written by Mark T. Holder.

//	This program is for internal use by the lab of Dr. Tandy Warnow only.
//	Do not redistribute the code.  It is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//
//	Some of the code is from publically available source by Paul Lewis, Ziheng Yang, 
//	John Huelsenbeck, David Swofford , and others  (as noted in the code).
//	In fact the main structure of the program was created by modifying Paul Lewis' 
//	basiccmdline.cpp from his NCL
//
//	This code was used in Mark's dissertation, some changes were made in order to 
//	get it to compile on gcc.  It is possible that this porting introduced bugs (very
//	little debugging has been done on UNIX platforms).	I would suggest checking 
//	the simulator by generating data on trees with short branches, etc.
 

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "matrices.hpp"
#include "jph.hpp"
#include <cassert>
//From John's Mr Bayes

/*--------------------------------------------------------------------------------------------------
|
|	AllocMatrix
|
|	Allocate and set up a 'rows' x 'cols' matrix of arbitrary type.	 Storage is allocated for both
|	a vector of row pointers and the matrix itself.	 The matrix is allocated as a single block, so
|	that its elements may be referenced either as a[i][j] or as (*a)[cols*i + j].
*/

int AllocMatrix (VoidPtr pA, size_t elSize, unsigned rows, unsigned cols)

{
	unsigned			i;
	unsigned	rowBytes;
	char		*p, **a;

	if (*(VoidPtr *)pA != NULL)
		DeallocMatrix(pA);

	if ((a = (char **)calloc(rows, sizeof(VoidPtr))) == NULL)
		return ERROR;
	*(VoidPtr *)pA = a;

	rowBytes = (unsigned)(cols * elSize);

	if ((a[0] = (char *)calloc(rows, rowBytes)) == NULL)
		{
		free(a);
		return ERROR;
		}

	for (i = 0, p = *a; i < rows; i++, p += rowBytes)
		a[i] = p;
		
	return NO_ERROR;
}

/*--------------------------------------------------------------------------------------------------
|
|	DeallocMatrix
|
|	Deallocate memory for matrix allocated by AllocMatrix.
*/

void DeallocMatrix (VoidPtr pA)

{
	char		**a;
	
	if (*(VoidPtr *)pA != NULL)
		{
		a = (char **)(*(VoidPtr *)pA);
		free(a[0]);
		free(a);
		//if(a[0])	delete[]a[0];
		//delete []a;
		*(VoidPtr *)pA = NULL;
		}
}

/*--------------------------------------------------------------------------------------------------
|
|	ShortMatrix
|
|	Set row pointers into a 'rows' x 'cols' matrix whose storage is at 'buffer'.  This storage
|	can be defined as matrix[M][N] in the caller or can be any buffer of size at least rows x cols.
*/

short **ShortMatrix (unsigned rows, unsigned cols, short **a, short *buffer)

{
	unsigned			i;
	short		*p;

	for (i = 0, p = buffer; i < rows; i++, p += cols)
		a[i] = p;
		
	return a;
}

/*--------------------------------------------------------------------------------------------------
|
|	DoubleMatrix
|
|	Set row pointers into a 'rows' x 'cols' matrix whose storage is at 'buffer'.  This storage
|	can be defined as matrix[M][N] in the caller or can be any buffer of size at least rows x cols.
*/

double **DoubleMatrix (unsigned rows, unsigned cols, double **a, double *buffer)
{
	unsigned			i;
	double		*p;

	for (i = 0, p = buffer; i < rows; i++, p += cols)
		a[i] = p;
		
	return a;
}

/*--------------------------------------------------------------------------------------------------
|
|	SetIdentityMatrix
|
|	Initialize a matrix to the identity matrix.
*/

void SetIdentityMatrix (double **a, unsigned n)

{
	unsigned			i, j;

	for (i = 0; i < n; i++)
		{
		for (j = 0; j < n; j++)
			a[i][j] = 0.0;
		a[i][i] = 1.0;
		}
}


/*--------------------------------------------------------------------------------------------------
|
|	CopyDoubleMatrix
|
|	Copy matrix 'a' to 'b'.
*/

void CopyDoubleMatrix (double **a, double **b, unsigned m, unsigned n)

{
	(void)memcpy(b[0], a[0], sizeof(double) * m * n);
}


void ListDoubleMatrix (char *title, double **a, unsigned m, unsigned n)

{
	unsigned		i, j;
	
	printf("%s\n", title);
	for (i = 0; i < m; i++)
		{
		for (j = 0; j < n; j++)
			printf("%18.10g", a[i][j]);
		printf("\n");
		}
}


void ListDoubleVector (char *title, double *a, unsigned n)

{
	unsigned		i;
	
	printf("%s", title);
	for (i = 0; i < n; i++)
		printf("%18.10g", a[i]);
	printf("\n");
}





double **psdmatrix (unsigned dim)
// allocate a complex square matrix with subscript 
// range m[0..dim-1][0..dim-1]
{
	if (dim == 0)
		return NULL;
	double	**m;	
	m = new double *[dim];
	*m = new double[dim*dim];
	for (unsigned i = 1; i < dim; i++)
		*(m+i)	=*(m+i-1) + dim;
	return m;
}



void free_psdmatrix (double **m) 
{
	if (m == NULL)
		return;
	delete [] *m;
	delete [] m;	
}

double ***psdmatrices (unsigned dim,unsigned num)
// allocate a double square matrix with subscript 
// range m[0..dim-1][0..dim-1]
{	double	***m;
	
	m=new double **[num];
	*m=new double *[num*dim];
	**m=new double [num*dim*dim];
	for (unsigned i=1; i < num*dim; i++)
		m[0][i]=m[0][i-1]+dim;
	for (unsigned i=1; i < num; i++)
		m[i]=m[i-1]+dim;
	return m; // return pointer to array of pointers to rows
}

void free_psdmatrices (double ***m,unsigned num)
{	
	assert(num>0);
	if (m == NULL)
		return;
	if (*m) {
		delete []**m;
		delete [] *m;
	}
	delete [] m;
	
}

double **new_RectMats(unsigned rows,unsigned cols)
{	assert(rows>0 && cols>0);
	double **temp;
	temp=new double *[rows];
	temp[0]=new double[rows*cols];
	for (unsigned i=1; i < rows; i++)
		temp[i]=&(temp[i-1][cols]);
	return temp;
}
void free_RectMats(double **arr)
{	if (arr)	delete [] *arr;
	delete [] arr;
}

void copy_psdmatrix (double **from, double **to, unsigned dim)

{
	unsigned		row, col;
	for (row = 0; row < dim; row++)
		{for (col = 0; col < dim; col++)
			{to[row][col] = from[row][col];
			}
		}

}

void copy_psdmatrix (double *from, double *to, unsigned dim)
{	//MTH version
	assert(dim < 180);
	unsigned dimt=dim*dim;
	for (unsigned row=0; row < dimt; row++)
		*to++=*from++;
		
}


