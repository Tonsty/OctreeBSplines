/*
 * C code from the article
 * "An Implicit Surface Polygonizer"
 * http::www.unchainedgeometry.com/jbloom/papers/polygonizer.pdf
 * by Jules Bloomenthal, jules@bloomenthal.com
 * in "Graphics Gems IV", Academic Press, 1994 */

/* implicit.c
 *     an implicit surface polygonizer, translated from Mesa
 *     applications should call polygonize()
 *
 * To compile a test program for ASCII output:
 *     cc implicit.c -o implicit -lm
 *
 * To compile a test program for display on an SGI workstation:
 *     cc -DSGIGFX implicit.c -o implicit -lgl_s -lm
 *
 * Authored by Jules Bloomenthal, Xerox PARC.
 * Copyright (c) Xerox Corporation, 1991.  All rights reserved.
 * Permission is granted to reproduce, use and distribute this code for
 * any and all purposes, provided that this notice appears in all copies.  */

/* A Brief Explanation
 * The code consists of a test program and the polygonizer itself.
 *
 * In the test program:
 * torus(), sphere() and blob() are primitive functions used to calculate
 *    the implicit value at a point (x,y,z)
 * triangle() is a 'callback' routine that is called by the polygonizer
 *    whenever it produces a new triangle
 * to select torus, sphere or blob, change the argument to polygonize()
 * if openGL supported, open window, establish perspective and viewport,
 *    create closed line loops, during polygonization, and slowly spin object
 * if openGL not supported, output vertices and triangles to stdout
 *
 * The main data structures in the polygonizer represent a hexahedral lattice,
 * ie, a collection of semi-adjacent cubes, represented as cube centers, corners,
 * and edges. The centers and corners are three-dimensional indices rerpesented
 * by integer i,j,k. The edges are two three-dimensional indices, represented
 * by integer i1,j1,k1,i2,j2,k2. These indices and associated data are stored
 * in hash tables.
 *
 * The client entry to the polygonizer is polygonize(), called from main().
 * This routine first allocates memory for the hash tables for the cube centers,
 * corners, and edges that define the polygonizing lattice. It then finds a start
 * point, ie, the center of the first lattice cell. It pushes this cell onto an
 * initially empty stack of cells awaiting processing. It creates the first cell
 * by computing its eight corners and assigning them an implicit value.
 *
 * polygonize() then enters a loop in which a cell is popped from the stack,
 * becoming the 'active' cell c. c is (optionally) decomposed (ie, subdivided)
 * into six tetrahedra; within each transverse tetrahedron (ie, those that
 * intersect the surface), one or two triangles are produced.
 *
 * The six faces of c are tested for intersection with the implicit surface; for
 * a transverse face, a new cube is generated and placed on the stack.

 * Some of the more important routines include:
 *
 * testface (called by polygonize): test given face for surface intersection;
 *    if transverse, create new cube by creating four new corners.
 * setcorner (called by polygonize, testface): create new cell corner at given
 *    (i,j,k), compute its implicit value, and add to corners hash table.
 * find (called by polygonize): search for point with given polarity
 * dotet (called by polygonize) set edge vertices, output triangle by
 *    invoking callback
 *
 * The section Cubical Polygonization contains routines to polygonize directly
 * from the lattice cell rather than first decompose it into tetrahedra;
 * dotet, however, is recommended over docube.
 *
 * The section Storage provides routines to handle the linked lists
 * in the hash tables.
 *
 * The section Vertices contains the following routines.
 * vertid (called by dotet): given two corner indices defining a cell edge,
 *    test whether the edge has been stored in the hash table; if so, return its
 *    associated vertex index. If not, compute intersection of edge and implicit
 *    surface, compute associated surface normal, add vertex to mesh array, and
 *    update hash tables
 * converge (called by polygonize, vertid): find surface crossing on edge */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/types.h>

#include "polygonizer.h"

#define TET	0   /* use tetrahedral decomposition */
#define NOTET 1 /* no tetrahedral decomposition  */

#define RES	10  /* # converge iterations    */

#define L	0   /* left direction:	-x, -i  */
#define R	1   /* right direction:	+x, +i  */
#define B	2   /* bottom direction: -y, -j */
#define T	3   /* top direction:	+y, +j  */
#define N	4   /* near direction:	-z, -k  */
#define F	5   /* far direction:	+z, +k  */
#define LBN	0   /* left bottom near corner  */
#define LBF	1   /* left bottom far corner   */
#define LTN	2   /* left top near corner     */
#define LTF	3   /* left top far corner      */
#define RBN	4   /* right bottom near corner */
#define RBF	5   /* right bottom far corner  */
#define RTN	6   /* right top near corner    */
#define RTF	7   /* right top far corner     */

/* the LBN corner of cube (i, j, k), corresponds with location
 * (start.x+(i-.5)*size, start.y+(j-.5)*size, start.z+(k-.5)*size) */

#define RAND()	    ((rand()&32767)/32767.)    /* random number between 0 and 1 */
#define HASHBIT	    (5)
#define HASHSIZE    (size_t)(1<<(3*HASHBIT))   /* hash table size (32768) */
#define MASK	    ((1<<HASHBIT)-1)
#define HASH(i,j,k) ((((((i)&MASK)<<HASHBIT)|((j)&MASK))<<HASHBIT)|((k)&MASK))
#define BIT(i, bit) (((i)>>(bit))&1)
#define FLIP(i,bit) ((i)^1<<(bit)) /* flip the given bit of i */

char *mycalloc (int nitems, int nbytes);
void makecubetable ();
TEST find (int sign, PROCESS *p, double x, double y, double z);
void converge (POINT *p1, POINT *p2, double v, double (*function)(double,double,double), POINT *p);
CORNER *setcorner (PROCESS *p, int i, int j, int k);
int setcenter(CENTERLIST *table[], int i, int j, int k);
int dotet (CUBE *cube, int c1, int c2, int c3, int c4, PROCESS *p);
int docube (CUBE *cube, PROCESS *p);
void testface (int i, int j, int k, CUBE *old, int face, int c1, int c2, int c3, int c4, PROCESS *p);
int vertid (CORNER *c1, CORNER *c2, PROCESS *p);
void vnormal (POINT *point, PROCESS *p, POINT *v);
void addtovertices (VERTICES *vertices, VERTEX v);
int addtotriangles (TRIANGLES *triangles, TRIANGLE t);
void converge2 (POINT *p1, POINT *p2, double v1, double v2, POINT *p);

/**** A Test Program ****/


/* torus: a torus with major, minor radii = 0.5, 0.1, try size = .05 */

double torus (double x, double y, double z)
{
    double x2 = x*x, y2 = y*y, z2 = z*z;
    double a = x2+y2+z2+(0.5*0.5)-(0.1*0.1);
    return a*a-4.0*(0.5*0.5)*(y2+z2);
}


/* sphere: an inverse square function (always positive) */

double sphere (double x, double y, double z)
{
    double rsq = x*x+y*y+z*z;
    return 1.0/(rsq < 0.00001? 0.00001 : rsq);
}


/* blob: a three-pole blend function, try size = .1 */

double blob (double x, double y, double z)
{
    return 4.0-sphere(x+1.0,y,z)-sphere(x,y+1.0,z)-sphere(x,y,z+1.0);
}

/***********************************************************************/

int gntris;	     /* global needed by application */
VERTICES gvertices;  /* global needed by application */
TRIANGLES gtriangles; /* global needed by application */


/* triangle: called by polygonize() for each triangle; write to stdout */

int triangle (int i1, int i2, int i3, VERTICES vertices)
{
    gvertices = vertices;
    gntris++;
	fprintf(stdout, "%d %d %d\n", i1, i2, i3);
    return 1;
}

/* triangle2: called by polygonize() for each triangle; write to triangles buffer */

int triangle2 (int i1, int i2, int i3, VERTICES vertices)
{
	TRIANGLE t = {i1, i2, i3};
	gvertices = vertices;
	gntris++;
	addtotriangles(&gtriangles, t);
	return 1;
}

/**********************************************************************/


/**** An Implicit Surface Polygonizer ****/


/* polygonize: polygonize the implicit surface function
 *   arguments are:
 *	 double function (x, y, z)
 *		 double x, y, z (an arbitrary 3D point)
 *	     the implicit surface function
 *	     return negative for inside, positive for outside
 *	 double size
 *	     width of the partitioning cube
 *	 int bounds
 *	     max. range of cubes (+/- on the three axes) from first cube
 *	 double x, y, z
 *	     coordinates of a starting point on or near the surface
 *	     may be defaulted to 0., 0., 0.
 *	 int triproc (i1, i2, i3, vertices)
 *		 int i1, i2, i3 (indices into the vertex array)
 *		 VERTICES vertices (the vertex array, indexed from 0)
 *	     called for each triangle
 *	     the triangle coordinates are (for i = i1, i2, i3):
 *		 vertices.ptr[i].position.x, .y, and .z
 *	     vertices are ccw when viewed from the out (positive) side
 *		 in a left-handed coordinate system
 *	     vertex normals point outwards
 *	     return 1 to continue, 0 to abort
 *	 int mode
 *	     TET: decompose cube and polygonize six tetrahedra
 *	     NOTET: polygonize cube directly
 *   returns error or NULL
 */

char *polygonize (double(*function)(double,double,double), double size, int bounds, double x, double y, double z, int (*triproc)(int,int,int,VERTICES), int mode)
{
    PROCESS p;
    int n, noabort;
    CORNER *setcorner();
    TEST in, out;

    p.function = function;
    p.triproc = triproc;
    p.size = size;
    p.bounds = bounds;
    p.delta = size/(double)(RES*RES);

    /* allocate hash tables and build cube polygon table: */
    p.centers = (CENTERLIST **) mycalloc(HASHSIZE,sizeof(CENTERLIST *));
    p.corners = (CORNERLIST **) mycalloc(HASHSIZE,sizeof(CORNERLIST *));
    p.edges =	(EDGELIST   **) mycalloc(2*HASHSIZE,sizeof(EDGELIST *));
    makecubetable();

    /* find point on surface, beginning search at (x, y, z): */
    srand(1);
    in = find(1, &p, x, y, z);
    out = find(0, &p, x, y, z);
    if (!in.ok || !out.ok) return "can't find starting point";
    converge(&in.p, &out.p, in.value, p.function, &p.start);

    /* push initial cube on stack: */
    p.cubes = (CUBES *) mycalloc(1, sizeof(CUBES)); /* list of 1 */
    p.cubes->cube.i = p.cubes->cube.j = p.cubes->cube.k = 0;
    p.cubes->next = NULL;

    /* set corners of initial cube: */
    for (n = 0; n < 8; n++)
	p.cubes->cube.corners[n] = setcorner(&p, BIT(n,2), BIT(n,1), BIT(n,0));

    p.vertices.count = p.vertices.max = 0; /* no vertices yet */
    p.vertices.ptr = NULL;

    setcenter(p.centers, 0, 0, 0);

    while (p.cubes != NULL) { /* process active cubes till none left */
	CUBE c;
	CUBES *temp = p.cubes;
	c = p.cubes->cube;

	noabort = mode == TET?
	       /* either decompose into tetrahedra and polygonize: */
	       dotet(&c, LBN, LTN, RBN, LBF, &p) &&
	       dotet(&c, RTN, LTN, LBF, RBN, &p) &&
	       dotet(&c, RTN, LTN, LTF, LBF, &p) &&
	       dotet(&c, RTN, RBN, LBF, RBF, &p) &&
	       dotet(&c, RTN, LBF, LTF, RBF, &p) &&
	       dotet(&c, RTN, LTF, RTF, RBF, &p)
	       :
	       /* or polygonize the cube directly: */
	       docube(&c, &p);
	if (! noabort) return "aborted";

	/* pop current cube from stack */
	p.cubes = p.cubes->next;
	free((char *) temp);
	/* test six face directions, maybe add to stack: */
	testface(c.i-1, c.j, c.k, &c, L, LBN, LBF, LTN, LTF, &p);
	testface(c.i+1, c.j, c.k, &c, R, RBN, RBF, RTN, RTF, &p);
	testface(c.i, c.j-1, c.k, &c, B, LBN, LBF, RBN, RBF, &p);
	testface(c.i, c.j+1, c.k, &c, T, LTN, LTF, RTN, RTF, &p);
	testface(c.i, c.j, c.k-1, &c, N, LBN, LTN, RBN, RTN, &p);
	testface(c.i, c.j, c.k+1, &c, F, LBF, LTF, RBF, RTF, &p);
    }
    return NULL;
}


/* testface: given cube at lattice (i, j, k), and four corners of face,
 * if surface crosses face, compute other four corners of adjacent cube
 * and add new cube to cube stack */

void testface (int i, int j, int k, CUBE *old, int face, int c1, int c2, int c3, int c4, PROCESS *p)
{
    CUBE newcube;
    CUBES *oldcubes = p->cubes;
    CORNER *setcorner();
    static int facebit[6] = {2, 2, 1, 1, 0, 0};
    int n, pos = old->corners[c1]->value > 0.0 ? 1 : 0, bit = facebit[face];

    /* test if no surface crossing, cube out of bounds, or already visited: */
    if ((old->corners[c2]->value > 0) == pos &&
	(old->corners[c3]->value > 0) == pos &&
	(old->corners[c4]->value > 0) == pos) return;
    if (abs(i) > p->bounds || abs(j) > p->bounds || abs(k) > p->bounds) return;
    if (setcenter(p->centers, i, j, k)) return;

    /* create new cube: */
    newcube.i = i;
    newcube.j = j;
    newcube.k = k;
    for (n = 0; n < 8; n++) newcube.corners[n] = NULL;
    newcube.corners[FLIP(c1, bit)] = old->corners[c1];
    newcube.corners[FLIP(c2, bit)] = old->corners[c2];
    newcube.corners[FLIP(c3, bit)] = old->corners[c3];
    newcube.corners[FLIP(c4, bit)] = old->corners[c4];
    for (n = 0; n < 8; n++)
	if (newcube.corners[n] == NULL)
	    newcube.corners[n] = setcorner(p, i+BIT(n,2), j+BIT(n,1), k+BIT(n,0));

    /*add cube to top of stack: */
    p->cubes = (CUBES *) mycalloc(1, sizeof(CUBES));
    p->cubes->cube = newcube;
    p->cubes->next = oldcubes;
}


/* setcorner: return corner with the given lattice location
   set (and cache) its function value */

CORNER *setcorner (PROCESS *p, int i, int j, int k)
{
    /* for speed, do corner value caching here */
    CORNER *c = (CORNER *) mycalloc(1, sizeof(CORNER));
    int index = HASH(i, j, k);
    CORNERLIST *l = p->corners[index];
    c->i = i; c->x = p->start.x+((double)i-.5)*p->size;
    c->j = j; c->y = p->start.y+((double)j-.5)*p->size;
    c->k = k; c->z = p->start.z+((double)k-.5)*p->size;
    for (; l != NULL; l = l->next)
	if (l->i == i && l->j == j && l->k == k) {
	    c->value = l->value;
	    return c;
	    }
    l = (CORNERLIST *) mycalloc(1, sizeof(CORNERLIST));
    l->i = i; l->j = j; l->k = k;
    l->value = c->value = p->function(c->x, c->y, c->z);
    l->next = p->corners[index];
    p->corners[index] = l;
    return c;
}


/* find: search for point with value of given sign (0: neg, 1: pos) */

TEST find (int sign, PROCESS *p, double x, double y, double z)
{
    int i;
    TEST test;
    double range = p->size;
    test.ok = 1;
    for (i = 0; i < 10000; i++) {
	test.p.x = x+range*(RAND()-0.5);
	test.p.y = y+range*(RAND()-0.5);
	test.p.z = z+range*(RAND()-0.5);
	test.value = p->function(test.p.x, test.p.y, test.p.z);
	if (sign == (test.value > 0.0)) return test;
	range = range*1.0005; /* slowly expand search outwards */
    }
    test.ok = 0;
    return test;
}


/**** Tetrahedral Polygonization ****/


/* dotet: triangulate the tetrahedron
 * b, c, d should appear clockwise when viewed from a
 * return 0 if client aborts, 1 otherwise */

int dotet (CUBE *cube, int c1, int c2, int c3, int c4, PROCESS *p)
{
    CORNER *a = cube->corners[c1];
    CORNER *b = cube->corners[c2];
    CORNER *c = cube->corners[c3];
    CORNER *d = cube->corners[c4];
    int index = 0, apos, bpos, cpos, dpos, e1, e2, e3, e4, e5, e6;
    if (apos = (a->value > 0.0)) index += 8;
    if (bpos = (b->value > 0.0)) index += 4;
    if (cpos = (c->value > 0.0)) index += 2;
    if (dpos = (d->value > 0.0)) index += 1;
    /* index is now 4-bit number representing one of the 16 possible cases */
    if (apos != bpos) e1 = vertid(a, b, p);
    if (apos != cpos) e2 = vertid(a, c, p);
    if (apos != dpos) e3 = vertid(a, d, p);
    if (bpos != cpos) e4 = vertid(b, c, p);
    if (bpos != dpos) e5 = vertid(b, d, p);
    if (cpos != dpos) e6 = vertid(c, d, p);
    /* 14 productive tetrahedral cases (0000 and 1111 do not yield polygons */
	switch (index) {
	case 1:	 return p->triproc(e5, e6, e3, p->vertices);
	case 2:	 return p->triproc(e2, e6, e4, p->vertices);
	case 3:	 return p->triproc(e3, e5, e4, p->vertices) &&
					p->triproc(e3, e4, e2, p->vertices);
	case 4:	 return p->triproc(e1, e4, e5, p->vertices);
	case 5:	 return p->triproc(e3, e1, e4, p->vertices) &&
					p->triproc(e3, e4, e6, p->vertices);
	case 6:	 return p->triproc(e1, e2, e6, p->vertices) &&
					p->triproc(e1, e6, e5, p->vertices);
	case 7:	 return p->triproc(e1, e2, e3, p->vertices);
	case 8:	 return p->triproc(e1, e3, e2, p->vertices);
	case 9:	 return p->triproc(e1, e5, e6, p->vertices) &&
					p->triproc(e1, e6, e2, p->vertices);
	case 10: return p->triproc(e1, e3, e6, p->vertices) &&
					p->triproc(e1, e6, e4, p->vertices);
	case 11: return p->triproc(e1, e5, e4, p->vertices);
	case 12: return p->triproc(e3, e2, e4, p->vertices) &&
					p->triproc(e3, e4, e5, p->vertices);
	case 13: return p->triproc(e6, e2, e4, p->vertices);
	case 14: return p->triproc(e5, e3, e6, p->vertices);
	}
    return 1;
}


/**** Cubical Polygonization (optional) ****/


#define LB	0  /* left bottom edge	*/
#define LT	1  /* left top edge	*/
#define LN	2  /* left near edge	*/
#define LF	3  /* left far edge	*/
#define RB	4  /* right bottom edge */
#define RT	5  /* right top edge	*/
#define RN	6  /* right near edge	*/
#define RF	7  /* right far edge	*/
#define BN	8  /* bottom near edge	*/
#define BF	9  /* bottom far edge	*/
#define TN	10 /* top near edge	*/
#define TF	11 /* top far edge	*/

static INTLISTS *cubetable[256];

/*			edge: LB, LT, LN, LF, RB, RT, RN, RF, BN, BF, TN, TF */
static int corner1[12]	   = {LBN,LTN,LBN,LBF,RBN,RTN,RBN,RBF,LBN,LBF,LTN,LTF};
static int corner2[12]	   = {LBF,LTF,LTN,LTF,RBF,RTF,RTN,RTF,RBN,RBF,RTN,RTF};
static int leftface[12]	   = {B,  L,  L,  F,  R,  T,  N,  R,  N,  B,  T,  F};
			     /* face on left when going corner1 to corner2 */
static int rightface[12]   = {L,  T,  N,  L,  B,  R,  R,  F,  B,  F,  N,  T};
			     /* face on right when going corner1 to corner2 */


/* docube: triangulate the cube directly, without decomposition */

int docube (CUBE *cube, PROCESS *p)
{
    INTLISTS *polys;
    int i, index = 0;
    for (i = 0; i < 8; i++) if (cube->corners[i]->value > 0.0) index += (1<<i);
    for (polys = cubetable[index]; polys; polys = polys->next) {
	INTLIST *edges;
	int a = -1, b = -1, count = 0;
	for (edges = polys->list; edges; edges = edges->next) {
	    CORNER *c1 = cube->corners[corner1[edges->i]];
	    CORNER *c2 = cube->corners[corner2[edges->i]];
	    int c = vertid(c1, c2, p);
	    if (++count > 2 && ! p->triproc(a, b, c, p->vertices)) return 0;
	    if (count < 3) a = b;
	    b = c;
	}
    }
    return 1;
}


/* nextcwedge: return next clockwise edge from given edge around given face */

int nextcwedge (int edge, int face)
{
    switch (edge) {
	case LB: return (face == L)? LF : BN;
	case LT: return (face == L)? LN : TF;
	case LN: return (face == L)? LB : TN;
	case LF: return (face == L)? LT : BF;
	case RB: return (face == R)? RN : BF;
	case RT: return (face == R)? RF : TN;
	case RN: return (face == R)? RT : BN;
	case RF: return (face == R)? RB : TF;
	case BN: return (face == B)? RB : LN;
	case BF: return (face == B)? LB : RF;
	case TN: return (face == T)? LT : RN;
	case TF: return (face == T)? RT : LF;
    }
	return -1;
}


/* otherface: return face adjoining edge that is not the given face */

int otherface (int edge, int face)
{
    int other = leftface[edge];
    return face == other? rightface[edge] : other;
}


/* makecubetable: create the 256 entry table for cubical polygonization */

void makecubetable ()
{
    int i, e, c, done[12], pos[8];
    for (i = 0; i < 256; i++) {
	for (e = 0; e < 12; e++) done[e] = 0;
	for (c = 0; c < 8; c++) pos[c] = BIT(i, c);
	for (e = 0; e < 12; e++)
	    if (!done[e] && (pos[corner1[e]] != pos[corner2[e]])) {
		INTLIST *ints = 0;
		INTLISTS *lists = (INTLISTS *) mycalloc(1, sizeof(INTLISTS));
		int start = e, edge = e;
		/* get face that is to right of edge from pos to neg corner: */
		int face = pos[corner1[e]]? rightface[e] : leftface[e];
		while (1) {
		    edge = nextcwedge(edge, face);
		    done[edge] = 1;
		    if (pos[corner1[edge]] != pos[corner2[edge]]) {
			INTLIST *tmp = ints;
			ints = (INTLIST *) mycalloc(1, sizeof(INTLIST));
			ints->i = edge;
			ints->next = tmp; /* add edge to head of list */
			if (edge == start) break;
			face = otherface(edge, face);
		    }
		}
		lists->list = ints; /* add ints to head of table entry */
		lists->next = cubetable[i];
		cubetable[i] = lists;
	    }
    }
}


/**** Storage ****/


/* mycalloc: return successful calloc or exit program */

char *mycalloc (int nitems, int nbytes)
{
   char *ptr = (char *) calloc(nitems, nbytes);
   if (ptr != NULL) return ptr;
   fprintf(stderr, "can't calloc %d bytes\n", nitems*nbytes);
   exit(1);
}


/* setcenter: set (i,j,k) entry of table[]
 * return 1 if already set; otherwise, set and return 0 */

int setcenter(CENTERLIST *table[], int i, int j, int k)
{
    int index = HASH(i, j, k);
    CENTERLIST *newlist, *l, *q = table[index];
    for (l = q; l != NULL; l = l->next)
	if (l->i == i && l->j == j && l->k == k) return 1;
    newlist = (CENTERLIST *) mycalloc(1, sizeof(CENTERLIST));
    newlist->i = i; newlist->j = j; newlist->k = k; newlist->next = q;
    table[index] = newlist;
    return 0;
}


/* setedge: set vertex id for edge */

void setedge (EDGELIST *table[], int i1, int j1, int k1, int i2, int j2, int k2, int vid)
{
    unsigned int index;
    EDGELIST *newlist;
    if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
	int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
    }
    index = HASH(i1, j1, k1) + HASH(i2, j2, k2);
    newlist = (EDGELIST *) mycalloc(1, sizeof(EDGELIST));
    newlist->i1 = i1; newlist->j1 = j1; newlist->k1 = k1;
    newlist->i2 = i2; newlist->j2 = j2; newlist->k2 = k2;
    newlist->vid = vid;
    newlist->next = table[index];
    table[index] = newlist;
}


/* getedge: return vertex id for edge; return -1 if not set */

int getedge (EDGELIST *table[], int i1, int j1, int k1, int i2, int j2, int k2)
{
    EDGELIST *q;
    if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
	int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
    };
    q = table[HASH(i1, j1, k1)+HASH(i2, j2, k2)];
    for (; q != NULL; q = q->next)
	if (q->i1 == i1 && q->j1 == j1 && q->k1 == k1 &&
	    q->i2 == i2 && q->j2 == j2 && q->k2 == k2)
	    return q->vid;
    return -1;
}


/**** Vertices ****/


/* vertid: return index for vertex on edge:
 * c1->value and c2->value are presumed of different sign
 * return saved index if any; else compute vertex and save */

int vertid (CORNER *c1, CORNER *c2, PROCESS *p)
{
    VERTEX v;
    POINT a, b;
    int vid = getedge(p->edges, c1->i, c1->j, c1->k, c2->i, c2->j, c2->k);
    if (vid != -1) return vid;			     /* previously computed */
    a.x = c1->x; a.y = c1->y; a.z = c1->z;
    b.x = c2->x; b.y = c2->y; b.z = c2->z;
    //converge(&a, &b, c1->value, p->function, &v.position); /* position */
    converge2(&a, &b, c1->value, c2->value, &v.position); /* position */
    vnormal(&v.position, p, &v.normal);			   /* normal */
    addtovertices(&p->vertices, v);			   /* save vertex */
    vid = p->vertices.count-1;
    setedge(p->edges, c1->i, c1->j, c1->k, c2->i, c2->j, c2->k, vid);
    return vid;
}


/* addtovertices: add v to sequence of vertices */

void addtovertices (VERTICES *vertices, VERTEX v)
{
    if (vertices->count == vertices->max) {
	int i;
	VERTEX *newvertex;
	vertices->max = vertices->count == 0 ? 10 : 2*vertices->count;
	newvertex = (VERTEX *) mycalloc(vertices->max, sizeof(VERTEX));
	for (i = 0; i < vertices->count; i++) newvertex[i] = vertices->ptr[i];
	if (vertices->ptr != NULL) free((char *) vertices->ptr);
	vertices->ptr = newvertex;
    }
    vertices->ptr[vertices->count++] = v;
}


/* addtotriangles: add t to sequence of triangles */

int addtotriangles (TRIANGLES *triangles, TRIANGLE t)
{
	if (triangles->count == triangles->max) {
		int i;
		TRIANGLE *newtriangle;
		triangles->max = triangles->count == 0 ? 10 : 2*triangles->count;
		newtriangle = (TRIANGLE *) mycalloc(triangles->max, sizeof(TRIANGLE));
		for (i = 0; i < triangles->count; i++) newtriangle[i] = triangles->ptr[i];
		if (triangles->ptr != NULL) free((char *) triangles->ptr);
		triangles->ptr = newtriangle;
	}
	triangles->ptr[triangles->count++] = t;
	return 1;
}


/* vnormal: compute unit length surface normal at point */

void vnormal (POINT *point, PROCESS *p, POINT *v)
{
    double f = p->function(point->x, point->y, point->z);
    v->x = p->function(point->x+p->delta, point->y, point->z)-f;
    v->y = p->function(point->x, point->y+p->delta, point->z)-f;
    v->z = p->function(point->x, point->y, point->z+p->delta)-f;
    f = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
    if (f != 0.0) {v->x /= f; v->y /= f; v->z /= f;}
}


/* converge: from two points of differing sign, converge to zero crossing */

void converge (POINT *p1, POINT *p2, double v, double (*function)(double,double,double), POINT *p)
{
    int i = 0;
    POINT pos, neg;
    if (v < 0) {
	pos.x = p2->x; pos.y = p2->y; pos.z = p2->z;
	neg.x = p1->x; neg.y = p1->y; neg.z = p1->z;
    }
    else {
	pos.x = p1->x; pos.y = p1->y; pos.z = p1->z;
	neg.x = p2->x; neg.y = p2->y; neg.z = p2->z;
    }
    while (1) {
	p->x = 0.5*(pos.x + neg.x);
	p->y = 0.5*(pos.y + neg.y);
	p->z = 0.5*(pos.z + neg.z);
	if (i++ == RES) return;
	if ((function(p->x, p->y, p->z)) > 0.0)
	     {pos.x = p->x; pos.y = p->y; pos.z = p->z;}
	else {neg.x = p->x; neg.y = p->y; neg.z = p->z;}
    }
}

void converge2 (POINT *p1, POINT *p2, double v1, double v2, POINT *p)
{
	double t = v1/(v1-v2);
	p->x=p1->x*(1-t)+p2->x*t;
	p->y=p1->y*(1-t)+p2->y*t;
	p->z=p1->z*(1-t)+p2->z*t;
}
