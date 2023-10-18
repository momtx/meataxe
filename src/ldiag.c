////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Functions for drawing lattices
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"
#include <string.h>



/// @defgroup ldiag Lattice drawing
/// @{
/// @details
/// The lattice drawing functions can be used to draw a modular lattice. The 
/// algorithm calculates x and y coordinates for each node of the lattice 
/// specifying the point where the node should be drawn. Note: this algorithm
/// is far from perfect. It tries in some way to minimize the number of crossings
/// between incidence lines, but the result should always be understood as a first
/// approximation to a `beautiful' diagram. Also note that the drawing of nodes and
/// lines is left to the application.



/// @class LdNode_t
/// Lattice drawing node data
/// The LdNode_t holds all per-node data used internally by the lattice 
/// drawing algorithms. Each node has a single number (unsigned long) of 
/// user-defined data. 
/// This field may be used by the application to attach additional 
/// information to the nodes. It is not used by the drawing algorithm. 
/// PosX and PosY contain the x and y position of the node an may be
/// read by the application. All other fields are for internal use only.

/// @class LdLattice_t
/// Lattice drawing data structure
/// The |LdLattice| holds all data used internally by the lattice drawing
/// algorithms. |Node| is a list of the nodes. The nodes may appear in any 
/// order, they need not be sorted, and node 0 need not be the bottom node.

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Create a lattice drawing structure.
/// This function allocates and initializes a new LdLattice structure 
/// with a given number of nodes. The number of nodes of an existing 
/// LdLattice structure cannot be changed. When it is no longer needed, 
/// the data structure must be freed with ldFree(). 
///
/// Initially the lattice has no incidences. Before node positions are
/// calculated with ldSetPositions(), all incidences must be entered
/// using ldAddIncidence().
/// @param num_nodes Number of nodes in the lattice.
/// @return Lattice data structure or NULL on error.

LdLattice_t *ldAlloc(int num_nodes)
{
    LdLattice_t *x;

    /* Allocate memory
       --------------- */
    x = ALLOC(LdLattice_t);
    if (x == NULL)
    {
	mtxAbort(MTX_HERE,"Cannot allocate lattice structure");
	return NULL;
    }
    x->Nodes = NALLOC(LdNode_t,num_nodes);
    if (x->Nodes == NULL)
    {
	sysFree(x);
	mtxAbort(MTX_HERE,"Cannot allocate <Nodes>");
	return NULL;
    }
    x->IsSub = NALLOC(int,num_nodes * num_nodes);
    if (x->IsSub == NULL)
    {
	sysFree(x->Nodes);
	sysFree(x);
	mtxAbort(MTX_HERE,"Cannot allocate <IsSub>");
	return NULL;
    }

    /* Initialize everything
       --------------------- */
    x->NNodes = num_nodes;
    memset(x->Nodes,0,sizeof(LdNode_t) * num_nodes);
    memset(x->IsSub,0,sizeof(int) * num_nodes * num_nodes);

    return x;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Free a lattice drawing structure.
/// @param l Pointer to the data structure.
/// @return $0$ on success, $-1$ on error.
/// This function frees an |LdLattice| structure including any internally
/// allocated memory.

int ldFree(LdLattice_t *l)
{
    if (l->Nodes != NULL)
    {
        memset(l->Nodes,0,sizeof(LdNode_t) * l->NNodes);
	sysFree(l->Nodes);
    }
    if (l->IsSub != NULL)
    {
	sysFree(l->IsSub);
        memset(l->IsSub,0,sizeof(int) * l->NNodes * l->NNodes);
    }
    memset(l,0,sizeof(*l));
    sysFree(l);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add an incidence relation.
/// This function adds an incidence relation between two nodes to a given 
/// lattice. Both @a sub and @a sup must be valid node numbers, i.e. greater 
/// or equal to zero and less than the number of nodes. Apart from this 
/// range check, no further plausibility tests are performed.
/// @param lat Pointer to the lattice data structure.
/// @param sub Number of the `lower' node (contained in @a sup).
/// @param sup Number of the `upper' node (containing @a sub).
/// @return 0 on success, -1 on error.

int ldAddIncidence(LdLattice_t *lat, int sub, int sup)
{
    if (sub < 0 || sub >= lat->NNodes)
    {
	mtxAbort(MTX_HERE,"sub = %d: %s",sub,MTX_ERR_BADARG);
	return -1;
    }
    if (sup < 0 || sup >= lat->NNodes)
    {
	mtxAbort(MTX_HERE,"sup = %d: %s",sup,MTX_ERR_BADARG);
	return -1;
    }
    LD_ISSUB(lat,sub,sup) = 1;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

LdLattice_t *ldFactor(LdLattice_t *l, int min, int max)
{
    min = max = 0;
    return l;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Find the bottom node
   
static int FindBottom(LdLattice_t *l)
{
    int i, k;

    if (l->NNodes == 0)
	return -1;

    /* Search bottom, starting with node 0
       ----------------------------------- */
    i = 0;
    do 
    {
	for (k = l->NNodes - 1; k > 0 && !LD_ISSUB(l,k,i); --k);
	if (k >= 0)
	    i = k;
    } while (k > 0);
    
    return i;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Set Layer numbers 

static int FindLayers(LdLattice_t *l)
{
    int i;
    int finished;

    /* Clear layer numbers
       ------------------- */
    for (i = 0; i < l->NNodes; ++i)
	l->Nodes[i].Layer = -1;
    l->NLayers = -1;

    /* Find the bottom node and assign layer number 0
       ---------------------------------------------- */
    i = FindBottom(l);
    if (i < 0)
    {
	mtxAbort(MTX_HERE,"Cannot find bottom node");
	return -1;
    }
    l->Nodes[i].Layer = 0;

    /* Set the remaining layer numbers
       ------------------------------- */
    for (i = 0, finished = 0; !finished; ++i)
    {
	int k;
	finished = 1;
	for (k = 0; k < l->NNodes; ++k)
	{
	    int m;
	    if (l->Nodes[k].Layer != i)
		continue;
	    for (m = 0; m < l->NNodes; ++m)
	    {
		if (LD_ISSUB(l,k,m))
		{
		    if (l->Nodes[m].Layer >= 0 && l->Nodes[m].Layer != i + 1)
			mtxAbort(MTX_HERE,"Inconsistent layer numbers - "
			    "lattice is not modular!");
		    finished = 0;
		    l->Nodes[m].Layer = i + 1;
		}
	    }
	}
    }
    l->NLayers = i;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/* Set all y positions */

static int setYPositions(LdLattice_t *l)

{
    int i;
    double step, offset;

    if (l->NNodes == 0)
	return 0;
    if (l->NLayers == 1)
    {
	offset = 0.5;
	step = 0;
    }
    else
    {
	offset = 0.0;
	step = 1.0 / (l->NLayers - 1);
    }

    /* Set all Y Positions
       ------------------- */
    for (i = 0; i < l->NNodes; ++i)
	l->Nodes[i].PosY = offset + l->Nodes[i].Layer * step;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Set all x positions initially. Places nodes of one layer at equidistant x positons.

static int setInitialXPositions(LdLattice_t *l)
{
    int layer;

    for (layer = 0; layer < l->NLayers; ++layer)
    {
	int i;
	int nnodes;
	int count;
	double offset, step;

	for (nnodes = 0, i = 0; i < l->NNodes; ++i)
	{
	    if (l->Nodes[i].Layer == layer)
		++nnodes;
	}
	if (nnodes == 0)
	{
	    mtxAbort(MTX_HERE,"No nodes in layer %d - invalid lattice",layer);
	    return -1;
	}
	step = 1.0 / nnodes;
	offset = step / 2 ;
	for (count = 0, i = 0; i < l->NNodes; ++i)
	{
	    if (l->Nodes[i].Layer == layer)
		l->Nodes[i].PosX = offset + step * count++;
	}
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculate scores for x optimization

static void CalcScores(LdLattice_t* l, int direction)
{
   for (int i = 0; i < l->NNodes; ++i) {
      l->Nodes[i].Score = 0.0;
      l->Nodes[i].ScoreCount = 0;
   }

   for (int i = 0; i < l->NNodes; ++i) {
      for (int k = 0; k < l->NNodes; ++k) {
         if (!LD_ISSUB(l, i, k)) { continue; }

         if (direction > 0) {
            l->Nodes[i].Score += l->Nodes[k].PosX;
            ++l->Nodes[i].ScoreCount;
         } else {
            l->Nodes[k].Score += l->Nodes[i].PosX;
            ++l->Nodes[k].ScoreCount;
         }

         if (l->Nodes[i].Layer <= l->NLayers / 2) {
            l->Nodes[i].Score += 2 * l->Nodes[k].PosX;
            l->Nodes[i].ScoreCount += 2;
         }
         else {
            l->Nodes[k].Score += 2 * l->Nodes[i].PosX;
            l->Nodes[k].ScoreCount += 2;
         }
      }
   }

   for (int i = 0; i < l->NNodes; ++i) {
      if (l->Nodes[i].ScoreCount != 0)
         l->Nodes[i].Score /= l->Nodes[i].ScoreCount;
   }

}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Optimize X positions based on current scores
///
/// This function orders all nodes within each layer by their score, as found
/// in the Score field of the LdNode_t structure.
///
/// The return value is the number of changes made.

static int ReOrder(LdLattice_t *l)
{
    int i;
    int num_changes = 0;

    for (i = 0; i < l->NNodes; ++i)
    {
	int layer = l->Nodes[i].Layer;
	int k;
	for (k = 0; k < l->NNodes; ++k)
	{
	    double pos_k = l->Nodes[k].PosX;
	    double pos_i = l->Nodes[i].PosX;
	    double score_k = l->Nodes[k].Score;
	    double score_i = l->Nodes[i].Score;

	    if (l->Nodes[k].Layer != layer)
		continue;
	    if (((score_k < score_i && pos_k > pos_i) ||
	        (score_k > score_i && pos_k < pos_i))
		&& mtxRandomInt(100) > 50

		)
	    {
		double tmp;
		++num_changes;
/*printf("exchang %d(pos=%.2f, score=%.2f) <--> %d(pos=%.2f, score=%.2f)\n",
  i,pos_i,score_i,k,pos_k,score_k);*/
		tmp = l->Nodes[i].PosX;
		l->Nodes[i].PosX = l->Nodes[k].PosX;
		l->Nodes[k].PosX = tmp;
/*
		tmp = l->Nodes[i].Score;
		l->Nodes[i].Score = l->Nodes[k].Score;
		l->Nodes[k].Score = tmp;
*/
/*
printf("        %d(pos=%.2f, score=%.2f)      %d(pos=%.2f, score=%.2f)\n",
  i,pos_i,score_i,k,pos_k,score_k);
*/
	    }
	}
    }
    return num_changes;
}



/* --------------------------------------------------------------------------
   setXPositions() - Set all X positions

   This function calculates the x positions of all nodes. We repeat the 
   basic step - calculate scores, then reorder - until the configuration is
   stable or the maximum number of rounds is reached.
   -------------------------------------------------------------------------- */

static int setXPositions(LdLattice_t* l)
{
   if (setInitialXPositions(l) != 0) {
      return -1;
   }
   for (int round = 0; round < 50; ++round) {
      CalcScores(l, 0);
      int numA = ReOrder(l);
      //printf("Round %da: %d changes\n", round, numA);
      CalcScores(l, 1);
      int numB = ReOrder(l);
      //printf("Round %db: %d changes\n", round, numB);
      if (numA == 0 && numB == 0) { break; }
   }

   return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculate node positions.
/// This function calculates the x and y coordinates for a lattice drawing.
/// |l| is a pointer to a lattice drawing data structure. All incidences of
/// the lattice must have been defined using |ldAddIncidence()|. On successful
/// return (return value $0$) the node positions are stored in the |PosX| and
/// |PosY| fields of the |LdNode_t| structures contained in |l|.
/// Both x and y coordinates are normalized to the interval $[0,1]$.
/// @param l Pointer to the lattice data structure.
/// @return $0$ on success, $-1$ on error.

int ldSetPositions(LdLattice_t *l)
{
    if (FindLayers(l) != 0)
    {
	mtxAbort(MTX_HERE,"Cannot set layers");
	return -1;
    }
    if (setYPositions(l))
    {
	mtxAbort(MTX_HERE,"Error setting y positions");
	return -1;
    }
    if (setXPositions(l))
    {
	mtxAbort(MTX_HERE,"Error setting x positions");
	return -1;
    }
    return 0;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
