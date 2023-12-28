////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Functions for reading and writing the .cfinfo file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <string.h>
#include <stdlib.h>
#include <ctype.h>

/// @defgroup cfinfo Module Structure
/// @{

/// @class CfInfo
/// Constituent data structure.
/// The CFInfo data structure contains all information about one 
/// irreducible constituent of a module. Most of this information is
/// useful only in the context of a LatInfo_t structure.
///
/// @c dim is the dimension of the constituent. If there are more than one non-isomorphic
/// constituents of the same dimension they are numbered in order of discovery with
/// @c num = 1, 2, 3, ...
///
/// @c mult is the multiplicity of this in the module.
///
/// @c idWord and @c idPol specify the identifying word for this module.
/// The actual algebra element defined by these two inputs is a = p(w), where w is the
/// algebra element with number @c idWord produced by the @ref wgen "word generator",
/// and p(X) is the polynomial @c idPol.
///
/// Analogously, @c peakWord and @c peakPol specify the peak word for this module.
/// See the @ref prog_pwkond "pwkond program" description for details.

/// @class LatInfo_t
/// Module Structure Information.
/// This data structure contains all information about a module which is
/// gathered during the submodule lattice calculation.

////////////////////////////////////////////////////////////////////////////////////////////////////

static void WriteWord(StfData *f, uint32_t w, Poly_t *p)
{
    int i;

    stfPut(f,"[");
    stfPutU32(f,w);
    stfPut(f,",");
    stfPutInt(f,ffOrder);
    if (p == NULL)
    {
	stfPut(f,",-1");
    }
    else
    {
	stfPut(f,",");
	stfPutInt(f,p->degree);
	for (i = 0; i <= p->degree; ++i)
	{
	    stfPut(f,",");
	    stfPutInt(f,ffToInt(p->data[i]));
	}
    }
    stfPut(f,"]");
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a word number and the associated polynomial.
/// @return 0 on success, -1 on error

static void ReadWord(StfData *f, uint32_t *w, Poly_t **p, const char *fn)
{
    int fl, deg;  
    int i;

    if (stfMatch(f," [")) {
	mtxAbort(MTX_HERE,"%s: missing '['",fn);
    }
    stfGetInt(f,&i);
    *w = i;
    if (stfMatch(f,",")) {
	mtxAbort(MTX_HERE,"%s: missing ','",fn);
    }
    stfGetInt(f,&fl);
    if (stfMatch(f,",")) {
	mtxAbort(MTX_HERE,"%s: missing ','",fn);
    }
    stfGetInt(f,&deg);
    if (deg == -1) {
	*p = NULL;
    } else {
        *p = polAlloc(fl,deg);
    	for (i = 0; i <= deg; ++i)
    	{
	    int coeff;
            if (stfMatch(f,",")) 
	    {
		mtxAbort(MTX_HERE,"%s: missing ','",fn);
	    }
	    stfGetInt(f,&coeff);
	    (*p)->data[i] = ffFromInt(coeff);
    	}
    }
    if (stfMatch(f,"]")) {
	mtxAbort(MTX_HERE,"%s: missing ']'",fn);
    }
}


#define RDVEC(name,fld)\
	else if (!strcmp(c,name))\
	{\
	    if (stfMatch(f," ["))\
	    { mtxAbort(MTX_HERE,"Error in cfinfo file: Missing '['"); }\
	    for (int i = 0; i < li->nCf; ++i) \
	    {\
		int val = 0;\
		if (i > 0)\
		    stfMatch(f,",");\
		stfGetInt(f,&val);\
		li->Cf[i].fld = val;\
	    }\
	    if (stfMatch(f,"]"))\
	    { mtxAbort(MTX_HERE,"Error in cfinfo file: Missing ']'"); }\
	}


static void readCfFile(StfData* f, const char* fn, LatInfo_t* li)
{
    // Read header
    if (stfReadLine(f) || strcmp(stfGetName(f),"CFInfo"))
	mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);

    // Read data
    while (stfReadLine(f) == 0)
    {
	const char *c = stfGetName(f);
	if (!strcmp(c,"CFInfo.NCF"))
	    stfGetInt(f,&li->nCf);
	else if (!strcmp(c,"CFInfo.ConstituentNames"))
	    ;	/* Ignore */
	else if (!strcmp(c,"CFInfo.Field"))
	{   
	    stfGetInt(f,&li->field);
	    ffSetField(li->field);
	}
	else if (!strcmp(c,"CFInfo.NGen"))
	    stfGetInt(f,&li->NGen);
	RDVEC("CFInfo.Dimension",dim)
	RDVEC("CFInfo.Number",num)
	RDVEC("CFInfo.Multiplicity",mult)
	RDVEC("CFInfo.SplittingField",spl)
	RDVEC("CFInfo.NMountains",nmount)
	RDVEC("CFInfo.NDottedLines",ndotl)
	else if (!strcmp(c,"CFInfo.IdWord"))
	{
	    if (stfMatch(f," [") != 0) 
		mtxAbort(MTX_HERE,"%s: Missing '['",fn);
	    for (int i = 0; i < li->nCf; ++i)
	    {
		ReadWord(f,&(li->Cf[i].idWord),&(li->Cf[i].idPol),fn);
		if (stfMatch(f,i < li->nCf - 1 ? "," : "];") != 0)
		    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
	    }
	}
	else if (!strcmp(c,"CFInfo.PeakWord"))
	{
	    if (stfMatch(f," [") != 0) 
		mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
	    for (int i = 0; i < li->nCf; ++i)
	    {
		ReadWord(f,&(li->Cf[i].peakWord),&(li->Cf[i].peakPol),fn);
		if (stfMatch(f,i < li->nCf - 1 ? "," : "];") != 0)
		    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
	    }
	}
	else if (!strcmp(c,"CFInfo.LoewyLength"))   /* for compatibility */
	    ;	
	else if (!strcmp(c,"CFInfo.NSocles"))
	    ;	/* Is set when reading "CFInfo.Socles" */
	else if (!strcmp(c,"CFInfo.Socles"))
	{
	    if (stfMatch(f," [") != 0) 
		mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
	    for (int i = 0; stfMatch(f,"];"); ++i)
	    {
		int mult[LAT_MAXCF];
		int count = LAT_MAXCF;
		if (i > 0) stfMatch(f,",");
		stfGetVector(f,&count,mult);
		if (count != li->nCf)
		    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		latAddSocle(li,mult);
	    }
	}
	else if (!strcmp(c,"CFInfo.NHeads"))
	    ;	/* Is set when reading "CFInfo.Heads" */
	else if (!strcmp(c,"CFInfo.Heads"))
	{
	    if (stfMatch(f," [") != 0) 
		mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
	    for (int i = 0; stfMatch(f,"];"); ++i)
	    {
		int mult[LAT_MAXCF];
		int count = LAT_MAXCF;
		if (i > 0) stfMatch(f,",");
		stfGetVector(f,&count,mult);
		if (count != li->nCf)
		    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		latAddHead(li,mult);
	    }
	}
	else
	    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

LatInfo_t* latCreate(const char* baseName)
{
   MTX_ASSERT(baseName != NULL);
   MTX_ASSERT(strlen(baseName) < LAT_MAXBASENAME - 1); // TODO: remove limit
   LatInfo_t* li = mmAlloc(MTX_TYPE_LATINFO, sizeof(LatInfo_t));
   strcpy(li->BaseName, baseName);
   li->NGen = 2;
   return li;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Destroys a LatInfo_t object and releases associated resources.

void latDestroy(LatInfo_t* li)
{
   for (int i = 0; i < li->nCf; ++i) {
      CfInfo* cf = li->Cf + i;
      if (cf->idPol != NULL) {
         polFree(cf->idPol);
         cf->idPol = NULL;
      }
      if (cf->peakPol != NULL) {
         polFree(cf->peakPol);
         cf->peakPol = NULL;
      }
   }
   mmFree(li, MTX_TYPE_LATINFO);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Reads a lattice information file.
///
/// This function reads a .cfinfo file and returns a LatInfo_t object containing the data read
/// fro the file.
///
/// @param basename Base name (without the trailing ".cfinfo").

LatInfo_t* latLoad(const char* baseName)
{
   LatInfo_t* li = latCreate(baseName);

   char fn[LAT_MAXBASENAME + 20];
   sprintf(fn, "%s.cfinfo", baseName);
   StfData* f = stfOpen(fn, "r");
   readCfFile(f, fn, li);
   stfClose(f);
   return li;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Writes a lattice information file.
///
/// This function writes the contents of a LatInfo_t structure into a file.
/// The file name is constructed from the @c BaseName field of the structure
/// by appending ".cfinfo".

void latSave(const LatInfo_t *li)
{
    StfData *f;
    int i;
    int tmp[LAT_MAXCF];
    char fn[LAT_MAXBASENAME + 20];

    MTX_ASSERT(li != NULL);

    // Open the file
    snprintf(fn, sizeof(fn), "%s.cfinfo", li->BaseName);
    f = stfOpen(fn,"w");

    /* Write data
       ---------- */
    stfWriteValue(f,"CFInfo","rec()");
    stfWriteInt(f,"CFInfo.NGen",li->NGen);
    stfWriteInt(f,"CFInfo.Field",li->field);
    stfWriteInt(f,"CFInfo.NCF",li->nCf);


    stfBeginEntry(f,"CFInfo.ConstituentNames");
    stfPut(f,"[");
    for (i = 0; i < li->nCf; ++i)
    {
	stfPutString(f,latCfName(li,i));
	if (i < li->nCf-1) stfPut(f,",");
    }
    stfPut(f,"]");
    stfEndEntry(f);

#define WRVEC(name,field)\
    for (i = 0; i < li->nCf; ++i)\
        tmp[i] = li->Cf[i].field;\
    stfWriteVector(f,"CFInfo." #name,li->nCf,tmp);

    WRVEC(Dimension,dim);
    WRVEC(Number,num);
    WRVEC(Multiplicity,mult);
    WRVEC(SplittingField,spl);
    WRVEC(NMountains,nmount);
    WRVEC(NDottedLines,ndotl);

    stfBeginEntry(f,"CFInfo.PeakWord");
    stfPut(f,"[");
    for (i = 0; i < li->nCf; ++i)
    {
        WriteWord(f,li->Cf[i].peakWord,li->Cf[i].peakPol);
	if (i < li->nCf-1) stfPut(f,",");
    }
    stfPut(f,"]");
    stfEndEntry(f);

    stfBeginEntry(f,"CFInfo.IdWord");
    stfPut(f,"[");
    for (i = 0; i < li->nCf; ++i)
    {
        WriteWord(f,li->Cf[i].idWord,li->Cf[i].idPol);
	if (i < li->nCf-1) stfPut(f,",");
    }
    stfPut(f,"]");
    stfEndEntry(f);

    stfWriteInt(f,"CFInfo.NSocles",li->NSocles);
    stfBeginEntry(f,"CFInfo.Socles");
    stfPut(f,"[");
    for (i = 0; i < li->NSocles; ++i)
    {
	if (i > 0) stfPut(f,",");
	stfPutVector(f,li->nCf,li->Socle + i * li->nCf);
    }
    stfPut(f,"]");
    stfEndEntry(f);

    stfWriteInt(f,"CFInfo.NHeads",li->NHeads);
    stfBeginEntry(f,"CFInfo.Heads");
    stfPut(f,"[");
    for (i = 0; i < li->NHeads; ++i)
    {
	if (i > 0) stfPut(f,",");
	stfPutVector(f,li->nCf,li->Head + i * li->nCf);
    }
    stfPut(f,"]");
    stfEndEntry(f);

    stfClose(f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Make Constituent Name.
///
/// This function returns the name of the @a cf-th constituent of a module. The constituent name
/// consists of the dimension and an appendix which is based on @c num field in the constituent's
/// data structure. Usually the appendix is a single letter ('a', 'b', ...). If there are more than 
/// 26 constituents with the same dimension, a two-letter appendix ('aa', 'ab', etc.) is used. 
///
/// Note: The returned pointer points to a member of the @c CfInfo structure associated with the
/// constituent.

const char* latCfName(const LatInfo_t* li, int cf)
{
   MTX_ASSERT(li != NULL);
   MTX_ASSERT(cf >= 0 && cf < li->nCf);

   const CfInfo* const ip = li->Cf + cf;
   char* const n = (char*) ip->name;
   const size_t ns = sizeof(ip->name);
   unsigned long const dim = (unsigned long) ip->dim;
   unsigned long const num = (unsigned long) ip->num;

   if (num < 26U) {
      snprintf(n, ns, "%ld%c", dim, (char)num + 'a');
   }
   else if (num < 26U * 26) {
      snprintf(n, ns, "%ld%c%c", dim, (char)(num / 26 - 1) + 'a', (char)(num % 26) + 'a');
   }
   else {
      snprintf(n, ns, "%ldcf%ld", dim, num);
   }
   return ip->name;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Add a Layer to the Socle Series.

int latAddSocle(LatInfo_t *li, int *mult)
{
    int i;
    int *ptr;

    li->Socle = NREALLOC(li->Socle,int,li->nCf * (li->NSocles + 1));
    ptr = li->Socle + li->nCf * li->NSocles;
    for (i = 0; i < li->nCf; ++i)
	ptr[i] = mult[i];
    ++li->NSocles;
    return li->NSocles;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

/// Add a Layer to the Radical Series.

int latAddHead(LatInfo_t *li, int *mult)
{
    int i;
    int *ptr;

    li->Head = NREALLOC(li->Head,int,li->nCf * (li->NHeads + 1));
    ptr = li->Head + li->nCf * li->NHeads;
    for (i = 0; i < li->nCf; ++i)
	ptr[i] = mult[i];
    ++li->NHeads;
    return li->NHeads;
}

/// @}

// vim:fileencoding=utf8:sw=3:ts=8:et:cin
