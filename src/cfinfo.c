////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Functions for reading and writing the .cfinfo file
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "meataxe.h"

#include <string.h>
#include <stdlib.h>
#include <ctype.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Local data



/// @defgroup cfinfo Module Structure
/// @{

/// @class CfInfo
/// Constituent data structure.
/// The CFInfo data structure contains all information about one 
/// irreducible constituent of a module. Most of this information is
/// useful only in the context of a Lat_Info structure.
///
/// @c dim is the dimension of the constituent, and @c num enumerates
/// constituents of same dimension. @c mult ist the multiplicity of this
/// constituent in the module.
/// @c idword and @c idpol specify the identifying word for this module.
/// The identifying word is used by the CHOP program.

/// @class Lat_Info
/// Module Structure Information.
/// This data structure contains all information about a module which is
/// gathered during the submodule lattice calculation.

/* ------------------------------------------------------------------
   writeword(), readword() - Write and read words
   ------------------------------------------------------------------ */

static void WriteWord(StfData *f, long w, Poly_t *p)
{
    int i;

    stfPut(f,"[");
    stfPutInt(f,w);
    stfPut(f,",");
    stfPutInt(f,ffOrder);
    if (p == NULL)
    {
	stfPut(f,",-1");
    }
    else
    {
	stfPut(f,",");
	stfPutInt(f,p->Degree);
	for (i = 0; i <= p->Degree; ++i)
	{
	    stfPut(f,",");
	    stfPutInt(f,ffToInt(p->Data[i]));
	}
    }
    stfPut(f,"]");
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reads a word number and the associated polynomial.
/// @return 0 on success, -1 on error

static int ReadWord(StfData *f, long *w, Poly_t **p, const char *fn)
{
    int fl, deg;  
    int i;

    if (stfMatch(f," [")) {
	mtxAbort(MTX_HERE,"%s: missing '['",fn);
	return -1;
    }
    stfGetInt(f,&i);
    *w = i;
    if (stfMatch(f,",")) {
	mtxAbort(MTX_HERE,"%s: missing ','",fn);
	return -1;
    }
    stfGetInt(f,&fl);
    if (stfMatch(f,",")) {
	mtxAbort(MTX_HERE,"%s: missing ','",fn);
	return -1;
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
		return -1;
	    }
	    stfGetInt(f,&coeff);
	    (*p)->Data[i] = ffFromInt(coeff);
    	}
    }
    if (stfMatch(f,"]")) {
	mtxAbort(MTX_HERE,"%s: missing ']'",fn);
	return -1;
    }
    return 0;
}


#define RDVEC(name,fld)\
	else if (!strcmp(c,name))\
	{\
	    if (stfMatch(f," ["))\
	    { mtxAbort(MTX_HERE,"Error in cfinfo file: Missing '['"); return -1; }\
	    for (int i = 0; i < li->NCf; ++i) \
	    {\
		int val = 0;\
		if (i > 0)\
		    stfMatch(f,",");\
		stfGetInt(f,&val);\
		li->Cf[i].fld = val;\
	    }\
	    if (stfMatch(f,"]"))\
	    { mtxAbort(MTX_HERE,"Error in cfinfo file: Missing ']'"); return -1; }\
	}


static int readCfFile(StfData* f, const char* fn, Lat_Info* li)
{

    /* Read header
       ----------- */
    if (stfReadLine(f) || strcmp(stfGetName(f),"CFInfo"))
    {
	mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
	return -1;
    }

    /* Read data
       --------- */
    while (stfReadLine(f) == 0)
    {
	const char *c = stfGetName(f);
	if (!strcmp(c,"CFInfo.NCF"))
	    stfGetInt(f,&li->NCf);
	else if (!strcmp(c,"CFInfo.ConstituentNames"))
	    ;	/* Ignore */
	else if (!strcmp(c,"CFInfo.Field"))
	{   
	    stfGetInt(f,&li->Field);
	    ffSetField(li->Field);
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
	    {
		mtxAbort(MTX_HERE,"%s: Missing '['",fn);
		return -1;
	    }
	    for (int i = 0; i < li->NCf; ++i)
	    {
		if (ReadWord(f,&(li->Cf[i].idword),&(li->Cf[i].idpol),fn) != 0) {
		    return -1;
		}
		if (stfMatch(f,i < li->NCf - 1 ? "," : "];") != 0)
		{
		    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		    return -1;
		}
	    }
	}
	else if (!strcmp(c,"CFInfo.PeakWord"))
	{
	    if (stfMatch(f," [") != 0) 
	    {
		mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		return -1;
	    }
	    for (int i = 0; i < li->NCf; ++i)
	    {
		if (ReadWord(f,&(li->Cf[i].peakword),&(li->Cf[i].peakpol),fn) != 0) {
		    return -1;
		}
		if (stfMatch(f,i < li->NCf - 1 ? "," : "];") != 0)
		{
		    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		    return -1;
		}
	    }
	}
	else if (!strcmp(c,"CFInfo.LoewyLength"))   /* for compatibility */
	    ;	
	else if (!strcmp(c,"CFInfo.NSocles"))
	    ;	/* Is set when reading "CFInfo.Socles" */
	else if (!strcmp(c,"CFInfo.Socles"))
	{
	    if (stfMatch(f," [") != 0) 
	    {
		mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		return -1;
	    }
	    for (int i = 0; stfMatch(f,"];"); ++i)
	    {
		int mult[LAT_MAXCF];
		int count = LAT_MAXCF;
		if (i > 0) stfMatch(f,",");
		stfGetVector(f,&count,mult);
		if (count != li->NCf)
		{
		    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		    return -1;
		}
		latAddSocle(li,mult);
	    }
	}
	else if (!strcmp(c,"CFInfo.NHeads"))
	    ;	/* Is set when reading "CFInfo.Heads" */
	else if (!strcmp(c,"CFInfo.Heads"))
	{
	    if (stfMatch(f," [") != 0) 
	    {
		mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		return -1;
	    }
	    for (int i = 0; stfMatch(f,"];"); ++i)
	    {
		int mult[LAT_MAXCF];
		int count = LAT_MAXCF;
		if (i > 0) stfMatch(f,",");
		stfGetVector(f,&count,mult);
		if (count != li->NCf)
		{
		    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
		    return -1;
		}
		latAddHead(li,mult);
	    }
	}
	else
	{
	    mtxAbort(MTX_HERE,"%s: %s",fn,MTX_ERR_FILEFMT);
	    return -1;
	}
    }
    return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reads a Lattice Information File.
/// This funktion reads a .cfinfo file and stores the data into a Lat_Info structure.
/// @param li Pointer to the data structure.
/// @param basename Base name (without the trailing ".cfinfo").
/// @return 0 on success, -1 on error.

int latReadInfo(Lat_Info *li, const char *basename)
{
    MTX_ASSERT(li != NULL, -1);
    MTX_ASSERT(basename != NULL, -1);
    MTX_ASSERT(strlen(basename) < LAT_MAXBASENAME - 1, -1);

    // Initialize the data structure.
    memset(li,0,sizeof(Lat_Info));
    strcpy(li->BaseName,basename);
    li->NGen = 2;

    // Open the file and read
    char fn[LAT_MAXBASENAME + 20];
    sprintf(fn,"%s.cfinfo",basename);
    int result = 0;
    StfData *f = stfOpen(fn,"r");
    if (f == NULL) {
	mtxAbort(MTX_HERE,"Cannot open %s",fn);
	result = -1;
    } else {
       result = readCfFile(f, fn, li);
    }
    stfClose(f);
    if (result == 0)
       mtxMessage(1,"Read %s: %d composition factors",fn,li->NCf);
    return result;
}



/// Write a Lattice Information File.
/// This funktion writes the contents of a Lat_Info structure into a file.
/// The file name is constructed from the @c BaseName field of the structure
/// by appending ".cfinfo".
/// @param li Pointer to the data structure.
/// @return 0 on success, -1 on error.

int latWriteInfo(const Lat_Info *li)
{
    StfData *f;
    int i;
    int tmp[LAT_MAXCF];
    char fn[LAT_MAXBASENAME + 20];

    MTX_ASSERT(li != NULL, -1);

    /* Open the file
       ------------- */
    strcpy(fn,li->BaseName);
    strcat(fn,".cfinfo");
    f = stfOpen(fn,"w");
    if (f == NULL)
	return -1;

    /* Write data
       ---------- */
    stfWriteValue(f,"CFInfo","rec()");
    stfWriteInt(f,"CFInfo.NGen",li->NGen);
    stfWriteInt(f,"CFInfo.Field",li->Field);
    stfWriteInt(f,"CFInfo.NCF",li->NCf);


    stfBeginEntry(f,"CFInfo.ConstituentNames");
    stfPut(f,"[");
    for (i = 0; i < li->NCf; ++i)
    {
	stfPutString(f,latCfName(li,i));
	if (i < li->NCf-1) stfPut(f,",");
    }
    stfPut(f,"]");
    stfEndEntry(f);

#define WRVEC(name,field)\
    for (i = 0; i < li->NCf; ++i)\
        tmp[i] = li->Cf[i].field;\
    stfWriteVector(f,"CFInfo." #name,li->NCf,tmp);

    WRVEC(Dimension,dim);
    WRVEC(Number,num);
    WRVEC(Multiplicity,mult);
    WRVEC(SplittingField,spl);
    WRVEC(NMountains,nmount);
    WRVEC(NDottedLines,ndotl);

    stfBeginEntry(f,"CFInfo.PeakWord");
    stfPut(f,"[");
    for (i = 0; i < li->NCf; ++i)
    {
        WriteWord(f,li->Cf[i].peakword,li->Cf[i].peakpol);
	if (i < li->NCf-1) stfPut(f,",");
    }
    stfPut(f,"]");
    stfEndEntry(f);

    stfBeginEntry(f,"CFInfo.IdWord");
    stfPut(f,"[");
    for (i = 0; i < li->NCf; ++i)
    {
        WriteWord(f,li->Cf[i].idword,li->Cf[i].idpol);
	if (i < li->NCf-1) stfPut(f,",");
    }
    stfPut(f,"]");
    stfEndEntry(f);

    stfWriteInt(f,"CFInfo.NSocles",li->NSocles);
    stfBeginEntry(f,"CFInfo.Socles");
    stfPut(f,"[");
    for (i = 0; i < li->NSocles; ++i)
    {
	if (i > 0) stfPut(f,",");
	stfPutVector(f,li->NCf,li->Socle + i * li->NCf);
    }
    stfPut(f,"]");
    stfEndEntry(f);

    stfWriteInt(f,"CFInfo.NHeads",li->NHeads);
    stfBeginEntry(f,"CFInfo.Heads");
    stfPut(f,"[");
    for (i = 0; i < li->NHeads; ++i)
    {
	if (i > 0) stfPut(f,",");
	stfPutVector(f,li->NCf,li->Head + i * li->NCf);
    }
    stfPut(f,"]");
    stfEndEntry(f);

    stfClose(f);
    MESSAGE(1,("Wrote %s: %d composition factors\n",fn,li->NCf));
    return 0;
}




/// Make Constituent Name.
/// This function returns the name of the @a cf-th constituent of a module.
/// The constituent name consists of the dimension and an appendix which is
/// built from the @c num field in the constituent's data structure. Usually
/// the appendix is a single letter ('a', 'b', ...). If there are more than 
/// 26 constituents with the same dimension, a two-letter appendix (`aa', 
/// `ab', etc.) is used. 
///
/// Note: The return value points to a static buffer which is overwritten 
/// at each call. The constituent data inside @a li must have been set up 
/// properly, i.e., the module must have been chopped.
/// @param li Constituent info structure.
/// @param cf Index of the constituent.
/// @return Pointer to the constituent name (without base name).

const char *latCfName(const Lat_Info *li, int cf)
{
    static char buf[20];
    int num, dim;

    buf[0] = 0;
    MTX_ASSERT(li != NULL, buf);
    MTX_ASSERT(cf >= 0 && cf < li->NCf, buf);

    /* Get dimension and number of the constituent
       ------------------------------------------- */
    num = li->Cf[cf].num;
    dim = li->Cf[cf].dim;

    /* Build the constituent's name
       ---------------------------- */
    if (num < 26)
	sprintf(buf,"%d%c",dim,(char)num+'a');
    else if (num < 26*26)
	sprintf(buf,"%d%c%c",dim,(char)(num/26-1)+'a',(char)(num%26)+'a');
    else
	sprintf(buf,"%dcf%d",dim,num);

    return buf;
}


/// Add a Layer to the Socle Series.

int latAddSocle(Lat_Info *li, int *mult)
{
    int i;
    int *ptr;

    li->Socle = NREALLOC(li->Socle,int,li->NCf * (li->NSocles + 1));
    ptr = li->Socle + li->NCf * li->NSocles;
    for (i = 0; i < li->NCf; ++i)
	ptr[i] = mult[i];
    ++li->NSocles;
    return li->NSocles;
}


/// Add a Layer to the Radical Series.

int latAddHead(Lat_Info *li, int *mult)
{
    int i;
    int *ptr;

    li->Head = NREALLOC(li->Head,int,li->NCf * (li->NHeads + 1));
    ptr = li->Head + li->NCf * li->NHeads;
    for (i = 0; i < li->NCf; ++i)
	ptr[i] = mult[i];
    ++li->NHeads;
    return li->NHeads;
}

/// @}
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
