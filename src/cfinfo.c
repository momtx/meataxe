////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Functions for reading and writing the .cfinfo file
//
// (C) Copyright 1998-2014 Michael Ringe, Lehrstuhl D fuer Mathematik, RWTH Aachen
//
// This program is free software; see the file COPYING for details.
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <meataxe.h>

#include <string.h>
#include <stdlib.h>
#include <ctype.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Local data

MTX_DEFINE_FILE_INFO


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

    StfPut(f,"[");
    StfPutInt(f,w);
    StfPut(f,",");
    StfPutInt(f,FfOrder);
    if (p == NULL)
    {
	StfPut(f,",-1");
    }
    else
    {
	StfPut(f,",");
	StfPutInt(f,p->Degree);
	for (i = 0; i <= p->Degree; ++i)
	{
	    StfPut(f,",");
	    StfPutInt(f,FfToInt(p->Data[i]));
	}
    }
    StfPut(f,"]");
}


static int ReadWord(StfData *f, long *w, Poly_t **p, char *fn)
{
    int fl, deg;  
    int i;

    if (StfMatch(f," ["))
    {
	MTX_ERROR1("%s: missing '['",fn);
	return 0;
    }
    StfGetInt(f,&i);
    *w = i;
    if (StfMatch(f,",")) 
    {
	MTX_ERROR1("%s: missing ','",fn);
	return 0;
    }
    StfGetInt(f,&fl);
    if (StfMatch(f,",")) 
    {
	MTX_ERROR1("%s: missing ','",fn);
	return 0;
    }
    StfGetInt(f,&deg);
    if (deg == -1)
	*p = NULL;
    else
    {
        *p = PolAlloc(fl,deg);
    	for (i = 0; i <= deg; ++i)
    	{
	    int coeff;
            if (StfMatch(f,",")) 
	    {
		MTX_ERROR1("%s: missing ','",fn);
		return 0;
	    }
	    StfGetInt(f,&coeff);
	    (*p)->Data[i] = FfFromInt(coeff);
    	}
    }
    if (StfMatch(f,"]")) 
    {
	MTX_ERROR1("%s: missing ']'",fn);
	return 0;
    }
    return 1;
}


#define RDVEC(name,fld)\
	else if (!strcmp(c,name))\
	{\
	    if (StfMatch(f," ["))\
	    { MTX_ERROR("Error in cfinfo file: Missing '['"); return -1; }\
	    for (i = 0; i < li->NCf; ++i) \
	    {\
		int val = 0;\
		if (i > 0)\
		    StfMatch(f,",");\
		StfGetInt(f,&val);\
		li->Cf[i].fld = val;\
	    }\
	    if (StfMatch(f,"]"))\
	    { MTX_ERROR("Error in cfinfo file: Missing ']'"); return -1; }\
	}

/// Read a Lattice Information File.
/// This funktion reads a .cfinfo file and stores the data into a Lat_Info structure.
/// @param li Pointer to the data structure.
/// @param basename Base name (without the trailing ".cfinfo").
/// @return 0 on success, -1 on error.

int Lat_ReadInfo(Lat_Info *li, const char *basename)
{
    int i;
    StfData *f;
    char fn[LAT_MAXBASENAME + 20];

    /* Check arguments
       --------------- */
    MTX_VERIFY(li != NULL);
    MTX_VERIFY(basename != NULL);
    MTX_VERIFY(strlen(basename) < LAT_MAXBASENAME - 1);

    /* Initialize the data structure
       ----------------------------- */
    memset(li,0,sizeof(Lat_Info));
    strcpy(li->BaseName,basename);
    li->NGen = 2;

    /* Open the file
       ------------- */
    sprintf(fn,"%s.cfinfo",basename);
    if ((f = StfOpen(fn,FM_READ)) == NULL)
    {
	MTX_ERROR1("Cannot open %s",fn);
	return -1;
    }

    /* Read header
       ----------- */
    if (StfReadLine(f) || strcmp(StfGetName(f),"CFInfo"))
    {
	MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
	return -1;
    }

    /* Read data
       --------- */
    while (StfReadLine(f) == 0)
    {
	const char *c = StfGetName(f);
	if (!strcmp(c,"CFInfo.NCF"))
	    StfGetInt(f,&li->NCf);
	else if (!strcmp(c,"CFInfo.ConstituentNames"))
	    ;	/* Ignore */
	else if (!strcmp(c,"CFInfo.Field"))
	{   
	    StfGetInt(f,&li->Field);
	    FfSetField(li->Field);
	}
	else if (!strcmp(c,"CFInfo.NGen"))
	    StfGetInt(f,&li->NGen);
	RDVEC("CFInfo.Dimension",dim)
	RDVEC("CFInfo.Number",num)
	RDVEC("CFInfo.Multiplicity",mult)
	RDVEC("CFInfo.SplittingField",spl)
	RDVEC("CFInfo.NMountains",nmount)
	RDVEC("CFInfo.NDottedLines",ndotl)
	else if (!strcmp(c,"CFInfo.IdWord"))
	{
	    if (StfMatch(f," [") != 0) 
	    {
		MTX_ERROR1("%s: Missing '['",fn);
		return -1;
	    }
	    for (i = 0; i < li->NCf; ++i)
	    {
		ReadWord(f,&(li->Cf[i].idword),&(li->Cf[i].idpol),fn);
		if (StfMatch(f,i < li->NCf - 1 ? "," : "];") != 0)
		{
		    MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
		    return -1;
		}
	    }
	}
	else if (!strcmp(c,"CFInfo.PeakWord"))
	{
	    if (StfMatch(f," [") != 0) 
	    {
		MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
		return -1;
	    }
	    for (i = 0; i < li->NCf; ++i)
	    {
		ReadWord(f,&(li->Cf[i].peakword),&(li->Cf[i].peakpol),fn);
		if (StfMatch(f,i < li->NCf - 1 ? "," : "];") != 0)
		{
		    MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
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
	    if (StfMatch(f," [") != 0) 
	    {
		MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
		return -1;
	    }
	    for (i = 0; StfMatch(f,"];"); ++i)
	    {
		int mult[LAT_MAXCF];
		int count = LAT_MAXCF;
		if (i > 0) StfMatch(f,",");
		StfGetVector(f,&count,mult);
		if (count != li->NCf)
		{
		    MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
		    return -1;
		}
		Lat_AddSocle(li,mult);
	    }
	}
	else if (!strcmp(c,"CFInfo.NHeads"))
	    ;	/* Is set when reading "CFInfo.Heads" */
	else if (!strcmp(c,"CFInfo.Heads"))
	{
	    if (StfMatch(f," [") != 0) 
	    {
		MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
		return -1;
	    }
	    for (i = 0; StfMatch(f,"];"); ++i)
	    {
		int mult[LAT_MAXCF];
		int count = LAT_MAXCF;
		if (i > 0) StfMatch(f,",");
		StfGetVector(f,&count,mult);
		if (count != li->NCf)
		{
		    MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
		    return -1;
		}
		Lat_AddHead(li,mult);
	    }
	}
	else
	{
	    MTX_ERROR2("%s: %E",fn,MTX_ERR_FILEFMT);
	    return -1;
	}
    }
    StfClose(f);
    MESSAGE(1,("Read %s: %d composition factors\n",fn,li->NCf));
    return 0;
}



/// Write a Lattice Information File.
/// This funktion writes the contents of a Lat_Info structure into a file.
/// The file name is constructed from the @c BaseName field of the structure
/// by appending ".cfinfo".
/// @param li Pointer to the data structure.
/// @return 0 on success, -1 on error.

int Lat_WriteInfo(const Lat_Info *li)
{
    StfData *f;
    int i;
    int tmp[LAT_MAXCF];
    char fn[LAT_MAXBASENAME + 20];

    /* Check arguments
       --------------- */
    MTX_VERIFY(li != NULL);

    /* Open the file
       ------------- */
    strcpy(fn,li->BaseName);
    strcat(fn,".cfinfo");
    f = StfOpen(fn,FM_CREATE);
    if (f == NULL)
	return -1;

    /* Write data
       ---------- */
    StfWriteValue(f,"CFInfo","rec()");
    StfWriteInt(f,"CFInfo.NGen",li->NGen);
    StfWriteInt(f,"CFInfo.Field",li->Field);
    StfWriteInt(f,"CFInfo.NCF",li->NCf);


    StfBeginEntry(f,"CFInfo.ConstituentNames");
    StfPut(f,"[");
    for (i = 0; i < li->NCf; ++i)
    {
	StfPutString(f,Lat_CfName(li,i));
	if (i < li->NCf-1) StfPut(f,",");
    }
    StfPut(f,"]");
    StfEndEntry(f);

#define WRVEC(name,field)\
    for (i = 0; i < li->NCf; ++i)\
        tmp[i] = li->Cf[i].field;\
    StfWriteVector(f,"CFInfo." #name,li->NCf,tmp);

    WRVEC(Dimension,dim);
    WRVEC(Number,num);
    WRVEC(Multiplicity,mult);
    WRVEC(SplittingField,spl);
    WRVEC(NMountains,nmount);
    WRVEC(NDottedLines,ndotl);

    StfBeginEntry(f,"CFInfo.PeakWord");
    StfPut(f,"[");
    for (i = 0; i < li->NCf; ++i)
    {
        WriteWord(f,li->Cf[i].peakword,li->Cf[i].peakpol);
	if (i < li->NCf-1) StfPut(f,",");
    }
    StfPut(f,"]");
    StfEndEntry(f);

    StfBeginEntry(f,"CFInfo.IdWord");
    StfPut(f,"[");
    for (i = 0; i < li->NCf; ++i)
    {
        WriteWord(f,li->Cf[i].idword,li->Cf[i].idpol);
	if (i < li->NCf-1) StfPut(f,",");
    }
    StfPut(f,"]");
    StfEndEntry(f);

    StfWriteInt(f,"CFInfo.NSocles",li->NSocles);
    StfBeginEntry(f,"CFInfo.Socles");
    StfPut(f,"[");
    for (i = 0; i < li->NSocles; ++i)
    {
	if (i > 0) StfPut(f,",");
	StfPutVector(f,li->NCf,li->Socle + i * li->NCf);
    }
    StfPut(f,"]");
    StfEndEntry(f);

    StfWriteInt(f,"CFInfo.NHeads",li->NHeads);
    StfBeginEntry(f,"CFInfo.Heads");
    StfPut(f,"[");
    for (i = 0; i < li->NHeads; ++i)
    {
	if (i > 0) StfPut(f,",");
	StfPutVector(f,li->NCf,li->Head + i * li->NCf);
    }
    StfPut(f,"]");
    StfEndEntry(f);

    StfClose(f);
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

const char *Lat_CfName(const Lat_Info *li, int cf)
{
    static char buf[20];
    int num, dim;

    /* Check arguments
       --------------- */
    MTX_VERIFY(li != NULL);
    MTX_VERIFY(cf >= 0 && cf < li->NCf);

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

int Lat_AddSocle(Lat_Info *li, int *mult)
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

int Lat_AddHead(Lat_Info *li, int *mult)
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
