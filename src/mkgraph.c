////////////////////////////////////////////////////////////////////////////////////////////////////
// C MeatAxe - Draw the submodule lattice
////////////////////////////////////////////////////////////////////////////////////////////////////

/* Use the PostScript output routine by default */

#include "meataxe.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>




#define LBUFSIZE 2000	/* Input line buffer */
#define MAXIRRED 20	/* Max number of irreducibles */


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


static MtxApplicationInfo_t AppInfo = { 
"mkgraph", "Plot Submodule Lattice",
"\n"
"SYNTAX\n"
"    mkgraph " MTX_COMMON_OPTIONS_SYNTAX " [-c <Colors>] [-b <Block>] "
    "<Name> [<Lower> <Upper>]\n"
"\n"
"OPTIONS\n"
MTX_COMMON_OPTIONS_DESCRIPTION
"    -G ...................... Produce GAP output\n"
"    -b ...................... Select block (Use with mksub -b)\n"
"    -c ...................... Set Colors. Format is `name=R/G/B', where\n"
"                              `name' is any of `std' (standard color),\n"
"                              `line' (lines), `sub' (submodule boxes),\n"
"                              `soc' (socle series), `rad' (radical series),\n"
"                              `mnt' (mountains). R,G,B are integers in the\n"
"                              range 0..99.\n"
"\n"
"FILES\n"
"    <Name>.gra    i  Lattice calculated by mksub\n"
"    <Name>.ps     o  Picture in Postscript format\n"
};

static MtxApplication_t *App = NULL;


static char name[100];
static char ifilename[130];
static char ofilename[130];
static long block = -1;
static Lat_Info LI;		/* Data from .cfinfo file */
static enum { O_PS, O_GAP } OutputMode = O_PS;

struct { char name[20]; long r, b, g; } ColorMap[] =
   {
   	{"std",0,0,0},
   	{"sub",0,0,0},
   	{"rad",0,0,0},
   	{"soc",0,0,0},
   	{"line",0,0,0},
   	{"mnt",0,0,0},
   	{"",0,0,0}	/* End of list marker */
   };


int upper, lower;


/* Data read from the input file
   ----------------------------- */
int nsub;		/* Number of submodules */
long *dim;		/* Dimension of submodules */
int **max;		/* Maximal submodules */
int **maxtype;		/* Types of irreducible factors */
char *issoc;		/* Socle series */
char *israd;		/* Radical series */
char *ismount;		/* Mountains */

/* The factor
   ---------- */
static LdLattice_t *Lattice;



/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(char *msg)

{
    mtxAbort(MTX_HERE,"*** Fatal error: %s",msg);
    exit(1);
}


/* ------------------------------------------------------------------
   readfile() - Read the .gra file
   ------------------------------------------------------------------ */

void readfile(void)

{
    FILE *f;
    char *buf, *c;
    int i, nmax, *lp, *kp;


    /* Open input file
       --------------- */
    f = sysFopen(ifilename,"r");
    if (f == NULL)
    {	perror(ifilename);
		err("Error opening input file");
    }
    buf = (char *)sysMalloc(LBUFSIZE);
    MESSAGE(1,("Reading %s\n",ifilename));

    /* Read number of submodules 
       ------------------------- */
    fgets(buf,LBUFSIZE,f);
    nsub = atoi(buf);
    MESSAGE(1,("%d submodules\n",nsub));

    /* Allocate some arrays
       -------------------- */
    dim = NALLOC(long,nsub);
    max     = NALLOC(int *,nsub);
    maxtype = NALLOC(int *,nsub);
    issoc = NALLOC(char,nsub);
    israd = NALLOC(char,nsub);
    ismount = NALLOC(char,nsub);
    memset(issoc,0,nsub);
    memset(israd,0,nsub);
    memset(ismount,0,nsub);

    /* Read the lattice
       ---------------- */
    for (i = 0; i < nsub; ++i)
    {
	fgets(buf,LBUFSIZE,f);
	for (c = strtok(buf," "); *c != 0; ++c)
	{	
	    if (*c == 'm') ismount[i] = 1; else
	    if (*c == 'r') israd[i] = 1; else
	    if (*c == 's') issoc[i] = 1;
	}
	c = strtok(NULL," ");
	nmax = atoi(c);
	lp = max[i] = NALLOC(int,nmax+1);
	kp = maxtype[i] = NALLOC(int,nmax+1);
	while (nmax-- > 0)
	{
	    *lp++ = atoi(strtok(NULL," "));
	    *kp++ = atoi(strtok(NULL," "));
	}
	*kp = *lp = -1;
    }
    fclose(f);
    free(buf);
/*
printf("max[3] = %d %d %d\n",max[3][0],max[3][1],max[3][2]);
printf("typ[3] = %d %d %d\n",maxtype[3][0],maxtype[3][1],maxtype[3][2]);
*/
}


/* ------------------------------------------------------------------
   SocLevel(), RadLevel() - Returns the level of a given submudule
   in the socle or radical series. Both functions return -1 if the
   submodule does not belong to the series.
   ------------------------------------------------------------------ */

static int SocLevel(int n)

{
    int i, lev=0;
    if (!issoc[n]) return -1;
    for (i = 0; i <= n; ++i)
        if (issoc[i]) ++lev;
    return lev;
}


static int RadLevel(int n)

{
    int i, lev=0;
    if (!israd[n]) return -1;
    for (i = nsub-1; i >= n; --i)
        if (israd[i]) ++lev;
    return lev;
}


/* ------------------------------------------------------------------
   buildroot()
 
   Input: nsub, max, lower, upper
   Output: Lattice
   ------------------------------------------------------------------ */

void buildroot(void)

{   char *flag = NALLOC(char,nsub);
    int *map = NALLOC(int,nsub);
    int finis;
    int i, k;
    int *ip;
    int xnsub;

    if (lower == -1) lower = 0;
    if (upper == -1) upper = nsub-1;

    /* Select all modules below <upper>
       -------------------------------- */
    memset(flag,0,nsub);
    flag[upper] = 1;
    for (finis = 0; !finis; )
    {	finis = 1;
	for (i = 0; i < nsub; ++i)
	{   if (flag[i] == 1)
	    {	for (ip = max[i]; *ip >= 0; ++ip)
		{   if (flag[*ip] == 0)
		    {	flag[*ip] = 1;
			finis = 0;
		    }
		}
		flag[i] = 2;
	    }
	}
    }

    /* Select all modules which are also above <lower>
       ----------------------------------------------- */
    if (flag[lower] != 2) err("Illegal limits\n");
    flag[lower] = 3;
    map[lower] = 0;
    xnsub = 1;
    for (finis = 0; !finis; )
    {	finis = 1;
	for (i = 0; i < nsub; ++i)
	{   if (flag[i] == 2)
	    {	for (ip = max[i]; *ip >= 0; ++ip)
		{   if (flag[*ip] == 3)
		    {	finis = 0;
			flag[i] = 3;
			map[i] = xnsub++;
			break;
		    }
		}
	    }
	}
    }
    if (lower > 0 || upper < nsub - 1)
	MESSAGE(1,("%d modules between %d and %d\n",xnsub,lower,upper));

    /* Calculate the factor lattice
       ---------------------------- */
    Lattice = ldAlloc(xnsub);
    k = 0;
    for (i = 0; i < nsub; ++i)
    {   
	if (flag[i] == 3)
        {   
	    Lattice->Nodes[k].UserData = i;
	    for (ip = max[i]; *ip >= 0; ++ip)
	    {   
		if (flag[*ip] == 3)
		    ldAddIncidence(Lattice,map[*ip],k);
	    }
	    ++k;
	}
    }
}





/* ------------------------------------------------------------------
   setColors() - Set colors from command line option.
   ------------------------------------------------------------------ */

static void setColors(const char *opt_text_ptr)

{
    int i;
    int r,b,g;
    char *c;

    while (*opt_text_ptr != 0)
    {
    	for (i = 0; *(c = ColorMap[i].name) != 0; ++i)
	    if (!strncmp(c,opt_text_ptr,strlen(c))) break;
	if (*c == 0)
	    mtxAbort(MTX_HERE,"-c: %s",MTX_ERR_OPTION);
	opt_text_ptr += strlen(c);
	if (*opt_text_ptr != '=' && *opt_text_ptr != ':')
	    mtxAbort(MTX_HERE,"-c: %s",MTX_ERR_OPTION);
	++opt_text_ptr;
	if (sscanf(opt_text_ptr,"%d/%d/%d",&r,&g,&b) !=3)
	    mtxAbort(MTX_HERE,"-c: %s",MTX_ERR_OPTION);
	while (*opt_text_ptr != 0 && *opt_text_ptr != ',') 
	    ++opt_text_ptr;
	if (*opt_text_ptr == ',') 
	    ++opt_text_ptr;
	if (r < 0 || g < 0 || b < 0 || r > 99 || g > 99 || b > 99) 
	    mtxAbort(MTX_HERE,"color value (-c): %s",MTX_ERR_RANGE);
	ColorMap[i].r = r;
	ColorMap[i].g = g;
	ColorMap[i].b = b;
	MESSAGE(2,("setColor(%s = %d/%d/%d)\n",ColorMap[i].name,
		r,g,b));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static void init(int argc, char **argv)
{	
    const char *c;

    App = appAlloc(&AppInfo,argc,argv);
    block = appGetIntOption(App,"-b",-1,0,-1);
    if (appGetOption(App,"-G"))
	OutputMode = O_GAP;
    if ((c = appGetTextOption(App,"-c",NULL)) != NULL)
	setColors(c);
    if (OutputMode == O_GAP)
	MtxMessageLevel = -1000;
    upper = lower = -1;
    appGetArguments(App,1,3);
    switch (App->ArgC)
    {	case 3:
		upper = atoi(App->ArgV[2]);
	case 2:
		lower = atoi(App->ArgV[1]);
	case 1:
		strcpy(name,App->ArgV[0]);
		break;
    }

    latReadInfo(&LI,name);
    if (block > 0)
    {
	snprintf(ifilename, sizeof(ifilename), "%s.gra.%ld",name,block);
	snprintf(ofilename, sizeof(ofilename), "%s.ps.%ld",name,block);
    }
    else
    {
	snprintf(ifilename, sizeof(ifilename), "%s.gra",name);
	snprintf(ofilename, sizeof(ofilename), "%s.ps",name);
    }
}




/* ------------------------------------------------------------------
   display() Postscript version
   ------------------------------------------------------------------ */

#define XMAP(x) ((x)*xsize + 10)
#define YMAP(y) ((y)*ysize)

FILE *psfile;
double xsize = 18.0 / 2.54 * 72.0;
double ysize = 26.0 / 2.54 * 72.0;
double xbox = 0.6 / 2.54 * 72.0;
double ybox = 0.6 / 2.54 * 72.0;

static char hdr1[] =
    "%%!PS-Adobe-2.0\n"
    "%%%%Creator: mkgraph (%s)\n"
    "%%%%Title: %s\n"
    "%%%%Pages: 1 1\n"
    "%%%%sndComments\n";
static char hdr2[] =
    "/NCols 1 def\n"
    "/NRows 1 def\n"
    "/ThisRow 1 def\n"
    "/ThisCol 1 def\n"
    "/Pagewidth %1.1f def\n"
    "/Pageheight %1.1f def\n"
    "/LeftClip Pagewidth NCols div ThisCol 1 sub mul def\n"
    "/BotClip Pageheight NRows div ThisRow 1 sub mul def\n"
    "NCols NRows scale\n"
    "LeftClip neg BotClip neg translate\n"
    "25 NCols div 25 NRows div translate\n"
    ;
    /* Clippath:
    "newpath LeftClip 50 NCols div add BotClip 50 NRows div add\n"
    "moveto Pagewidth NCols div 0 rlineto\n"
    "0 Pageheight NRows div rlineto Pagewidth NCols div neg 0 rlineto\n"
    "closepath clip\n";
    "[10 5] 0 setdash 0.7 setlinewidth clippath stroke\n";
    */

static char FontName[] = "Helvetica";
static char *FontAlias[] = {"Small","Norm","Big"};
static int FontSize[] = {5,8,12,0};

void writeheader(void)

{
    char mname[100], *c;
    int i;

    fprintf(psfile,hdr1,"ver0.0",ofilename);
    fprintf(psfile,hdr2,xsize,ysize);

    strcpy(mname,ofilename);
    for (c = mname; *c; ++c);
    c[-3] = 0;
    for (i = 0; FontSize[i] > 0; ++i)
    {
    	fprintf(psfile,"/%sFont { /%s findfont %d scalefont setfont } "
	     "def\n",FontAlias[i],FontName,FontSize[i]);
    }
    fprintf(psfile,"BigFont\n");
    fprintf(psfile,"%.1f %.1f moveto (", XMAP(0.0),YMAP(1.0));

    fprintf(psfile,"Module: %s",mname);
    if (block > 0) fprintf(psfile,", Block: %ld",block);
    if (lower != 0 || upper != nsub-1)
	fprintf(psfile,", Range: %d-%d",lower,upper);
    fprintf(psfile,") show \n");

    fprintf(psfile,"NormFont\n");

    fprintf(psfile,"/U { 0 %.1f rlineto } def\n",ybox);
    fprintf(psfile,"/D { 0 -%.1f rlineto } def\n",ybox);
    fprintf(psfile,"/L { -%.1f 0 rlineto } def\n",xbox);
    fprintf(psfile,"/R { %.1f 0 rlineto } def\n",xbox);
    fprintf(psfile,"/UR { %.1f %1.1f rlineto } def\n",xbox/2,ybox/2);
    fprintf(psfile,"/DR { %.1f -%.1f rlineto } def\n",xbox/2,ybox/2);
    fprintf(psfile,"/UL { -%.1f %.1f rlineto } def\n",xbox/2,ybox/2);
    fprintf(psfile,"/DL { -%.1f -%.1f rlineto } def\n",xbox/2,ybox/2);

    /* Definition of Sq, Di, Ci, and Lbl
       --------------------------------- */
    fprintf(psfile,
	"/Sq { subColor 2 copy newpath moveto -%.1f -%.1f rmoveto\n"
	"      U R D L closepath stroke } def\n",
	xbox/2,ybox/2);

    fprintf(psfile,
	"/Di { radColor 2 copy newpath moveto 0 -%1.1f rmoveto\n"
	"      UR UL DL DR closepath stroke } def\n",
	ybox/2);

    fprintf(psfile,
	"/Ci { socColor 2 copy newpath %1.1f 0 360 arc stroke } def\n",
     	ybox/2);

    fprintf(psfile,
	"/Lbl { stdColor newpath NormFont Thin moveto dup stringwidth pop\n"
        "       2 div neg -3 rmoveto show stroke } def\n");

    fprintf(psfile,
	"/RadLbl { stdColor newpath SmallFont Thin moveto -%1.1f %1.1f "
	"rmoveto dup stringwidth pop 2 add neg -3 rmoveto show stroke "
	"} def\n",
	xbox/2,ybox/2);
    fprintf(psfile,
	"/SocLbl { stdColor newpath SmallFont Thin moveto %1.1f 2 add "
	"-%1.1f rmoveto show stroke } def\n",
	xbox/2,ybox/2);

    /* Definition of Thin and Thick
       --------------------------- */
    fprintf(psfile,"/Thin { %1.1f setlinewidth } def\n",0.4);
    fprintf(psfile,"/Thick { %1.1f setlinewidth } def Thin\n",1.2);
    
    /* Definition of Colors
       --------------------------- */
    for (i = 0; *(c = ColorMap[i].name) != 0; ++i)
    	fprintf(psfile,"/%sColor {0.%2.2ld 0.%2.2ld 0.%2.2ld "
	    "setrgbcolor} def\n",
    	    ColorMap[i].name, ColorMap[i].r, ColorMap[i].g,
	    ColorMap[i].b);
}



void shownode(int i,double x,double y)

{   
    if (ismount[i]) fprintf(psfile,"Thick mntColor ");
    fprintf(psfile,"(%d) %1.1f %1.1f ",i,XMAP(x),YMAP(y));
    if (israd[i])
    {
    	fprintf(psfile,"Di ");
        fprintf(psfile,"(%d) %1.1f %1.1f RadLbl ",RadLevel(i),
		XMAP(x),YMAP(y));
    }
    if (issoc[i])
    {
    	fprintf(psfile,"Ci ");
        fprintf(psfile,"(%d) %1.1f %1.1f SocLbl ",SocLevel(i),
		XMAP(x),YMAP(y));
    }
    if (!issoc[i] && !israd[i]) fprintf(psfile,"Sq ");
    fprintf(psfile,"Lbl\n");
}


static char *linestyle[MAXIRRED] = {
"[] 0",
"[1 1] 0",
"[3 3] 0",
"[3 1 1 1] 0",
"[1 1 3 1 1 1] 0",
"[3 1 1 1 3 1] 0",
"[1 1 1 1 3 1 1 1] 0",
"[1 1 3 1 3 1 1 1] 0",
"[3 1 1 1 3 1 3 1] 0",
"[5 1] 0",
"[5 1 1 1] 0",
"[5 1 3 1] 0",
"[5 1 1 1 1 1] 0",
"[5 1 3 1 3 1] 0",
"[5 1 3 1 1 1] 0",
"[5 1 5 1 1 1] 0",
"[5 1 5 1 3 1] 0",
"[5 1 1 1 1 1 1 1] 0",
"[5 1 1 1 3 1 1 1] 0",
"[5 1 3 1 3 1 3 1] 0",
};


void showline(int i, int k,int type)

{	
    int t = type;
    if (t < 0) t = 0;
    if (type >= MAXIRRED) t = MAXIRRED-1;
    fprintf(psfile,"lineColor newpath %s setdash %% type=%d\n",
    	linestyle[t],type);
    fprintf(psfile,"%1.1f %1.1f moveto ",
	  XMAP(Lattice->Nodes[i].PosX),YMAP(Lattice->Nodes[i].PosY)+ybox/2);
    fprintf(psfile,"%1.1f %1.1f lineto\n",
	  XMAP(Lattice->Nodes[k].PosX),YMAP(Lattice->Nodes[k].PosY)-ybox/2);
    fprintf(psfile,"stroke [] 0 setdash\n");
}



void writelegend()

{
    int i;

    fprintf(psfile,"%% Legend\n%% -------\nnewpath\n");
    for (i = 0; i < LI.NCf; ++i)
    {
	fprintf(psfile,
	    "lineColor %s setdash %1.1f %1.1f moveto 60 0 rlineto stroke\n",
	    linestyle[i],XMAP(0.8),YMAP(1.0)-10.0*i);
	fprintf(psfile,
	    "stdColor [] 0 setdash %1.1f %1.1f moveto ",
	    XMAP(0.8)+65.0,YMAP(1.0)-10.0*i-3.0);
	fprintf(psfile,"(%s) show stroke\n",latCfName(&LI,i));
    }
    fprintf(psfile,"\n");
}


void display()

{
    int i, l;

    MESSAGE(0,("Writing lattice diagram to %s\n",ofilename));
    fflush(stdout);
    psfile = sysFopen(ofilename,"w");
    if (psfile == NULL)
    {	perror(ofilename);
	exit(1);
    }
    writeheader();
    writelegend();

    for (i = 0; i < Lattice->NNodes; ++i)
    {	
	fprintf(psfile,"1 { ");
	shownode(Lattice->Nodes[i].UserData,Lattice->Nodes[i].PosX,
	    Lattice->Nodes[i].PosY);
	fprintf(psfile,"newpath\n");
	for (l = 0; l < Lattice->NNodes; ++l)
	{
	    if (LD_ISSUB(Lattice,l,i))
	    {
		int ni = Lattice->Nodes[i].UserData;
		int nl = Lattice->Nodes[l].UserData;
		int m;
		for (m = 0; max[i][m] >= 0 && max[i][m] != nl; ++m);
		showline(l,i,maxtype[ni][m]);
if ((l == 3 && i < 3) || (l < 3 && i == 3))
{
printf("Show line %d-%d: orig=%d-%d, m=%d, type=%d\n",
    l,i,nl,ni,m,maxtype[ni][m]);
}
	    }
	}

	fprintf(psfile,"} repeat\n");
    }
    fprintf(psfile,"showpage\n");
    fprintf(psfile,"%%%%sOF\n");
    fclose(psfile);
}




/* ------------------------------------------------------------------
   DisplayGap() - Generate GAP output
   ------------------------------------------------------------------ */

const int GapXSize = 800;
const int GapYSize = 600;
const char *GapLatName = "MtxLattice";
const char *GapVlName = "MtxVertexList";

static void DisplayGap()

{
    int i, l;

    printf("# Generated by mkgraph\n");

    /* Generate Poset und Levels
       ------------------------- */
    printf("%s := GraphicMeatAxeLattice(\"%s\",%d,%d);\n",GapLatName,name,
	GapXSize,GapYSize);
    for (l = 0; l < Lattice->NLayers; ++l)
	printf("CreateLevel(%s,%d);\n",GapLatName,l);

    /* Insert vertices into the Poset
       ------------------------------ */
    printf("%s := [];\n",GapVlName);
    for (i = 0; i < Lattice->NNodes; ++i)
    {
	LdNode_t *n = Lattice->Nodes + i;
	char tmp[50];
	snprintf(tmp, sizeof(tmp), "%d",i);
	printf("Add(%s,Vertex(%s,rec(SubmoduleNumber:=%d),"
	    "rec(x:=%d,levelparam:=%d,label:=\"%s\",shape:=\"%s\")));\n",
	    GapVlName,GapLatName,i,(int)(n->PosX * GapXSize),n->Layer,
	    tmp,ismount[Lattice->Nodes[i].UserData] ? "diamond" : "circle");
    }

    /* Insert edges into the Poset
       --------------------------- */
    for (i = 0; i < Lattice->NNodes; ++i)
    {
	int k;
    	for (k = 0; k < Lattice->NNodes; ++k)
	{
	    if (!LD_ISSUB(Lattice,i,k))
		continue;
	    printf("Edge(%s,%s[%d],%s[%d]);\n",
		GapLatName,GapVlName,i+1,GapVlName,k+1);
	}
    }
/*    printf("delete(%s);\n",GapVlName);*/

/*
CreateLevel(p, level-nr, "level-label");

Fuer alle Vertices:
vertexliste := [];
inf := rec(x := X-Position [pixel], levelparam:=level-nr,
	color:=COLORS.red, shape:="circle/diamond/rectangle");
Add(vertexliste,Vertex(p,rec(flags...), inf));

Fuer alle Kanten:
Edge(p,vertexliste[from[i]],vertexliste[to[i]]);
*/
}


////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    init(argc,argv);
    readfile();
    buildroot();
    ldSetPositions(Lattice);
    switch (OutputMode)
    {
	case O_GAP: DisplayGap(); break;
	default: display(); break;
    }
    return 0;
}


/**
@page prog_mkgraph mkgraph - Draw a Submodule Lattice

@section mkgraph_syntax Command Line
<pre>
mkgraph [@em Options] [-G] [-b @em BlockNo] @em Name [@em Lower @em Upper]A
</pre>

@par @em Options
  Standard options, see @ref prog_stdopts
@par -G
  Generate (X)GAP output instead of postscript. Output is written to stdout.
@par -b
  Select block number @em BlockNo for drawing. Must be used if @ref prog_mksub "mksub"
  has been run in "block mode".
@par @em Name
  Name of the representation.
@par @em Lower
  Numer of the smallest submodule to be drawn.
@par @em Upper
  Numer of the greatest submodule to be drawn.

@section mkgraph_inp Input Files
@par @em Name.gra
  Lattice information.

@section mkgraph_out Output Files
@par @em Name.ps
  Lattice diagram.

@section mkgraph_desc Description
This program creates a graphical representation of a submodule lattice
in postscript or other formats. The first argument must be the module
name. @ref prog_mksub "mksub" must have been run, because @b mkgraph reads
the ".gra" file created by @ref prog_mksub "mksub".
If two additional arguments are specified, only the modules between 
@em Lower and @em Upper are drawn. If @ref prog_mksub "mksub" was used in "block mode", 
only single blocks can be drawn. In this case, the block number must 
be specified with "-b".

If no other format is specified, @b mkgraph produces PosScript output.
Submodules are represented by boxes, and each Submodule is linked with its maximal
submodules by a line. Different line styles are used to distinguish between irreducible
constituents.
Each node carries the submodule number as found in the @em Name.out file.
Local submodules are displayed with thick boxes. The socles are shown
as circles instead of boxes, and the radicals are shown in diamond shapes. 
The line style represents the irreducible constituent isomorphic to the 
factor module represented by the line.
The output is written to @em Name.ps unless "-G" is used.


If the option "-G" is used, @b mkgraph creates commands that can be read by xGAP.

@section mkgraph_impl Implementation Details
The algorithm used to position the nodes (submodules) is very simple. In
most cases the result is far from optimal, but also much better than a
random drawing. Submodules are grouped into layers according to their
composition length. All submodules in one layer are drawn at the same
y (vertical) position, and at equidistant x positions. The program tries
to optimize the order of submodules in each layer, such that lines between
the submodules become short.
**/
// vim:fileencoding=utf8:sw=3:ts=8:et:cin
