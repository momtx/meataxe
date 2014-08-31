#include <string.h>
#include <stdio.h>

static int IsBigEndian()
{
    union {
        unsigned char c[sizeof(int)];
	unsigned u;
    } x;
    x.u = 1;
    return x.c[0] == 0;
}


int main()
{
    FILE *f = 0;
    static char svnversion[50] = "UNKNOWN";

    printf("#ifndef CONFIG_H_INCLUDED\n");
    printf("#define CONFIG_H_INCLUDED\n");
    printf("#define MTX_CONFIG_BIG_ENDIAN %d\n",IsBigEndian());
    printf("#define MTX_CONFIG_LONG32 %d\n",sizeof(long) == 4);
    printf("#define MTX_CONFIG_LONG64 %d\n",sizeof(long) == 8);

    printf("#define MTX_CONFIG_STRING \"");
    if (sizeof(long) == 8) printf(" 64");
    if (sizeof(long) == 4) printf(" 32");
    if (IsBigEndian()) printf(" BIG_ENDIAN");
    printf("\"\n");

    if ((f = fopen("svnversion","r")) != 0) {
	char line[sizeof(svnversion) + 1];
	if (fgets(line,sizeof(line),f) != 0) {
	    char *c = line + strlen(line);
	    while (c > line && (unsigned char)c[-1] <= 32) --c;
	    if (c > line) {
	        *c = 0;
		for (c = line; *c; ++c) {
		    if (*c == '"' || *c == '\\')
			*c = '_';
		}
		strcpy(svnversion,line);
	    }
	}
	fclose(f);
    }
    printf("#define MTX_BUILD \"%s\"\n",svnversion);

    printf("#endif\n");
    return 0;
}

